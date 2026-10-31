"""Fixed-cell, ordinary-PA6 ordered-branch GPAW resource/scientific pilot.

No DFT is started by import, --help, --check-input, or --dry-run.  --dry-run
initializes the actual MPI calculator and estimates memory, without SCF.
Every run needs a fresh output directory and a matching independent geometry
validation report. Numerical parameters here are NOT PA6 production settings.
"""
from __future__ import annotations

import argparse
from collections import Counter
import fcntl
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import sys
import time
import traceback

import numpy as np
import spglib
from ase.calculators.mixing import SumCalculator
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import FixSymmetry
from ase.data import atomic_masses
from ase.geometry import find_mic
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.optimize import BFGS
from ase.parallel import broadcast
from ase.units import Bohr, GPa, Hartree
from gpaw import GPAW, PW
from gpaw.mpi import world
from gpaw.utilities.memory import MemNode
from scipy.optimize import linear_sum_assignment


ROOT = Path(__file__).resolve().parents[1]
ROTATIONS = np.array([np.eye(3), np.diag([-1, -1, 1]),
                      -np.eye(3), np.diag([1, 1, -1])], dtype=int)
TRANSLATIONS = np.array([[0, 0, 0], [.5, .5, .5], [0, 0, 0], [.5, .5, .5]])


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def encode(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def root_io(action):
    error = None
    if world.rank == 0:
        try:
            action()
        except Exception:
            error = traceback.format_exc()
    error = broadcast(error, root=0, comm=world)
    if error:
        raise RuntimeError(error)


def save(path, obj, append=False):
    def action():
        text = json.dumps(obj, default=encode, allow_nan=False,
                          indent=None if append else 2) + '\n'
        if append:
            with Path(path).open('a', encoding='utf-8') as handle:
                handle.write(text)
                handle.flush()
        else:
            tmp = Path(str(path) + '.tmp')
            tmp.write_text(text, encoding='utf-8')
            tmp.replace(path)
    root_io(action)


def resolve(value, input_dir):
    path = Path(value)
    return path.resolve() if path.is_absolute() else (input_dir / path).resolve()


def all_ranks(value):
    return [broadcast(value if world.rank == rank else None, root=rank, comm=world)
            for rank in range(world.size)]


def memory_state():
    mem = {}
    for line in Path('/proc/meminfo').read_text().splitlines():
        key, value = line.split(':', 1)
        if key in ('MemTotal', 'MemAvailable', 'SwapFree', 'SwapTotal'):
            mem[key + '_bytes'] = int(value.split()[0]) * 1024
    page_size = os.sysconf('SC_PAGE_SIZE')
    mem['rank_rss_bytes'] = int(Path('/proc/self/statm').read_text().split()[1]) * page_size
    mem['rank_peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    return mem


def operation_closure(atoms):
    """Declared c-unique operations, species-preserving bijection, true MIC."""
    scaled = atoms.get_scaled_positions(wrap=False)
    residuals = []
    for rotation, translation in zip(ROTATIONS, TRANSLATIONS):
        transformed = scaled @ rotation.T + translation
        worst = 0.0
        for number in set(atoms.numbers):
            ids = np.flatnonzero(atoms.numbers == number)
            delta = (transformed[ids, None] - scaled[None, ids]) @ atoms.cell.array
            _, dist = find_mic(delta.reshape(-1, 3), atoms.cell, pbc=True)
            dist = dist.reshape(len(ids), len(ids))
            row, col = linear_sum_assignment(dist)
            worst = max(worst, float(dist[row, col].max()))
        residuals.append(worst)
    return {'max_residual_A': max(residuals), 'per_operation_residual_A': residuals}


def symmetry(atoms, symprec):
    ds = spglib.get_symmetry_dataset((atoms.cell.array, atoms.get_scaled_positions(),
                                     atoms.numbers), symprec=symprec)
    return {'detection_succeeded': ds is not None,
            'number': int(ds.number) if ds else None,
            'international': ds.international if ds else None,
            'hall_number': int(ds.hall_number) if ds else None,
            'choice': ds.choice if ds else None,
            'declared_operations': operation_closure(atoms)}


def validate(spec, input_dir):
    required = ('branch_id', 'candidate_id', 'structure', 'structure_sha256',
                'candidate_metadata', 'validation_report', 'output_directory',
                'method', 'ecut_eV', 'kmesh', 'scf_energy_eV_electron',
                'scf_density_electron_electron', 'scf_eigenstates_eV2_electron',
                'max_scf_iterations', 'mode', 'resource_limits')
    for key in required:
        if key not in spec:
            raise ValueError('Explicit pilot input required: ' + key)
    if spec['branch_id'] != 'O' or spec['mode'] not in ('single_point', 'relax'):
        raise ValueError('Only ordered O fixed-cell single_point/relax is implemented')
    if spec['method'] not in ('PBE', 'PBE-D3(BJ)+ATM'):
        raise ValueError('This pilot implements PBE or PBE-D3(BJ)+ATM')
    if spec.get('primary_isotope', 'h10') != 'h10' or spec.get('variable_cell', False):
        raise ValueError('Only ordinary h10, fixed-cell pilot is implemented')
    if spec.get('production_parameters_validated', False):
        raise ValueError('Pilot runner cannot declare validated PA6 production parameters')
    for key in ('ecut_eV', 'scf_energy_eV_electron', 'scf_density_electron_electron',
                'scf_eigenstates_eV2_electron'):
        if not np.isfinite(spec[key]) or float(spec[key]) <= 0:
            raise ValueError(key + ' must be finite and positive')
    mesh = np.array(spec['kmesh'])
    if mesh.shape != (3,) or np.any(mesh != mesh.astype(int)) or np.any(mesh < 1):
        raise ValueError('kmesh must contain three positive integers')
    if not 1 <= int(spec['max_scf_iterations']) <= 500:
        raise ValueError('max_scf_iterations outside 1..500')
    if spec.get('gpts') is not None and spec.get('h_A') is not None:
        raise ValueError('Use gpts or h_A, not both')
    if spec['mode'] == 'relax':
        for key in ('fmax_eV_A', 'max_steps', 'maxstep_A'):
            if key not in spec:
                raise ValueError('Explicit relaxation input required: ' + key)
        if not 0 < spec['maxstep_A'] <= .1 or not 0 < spec['fmax_eV_A'] <= .1:
            raise ValueError('Invalid pilot fmax/maxstep')
        if not 1 <= spec['max_steps'] <= 200:
            raise ValueError('Pilot max_steps must be 1..200')
    limits = spec['resource_limits']
    for key in ('max_mpi_ranks', 'estimated_total_memory_limit_GiB',
                'minimum_available_memory_GiB', 'max_wall_seconds'):
        if key not in limits or limits[key] <= 0:
            raise ValueError('Missing/invalid resource limit: ' + key)
    if world.size > limits['max_mpi_ranks']:
        raise ValueError('MPI allocation exceeds declared resource limit')
    structure = resolve(spec['structure'], input_dir)
    metadata_path = resolve(spec['candidate_metadata'], input_dir)
    report_path = resolve(spec['validation_report'], input_dir)
    if sha(structure) != spec['structure_sha256']:
        raise ValueError('Input structure SHA256 mismatch')
    report = json.loads(report_path.read_text(encoding='utf-8-sig'))
    if not (report.get('passed') is True or report.get('status') == 'passed'):
        raise ValueError('Independent geometry validation has not passed')
    if report.get('structure_sha256') != sha(structure):
        raise ValueError('Independent validator report is not bound to this structure')
    metadata = json.loads(metadata_path.read_text(encoding='utf-8-sig'))
    if report.get('input_sha256') is not None and report['input_sha256'] != sha(metadata_path):
        raise ValueError('Independent validation report refers to different candidate metadata')
    atoms = read(structure)
    if not all(atoms.pbc) or not np.isfinite(atoms.positions).all() or atoms.get_volume() <= 0:
        raise ValueError('Finite, fully periodic input required')
    counts = Counter(atoms.get_chemical_symbols())
    nrepeat = counts.get('N', 0)
    if nrepeat < 1 or counts != Counter({'C': 6*nrepeat, 'H': 11*nrepeat,
                                       'N': nrepeat, 'O': nrepeat}):
        raise ValueError('Input is not integer C6H11NO repeats')
    if not np.allclose(atoms.get_masses(), atomic_masses[atoms.numbers], atol=.05, rtol=0):
        raise ValueError('Nonstandard isotope masses found; only ordinary h10 is allowed')
    roles = metadata.get('atom_roles', [])
    if len(roles) != len(atoms) or sorted(x.get('index', -1) for x in roles) != list(range(len(atoms))):
        raise ValueError('Complete atom role map indexed in structure order is required')
    if not metadata.get('bonds'):
        raise ValueError('Periodic covalent bond graph missing from candidate metadata')
    if metadata.get('symbols') is not None and metadata['symbols'] != atoms.get_chemical_symbols():
        raise ValueError('Metadata and structure atom ordering disagree')
    metadata_positions = metadata.get('positions_A', metadata.get('positions'))
    if metadata_positions is not None and not np.allclose(metadata_positions, atoms.positions, atol=1e-6):
        raise ValueError('Metadata geometry differs from structure')
    if metadata.get('branch_id', 'O') != 'O':
        raise ValueError('Candidate metadata is from another branch')
    closure = operation_closure(atoms)
    if closure['max_residual_A'] > spec.get('symprec_A', 1e-5):
        raise ValueError('Input does not satisfy declared P1121/n operations')
    return atoms, metadata, report, nrepeat, (structure, metadata_path, report_path)


def make_dft(spec, work, restart=None):
    kwargs = dict(mode=PW(float(spec['ecut_eV'])), xc='PBE', setups='paw',
                  kpts={'size': tuple(spec['kmesh']), 'gamma': True},
                  occupations={'name': 'fermi-dirac', 'width': 0.0},
                  nbands=spec.get('nbands', '120%'), spinpol=False,
                  convergence={'energy': spec['scf_energy_eV_electron'],
                               'density': spec['scf_density_electron_electron'],
                               'eigenstates': spec['scf_eigenstates_eV2_electron']},
                  symmetry={'point_group': spec.get('point_group', True), 'time_reversal': True},
                  parallel={'augment_grids': True}, maxiter=spec['max_scf_iterations'],
                  txt=str(work / 'gpaw.txt'), legacy_gpaw=True)
    if spec.get('h_A') is not None:
        kwargs['h'] = float(spec['h_A'])
    if spec.get('gpts') is not None:
        kwargs['gpts'] = tuple(spec['gpts'])
    if spec.get('xc_stencil') is not None:
        kwargs['xc'] = {'name': 'PBE', 'stencil': int(spec['xc_stencil'])}
    return GPAW(str(restart), **kwargs) if restart else GPAW(**kwargs)


def make_dispersion(spec):
    if spec['method'] == 'PBE':
        return None
    from dftd3.ase import DFTD3
    return DFTD3(damping='d3bj', params_tweaks={'method': 'PBE', 'atm': True},
                 realspace_cutoff={'disp2': 60*Bohr, 'disp3': 40*Bohr, 'cn': 60*Bohr,
                                   'width2': .05*Bohr, 'width3': .05*Bohr})


def datasets(dft):
    result = {}
    for setup in dft.wfs.setups:
        filename = getattr(setup, 'filename', None)
        result[setup.symbol] = {'filename': filename, 'type': setup.type,
                                'fingerprint': getattr(setup, 'fingerprint', None),
                                'sha256': sha(filename) if filename and Path(filename).is_file() else None}
    if set(result) != {'C', 'H', 'N', 'O'} or any(x['sha256'] is None for x in result.values()):
        raise ValueError('C/H/N/O PAW dataset provenance incomplete')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_json', type=Path)
    parser.add_argument('--check-input', action='store_true')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    path = args.input_json.resolve()
    spec = json.loads(path.read_text(encoding='utf-8-sig'))
    atoms, metadata, report, nrepeat, sources = validate(spec, path.parent)
    if args.check_input:
        if world.rank == 0:
            print(json.dumps({'status': 'input_checks_passed', 'natoms': len(atoms),
                              'repeat_units': nrepeat, 'symmetry': symmetry(atoms, 1e-5)}))
        return
    work = resolve(spec['output_directory'], path.parent)
    # flock releases on normal exit or process death; lock metadata is not proof
    # of a running task. This is one host-global lock shared by project attempts.
    lock_handle = None
    error = None
    if world.rank == 0:
        try:
            lock_path = Path('/tmp/pa6_ordered_local_pilot.lock')
            lock_handle = lock_path.open('a+')
            fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            lock_handle.seek(0)
            lock_handle.truncate()
            json.dump({'pid': os.getpid(), 'input': str(path), 'output': str(work)}, lock_handle)
            lock_handle.flush()
            work.mkdir(parents=True, exist_ok=False)
        except Exception:
            error = traceback.format_exc()
    error = broadcast(error, root=0, comm=world)
    if error:
        raise RuntimeError(error)
    world.barrier()
    started = time.time()
    trajectory = dft = None
    latest = {}
    state = {'status': 'initializing', 'input': spec, 'input_sha256': sha(path),
             'runner_sha256': sha(__file__), 'branch_id': 'O', 'candidate_id': spec['candidate_id'],
             'mpi_ranks': world.size, 'natoms': len(atoms), 'repeat_units': nrepeat,
             'formula': atoms.get_chemical_formula(), 'primary_isotope': 'h10',
             'start_unix': started, 'production_parameters_validated': False,
             'start_type': 'restart' if spec.get('restart_gpw') else 'cold',
             'parameter_scope': 'resource_and_numerical_pilot_not_production',
             'initial_symmetry': symmetry(atoms, spec.get('symprec_A', 1e-5)),
             'source_sha256': {str(p): sha(p) for p in sources},
             'units': {'energy': 'eV/cell and eV/C6H11NO repeat', 'force': 'eV/Angstrom',
                       'stress': 'eV/Angstrom^3, ASE tensile positive',
                       'stress_voigt_order': ['xx', 'yy', 'zz', 'yz', 'xz', 'xy']},
             'software': {'python': sys.version, 'platform': platform.platform(),
                          **{p: importlib.metadata.version(p) for p in
                             ('gpaw', 'ase', 'numpy', 'scipy', 'spglib', 'dftd3')}}}
    code = 1
    try:
        for source, name in [(path, 'input_snapshot.json'),
                             (Path(__file__), 'runner_snapshot.py'),
                             (sources[1], 'candidate_metadata_snapshot.json'),
                             (sources[2], 'validation_snapshot.json')]:
            root_io(lambda s=source, n=name: shutil.copy2(s, work / n))
        write(work / 'input_structure.extxyz', atoms)
        save(work / 'provenance.json', state)
        save(work / 'status.json', state)
        fixed_cell = atoms.cell.array.copy()
        if spec.get('constrain_symmetry', True):
            before = atoms.positions.copy()
            atoms.set_constraint(FixSymmetry(atoms, symprec=spec.get('symprec_A', 1e-5),
                                             adjust_cell=False, verbose=False))
            shift = np.max(np.linalg.norm(atoms.positions-before, axis=1))
            if shift > 1e-6:
                raise ValueError('Symmetry preparation moved input atoms by >1e-6 A')
            state['symmetry_preparation_max_displacement_A'] = float(shift)
        restart = resolve(spec['restart_gpw'], path.parent) if spec.get('restart_gpw') else None
        if restart and sha(restart) != spec.get('restart_sha256'):
            raise ValueError('Restart SHA256 mismatch')
        # Estimates are based on the input geometry, without loading wavefunctions.
        dft = make_dft(spec, work, restart=None if args.dry_run else restart)
        if restart and not args.dry_run:
            old = dft.get_atoms()
            if not (np.array_equal(old.numbers, atoms.numbers) and
                    np.allclose(old.cell, atoms.cell, atol=1e-10, rtol=0) and
                    np.allclose(old.positions, atoms.positions, atol=1e-8, rtol=0)):
                raise ValueError('Restart geometry/order differs from validated input')
            state['restart_sha256'] = sha(restart)
        dft.initialize(atoms)
        node = MemNode('Calculator')
        dft.estimate_memory(node)
        estimate = {'rank': world.rank, 'calculator_bytes': node.calculate_size(),
                    'process': memory_state()}
        estimates = all_ranks(estimate)
        estimated_bytes = sum(x['calculator_bytes'] + x['process']['rank_rss_bytes'] for x in estimates)
        state['memory_estimate'] = {'by_rank': estimates, 'sum_with_current_rss_bytes': estimated_bytes,
                                    'caveat': 'Estimate plus current RSS may double count some allocations; not a measured peak.'}
        state['datasets'] = datasets(dft)
        state['basis'] = {'bands': int(dft.wfs.bd.nbands), 'valence_electrons': float(dft.wfs.nvalence),
                          'grid': dft.wfs.gd.N_c.tolist(),
                          'bz_kpoints': dft.wfs.kd.bzk_kc.tolist(), 'ibz_kpoints': dft.wfs.kd.ibzk_kc.tolist()}
        save(work / 'memory_estimate.json', state['memory_estimate'])
        save(work / 'datasets.json', state['datasets'])
        save(work / 'provenance.json', state)
        if estimated_bytes > spec['resource_limits']['estimated_total_memory_limit_GiB'] * 2**30:
            raise MemoryError('Conservative memory estimate exceeds declared local budget')
        if args.dry_run:
            state['status'], code = 'passed_initialization_no_scf', 0
            return
        dispersion = make_dispersion(spec)
        calculations = [dft] + ([dispersion] if dispersion is not None else [])
        context = {'evaluation': -1, 'optimizer': None, 'scf_start': time.time()}
        limits = spec['resource_limits']

        def guard_and_scf_log():
            by_rank = all_ranks(memory_state())
            item = {'evaluation_id': context['evaluation'], 'SCF_iteration': int(dft.scf.niter),
                    'wall_seconds': time.time()-started,
                    'simultaneous_rank_rss_sum_bytes': sum(x['rank_rss_bytes'] for x in by_rank),
                    'host_mem_available_bytes': by_rank[0]['MemAvailable_bytes']}
            value = getattr(dft.hamiltonian, 'e_total_free', None)
            if value is not None and np.isfinite(value):
                item['dft_energy_free_eV'] = float(value*Hartree)
            save(work / 'scf_iterations.jsonl', item, append=True)
            if time.time()-started > limits['max_wall_seconds']:
                raise TimeoutError('Pilot wall budget reached at SCF boundary')
            if by_rank[0]['MemAvailable_bytes'] < limits['minimum_available_memory_GiB']*2**30:
                raise MemoryError('Available host memory below declared margin')

        dft.attach(guard_and_scf_log, 1)
        trajectory = Trajectory(work / 'trajectory.traj', 'w',
                                properties=['energy', 'free_energy', 'forces', 'stress'],
                                comm=world, master=world.rank == 0)

        class RecordedCalculator(SumCalculator):
            def calculate(self, atoms, properties, system_changes):
                current = atoms
                context['evaluation'] += 1
                point_started = time.time()
                super().calculate(current, ['energy', 'free_energy', 'forces', 'stress'], system_changes)
                raw = {p: np.array(self.results[p]).copy() for p in
                       ('energy', 'free_energy', 'forces', 'stress')}
                if not dft.scf.converged or not all(np.isfinite(x).all() for x in raw.values()):
                    raise RuntimeError('Nonfinite result or unconverged SCF')
                projected = raw['forces'].copy()
                for constraint in current.constraints:
                    constraint.adjust_forces(current, projected)
                op = operation_closure(current)
                item = {'evaluation_id': context['evaluation'],
                        'optimizer_step': getattr(context['optimizer'], 'nsteps', 0),
                        'evaluation_wall_seconds': time.time()-point_started,
                        'wall_seconds': time.time()-started,
                        'SCF_converged': bool(dft.scf.converged), 'SCF_iterations': int(dft.scf.niter),
                        'energy_eV_cell': float(raw['energy']),
                        'free_energy_eV_cell': float(raw['free_energy']),
                        'energy_eV_repeat': float(raw['energy'])/nrepeat,
                        'raw_forces_eV_A': raw['forces'], 'projected_forces_eV_A': projected,
                        'raw_fmax_eV_A': float(np.linalg.norm(raw['forces'], axis=1).max()),
                        'projected_fmax_eV_A': float(np.linalg.norm(projected, axis=1).max()),
                        'raw_stress_eV_A3': raw['stress'], 'raw_stress_kbar': raw['stress']/GPa*10,
                        'positions_A': current.positions.copy(), 'cell_A': current.cell.array.copy(),
                        'declared_operations': op,
                        'dft': {p: dft.results[p] for p in ('energy', 'free_energy', 'forces', 'stress')},
                        'dispersion': None if dispersion is None else
                        {p: dispersion.results[p] for p in ('energy', 'forces', 'stress')}}
                if not np.allclose(current.cell.array, fixed_cell, atol=1e-12, rtol=0):
                    raise RuntimeError('Fixed cell changed')
                if spec.get('constrain_symmetry', True) and op['max_residual_A'] > spec.get('symprec_A', 1e-5):
                    raise RuntimeError('Declared ordered symmetry was lost')
                save(work / 'evaluations.jsonl', item, append=True)
                save(work / f"evaluation_{context['evaluation']:04d}.json", item)
                snapshot = current.copy()
                snapshot.set_constraint([])
                snapshot.calc = SinglePointCalculator(snapshot, **raw)
                trajectory.write(snapshot)
                write(work / 'latest_structure.extxyz', snapshot, write_results=False)
                if spec.get('checkpoint_each_evaluation', True):
                    dft.write(str(work / 'latest.gpw.tmp'), mode='all')
                    root_io(lambda: (work / 'latest.gpw.tmp').replace(work / 'latest.gpw'))
                    world.barrier()
                latest.clear()
                latest.update(item)
                save(work / 'status.json', {'status': 'running_pilot', 'last_evaluation': context['evaluation'],
                                           'raw_fmax_eV_A': item['raw_fmax_eV_A'],
                                           'wall_seconds': item['wall_seconds']})

        atoms.calc = RecordedCalculator(calculations)
        if spec['mode'] == 'single_point':
            atoms.get_potential_energy()
            state['status'], code = 'completed_single_point_unvalidated_parameters', 0
        else:
            optimizer = BFGS(atoms, logfile=str(work / 'optimizer.log'),
                             restart=str(work / 'optimizer_restart.json'),
                             maxstep=spec['maxstep_A'], comm=world)
            context['optimizer'] = optimizer
            converged = bool(optimizer.run(fmax=spec['fmax_eV_A'], steps=spec['max_steps']))
            raw_ok = latest['raw_fmax_eV_A'] <= spec['fmax_eV_A']
            state['optimizer_converged'] = converged
            state['raw_force_target_met'] = raw_ok
            state['optimizer_steps'] = optimizer.nsteps
            state['status'] = ('pilot_force_target_met_unvalidated_parameters' if converged and raw_ok
                               else 'pilot_step_limit_or_raw_force_target_not_met')
            # A bounded short pilot that did not converge is still a completed
            # measurement; the scientific acceptance field remains false.
            code = 0
        state['final_evaluation'] = latest
        state['final_symmetry'] = symmetry(atoms, spec.get('symprec_A', 1e-5))
        final = atoms.copy()
        final.set_constraint([])
        write(work / 'final_structure.extxyz', final)
        dft.write(str(work / 'final.gpw'), mode='all')
        state['final_geometry_sha256'] = sha(work / 'final_structure.extxyz')
        state['checkpoint_sha256'] = sha(work / 'final.gpw')
        state['production_acceptance_passed'] = False
    except Exception as exc:
        state['status'], code = 'failed_pilot', 1
        state['error'], state['traceback'] = repr(exc), traceback.format_exc()
        state['last_completed_evaluation'] = latest
        if world.rank == 0:
            print(state['traceback'], file=sys.stderr, flush=True)
    finally:
        if trajectory is not None:
            trajectory.close()
        state['wall_seconds'], state['exit_code'] = time.time()-started, code
        state['final_memory_by_rank'] = all_ranks(memory_state())
        save(work / 'result.json', state)
        save(work / 'status.json', {k: state[k] for k in ('status', 'wall_seconds', 'exit_code')})
        if dft is not None:
            dft.close()
        if lock_handle is not None:
            lock_handle.close()
        if world.rank == 0:
            print(json.dumps({'status': state['status'], 'output': str(work),
                              'wall_seconds': state['wall_seconds'], 'exit_code': code}))
    if code:
        raise SystemExit(code)


if __name__ == '__main__':
    main()
