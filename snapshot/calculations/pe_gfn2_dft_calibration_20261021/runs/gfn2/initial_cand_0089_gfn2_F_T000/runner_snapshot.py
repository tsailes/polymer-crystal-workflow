"""Independent cold-start GFN2-xTB periodic single point; no geometry optimization.

Run with the isolated tblite/ASE environment. Native output is captured by the
parent process. Each input directory is one immutable attempt.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys
import time
import traceback


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, obj):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_json')
    args = parser.parse_args()
    inp = Path(args.input_json).resolve()
    spec = json.loads(inp.read_text(encoding='utf-8'))
    work = inp.parent
    # Set before numpy, ASE, BLAS, or the tblite extension is imported.
    threads = int(spec.get('threads', 4))
    os.environ.update(OMP_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS='1',
                      MKL_NUM_THREADS='1', OMP_MAX_ACTIVE_LEVELS='1')
    started = time.time()
    result = dict(schema_version=1, status='running', input=spec,
                  input_sha256=digest(inp), runner_sha256=digest(__file__),
                  start_unix=started, created_utc=datetime.now(timezone.utc).isoformat(),
                  geometry_optimized=False, cold_start=True,
                  environment_overrides={k: os.environ[k] for k in
                     ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_MAX_ACTIVE_LEVELS')})
    dump(work / 'status.json', result)
    try:
        import numpy as np
        import ase
        import tblite.ase
        import tblite.interface
        from ase.io import read, write
        from ase.units import Hartree, kB, GPa
        from tblite.ase import TBLite

        repeat = tuple(spec.get('repeat', [2, 3, 8]))
        if len(repeat) != 3 or any(not isinstance(v, int) or v < 1 for v in repeat):
            raise ValueError('repeat must contain three positive integers')
        temperature = float(spec['electronic_temperature_K'])
        if temperature < 0 or not np.isfinite(temperature):
            raise ValueError('Electronic temperature must be finite and nonnegative')
        if threads < 1 or int(spec.get('charge', 0)) != 0 or int(spec.get('multiplicity', 1)) != 1:
            raise ValueError('This PE worker accepts neutral singlet structures only')
        structure = Path(spec['structure']).resolve()
        primitive = read(structure)
        if len(primitive) != 12 or primitive.get_chemical_formula() != 'C4H8' or not primitive.pbc.all():
            raise ValueError('Expected a periodic twelve-atom C4H8 primitive cell')
        if primitive.constraints:
            raise ValueError('Frozen-geometry comparison requires a constraint-free source structure')
        primitive.calc = None
        primitive.set_initial_charges(np.zeros(len(primitive)))
        primitive.set_initial_magnetic_moments(np.zeros(len(primitive)))
        replicas = int(np.prod(repeat))
        nch2 = len(primitive) // 3
        supercell = primitive.repeat(repeat)
        settings = dict(method='GFN2-xTB', charge=0.0, multiplicity=1,
                        accuracy=float(spec.get('accuracy', .001)),
                        max_iterations=int(spec.get('max_iterations', 1000)),
                        electronic_temperature=temperature, verbosity=int(spec.get('verbosity', 1)),
                        cache_api=False)
        result.update(software={name: importlib.metadata.version(name) for name in ('tblite', 'ase', 'numpy')},
                      python=sys.version, platform=platform.platform(), settings=settings,
                      source_geometry=str(structure), source_geometry_sha256=digest(structure),
                      primitive_natoms=len(primitive), supercell_natoms=len(supercell),
                      n_CH2_per_primitive_cell=nch2, replica_count=replicas,
                      cell_A=primitive.cell.array.tolist(), positions_A=primitive.positions.tolist(),
                      atomic_numbers=primitive.numbers.tolist(), repeat=list(repeat),
                      adapter_source_path=tblite.ase.__file__, adapter_source_sha256=digest(tblite.ase.__file__),
                      interface_source_path=tblite.interface.__file__, interface_source_sha256=digest(tblite.interface.__file__))
        canonical = dict(numbers=primitive.numbers.tolist(),
                         cell_A=np.round(primitive.cell.array, 12).tolist(),
                         positions_A=np.round(primitive.positions, 12).tolist(), pbc=primitive.pbc.tolist())
        result['canonical_geometry_sha256'] = hashlib.sha256(json.dumps(canonical, sort_keys=True).encode()).hexdigest()
        write(work / 'input_structure.extxyz', primitive)
        write(work / 'supercell_input.extxyz', supercell)
        calc = TBLite(**settings)
        supercell.calc = calc
        # One explicit API calculation; reading calc.results prevents a second
        # singlepoint caused by ASE's free_energy property dispatch.
        calc.calculate(supercell, properties=['energy', 'forces', 'stress', 'charges'])
        energy = float(calc.results['energy'])
        free_energy = float(calc.results['free_energy'])
        all_forces = np.asarray(calc.results['forces'], dtype=float)
        stress = np.asarray(calc.results['stress'], dtype=float)
        all_charges = np.asarray(calc.results['charges'], dtype=float)
        folded = all_forces.reshape(replicas, len(primitive), 3)
        mean_force = folded.mean(axis=0)
        replica_force = float(np.linalg.norm(folded - mean_force, axis=2).max())
        folded_charges = all_charges.reshape(replicas, len(primitive))
        if not all(np.isfinite(value).all() for value in [energy, free_energy, all_forces, stress, all_charges]):
            raise ValueError('Non-finite quantum result')
        if replica_force > float(spec.get('replica_force_tolerance_eV_A', 1e-5)):
            raise ValueError(f'Replica force disagreement {replica_force:g} eV/A exceeds tolerance')
        occupations = np.asarray(calc._res.get('orbital-occupations'))
        orbital_eV = np.asarray(calc._res.get('orbital-energies')) * Hartree
        # Neutral singlet adapter supplies non-spin-polarized occupancies 0..2.
        if occupations.ndim != 1 or np.any(occupations < -1e-10) or np.any(occupations > 2 + 1e-10):
            raise ValueError('Unexpected orbital occupation convention')
        probability = np.clip(occupations / 2.0, 0.0, 1.0)
        middle = (probability > 0) & (probability < 1)
        entropy_kB = float(-2*np.sum(probability[middle]*np.log(probability[middle])
                           +(1-probability[middle])*np.log1p(-probability[middle])))
        # This is an occupation-derived diagnostic, not an independently
        # provided native internal-energy component.
        occupation_TS_eV = temperature * kB * entropy_kB
        result.update(status='passed_single_point', scf_converged=True,
                      scf_convergence_evidence='tblite singlepoint API returned without its SCC convergence exception; native cycle output is retained externally',
                      units={'energy': 'eV', 'forces': 'eV/angstrom', 'stress': 'eV/angstrom^3; ASE tensile positive',
                             'stress_voigt_order': ['xx', 'yy', 'zz', 'yz', 'xz', 'xy']},
                      energy_semantics='tblite ASE assigns energy and free_energy from the same native energy scalar; no independent native internal-energy scalar is exposed. Occupation-derived T*S is recorded separately.',
                      electronic_temperature_K=temperature,
                      supercell=dict(energy_eV=energy, free_energy_eV=free_energy,
                          forces_eV_A=all_forces.tolist(), stress_eV_A3=stress.tolist(),
                          stress_kbar=(stress / GPa * 10).tolist(), charges_e=all_charges.tolist()),
                      total=dict(energy_eV_cell=energy/replicas, energy_eV_per_cell=energy/replicas,
                          energy_free_eV_cell=free_energy/replicas,
                          energy_eV_CH2=energy/replicas/nch2, energy_free_eV_CH2=free_energy/replicas/nch2,
                          forces_eV_A=mean_force.tolist(), max_force_norm_eV_A=float(np.linalg.norm(mean_force,axis=1).max()),
                          stress_eV_A3=stress.tolist(), stress_kbar=(stress/GPa*10).tolist(),
                          mean_charges_e=folded_charges.mean(axis=0).tolist()),
                      replica_force_max_deviation_eV_A=replica_force,
                      replica_charge_max_deviation_e=float(np.abs(folded_charges-folded_charges.mean(axis=0)).max()),
                      occupation_diagnostic=dict(orbital_energies_eV=orbital_eV.tolist(), orbital_occupations_e=occupations.tolist(),
                          occupation_sum_e=float(occupations.sum()), entropy_supercell_kB=entropy_kB,
                          TS_supercell_eV=occupation_TS_eV, TS_eV_CH2=occupation_TS_eV/replicas/nch2),
                      wall_seconds=time.time()-started)
        dump(work/'result.json',result)
        dump(work/'status.json',result)
        print(json.dumps({'status':result['status'],'id':spec.get('id'),'energy_eV_CH2':result['total']['energy_eV_CH2'],
                          'replica_force_deviation_eV_A':replica_force,'temperature_K':temperature,'wall_seconds':result['wall_seconds']}),flush=True)
    except Exception as exc:
        result.update(status='failed_single_point',scf_converged=False,error=f'{type(exc).__name__}: {exc}',
                      traceback=traceback.format_exc(),wall_seconds=time.time()-started)
        dump(work/'result.json',result)
        dump(work/'status.json',result)
        print(result['traceback'],file=sys.stderr,flush=True)
        raise


if __name__ == '__main__':
    main()
