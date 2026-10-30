"""Read-only geometry/settings reuse audit; starts no electronic calculations."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone
import numpy as np
from ase.io import read
from ase.geometry import find_mic

BASE = Path('/mnt/e/PE_DFT/calculations')
OUT = BASE/'pe_gfn2_dft_calibration_20261021/audit/design_review.json'
OLD = BASE/'pe_gpaw_convergence_20261020'
GROUP = BASE/'pe_pnam_candidates_20261020/representatives/group_001/POSCAR'

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def compare(reference, positions, cell):
    cell = np.asarray(cell)
    assert cell.shape == (3, 3)
    delta = np.asarray(positions) - reference.positions
    _, lengths = find_mic(delta, reference.cell, pbc=True)
    return {
        'same_cell': bool(np.allclose(cell, reference.cell.array, rtol=0, atol=1e-12)),
        'cell_max_abs_difference_A': float(np.max(abs(cell-reference.cell.array))),
        'identity_atom_map': list(range(len(reference))),
        'periodic_identity_map_RMS_A': float(np.sqrt(np.mean(lengths**2))),
        'periodic_identity_map_max_displacement_A': float(np.max(lengths)),
        'unwrapped_max_abs_cartesian_difference_A': float(np.max(abs(delta))),
        'requires_strain_or_axis_change': bool(not np.allclose(cell, reference.cell.array, rtol=0, atol=1e-12)),
    }

def main():
    group = read(GROUP)
    old_poscar = OLD/'reference/POSCAR_GFN2'
    old = read(old_poscar, format='vasp')
    lock_path = OLD/'releases/20261021T043247_514339Z/parameters_lock.json'
    lock = json.loads(lock_path.read_text())
    rows = []
    for tag in ['pbe_d3bj', 'optb88_vdw', 'r2scan_d4']:
        run = OLD/f'runs/transfer_gfn2_{tag}_H'
        input_path, result_path = run/'input.json', run/'result.json'
        spec = json.loads(input_path.read_text())
        result = json.loads(result_path.read_text())
        geom = compare(group, result['positions_A'], result['cell_A'])
        settings = {
            'ecut_eV': spec['ecut_eV'] == 1800,
            'kmesh': spec['kmesh'] == [4, 8, 12],
            'gpts': spec['gpts'] == [112, 72, 44],
            'nbands': spec.get('nbands', 16) == 16,
            'point_group': spec['point_group'] is True,
            'mpi_ranks': spec['mpi_ranks'] == 4,
            'xc_stencil': spec['xc_stencil'] == 4,
            'scf_energy': spec['scf_energy_eV_electron'] == 1e-9,
            'scf_density': spec['scf_density_electron_electron'] == 1e-9,
            'scf_eigenstates': spec['scf_eigenstates_eV2_electron'] == 1e-11,
            'input_file_hash': sha(input_path) == result['input_sha256'],
            'source_geometry_hash': sha(old_poscar) == result['source_geometry_sha256'],
            'source_runner_hash': sha(OLD/'scripts/run_point.py') == result['runner_sha256'],
            'protocol_hash': sha(OLD/'protocol.json') == spec['protocol_sha256'],
            'PAW_C_H_hashes': all(item['file_sha256'] == lock['physics']['paw_sha256'][item['symbol']] for item in result['datasets']),
            'process_success': json.loads((run/'process.json').read_text())['exit_code'] == 0,
            'result_passed': result['status'] == 'passed_single_point' and result['SCF_converged'] is True,
        }
        reusable = all(settings.values()) and geom['same_cell'] and geom['periodic_identity_map_max_displacement_A'] < 1e-12
        rows.append({'method':spec['method'], 'run_directory':str(run), 'input_sha256':sha(input_path), 'result_sha256':sha(result_path),
                     'settings_checks': settings, 'geometry_comparison':geom, 'reusable_as_final_group001_H':reusable,
                     'free_energy_eV_cell':result['total']['energy_free_eV'],
                     'energy_definition':result['units']['energy_force_consistent']})
    out = {
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'status':'read_only_independent_design_audit',
        'quantum_calculations_started': False,
        'old_scientific_inputs_modified': False,
        'proposed_scope':'Pilot 21 fixed-cell PE geometries: 16 initial plus 4 final groups plus zero-angle final negative control',
        'group001_POSCAR':str(GROUP), 'group001_POSCAR_sha256':sha(GROUP),
        'old_reference_POSCAR':str(old_poscar), 'old_reference_POSCAR_sha256':sha(old_poscar),
        'old_reference_comparison':compare(group, old.positions, old.cell.array),
        'same_element_order':bool(np.array_equal(group.numbers,old.numbers)),
        'existing_H_results':rows,
        'all_three_H_points_reusable':all(row['reusable_as_final_group001_H'] for row in rows),
        'design_recommendations':[
            'Reuse exact matching final_group001 H results and provenance rather than rerun identical points.',
            'PBE-D3(BJ) is the primary DFT target; other two methods are separately reported cross-method checks.',
            'Initial structures and near-minimum final structures require separate residual summaries; do not pool all 21 as iid samples.',
            'Angle-block holdout is a stress test for interpolation/extrapolation inside one constrained angle family, not independent polymorph validation.',
            'Only four near-degenerate final groups: estimate descriptive residual spread but retain all unless supported by external validation; no 99.7 percent guarantee.',
            'Initial0, initial75 and final4 are useful H/U/L numerical transfer endpoints; add R0 U only if precision uncertainty relative to target gaps justifies it.',
            'Calibrate relative-energy residuals against fixed R0 directly; single-structure error sigma cannot be substituted without covariance treatment.',
            'Freeze an explicit low-energy window W, report it as a research decision, and only consider rejection when predicted relative energy minus validated margin exceeds W.',
            'Keep existing unstable zero-angle final as diagnostic negative control; do not silently relabel it a 3sigma rejection.'
        ],
        'electronic_temperature_semantics':{
            'old_GFN2_temperature_K':300.0, 'old_DFT_occupation_width_eV':0.0,
            'local_ASE_interface':'tblite.ase sets energy=native energy*Hartree and free_energy=energy; equality is interface aliasing, not proof of zero entropy.',
            'official_tblite_0p7_source':'https://raw.githubusercontent.com/tblite/tblite/v0.7.0/src/tblite/scf/iterator.f90',
            'source_evidence':'next_scf initializes energies with ts/natom; next_density evaluates electronic entropy and sums spin terms; native energy includes the electronic entropy contribution.',
            'required_gate':'Independent 0/300 K pilot on fixed geometries, finite E/F/stress, occupation/entropy or near-zero-temperature check; use same energy semantics thereafter.'
        },
        'reuse_code_audit':{
            'existing_formal_sigma_executor_found':False,
            'reusable_sources':['DFT frozen run_point.py and dispatcher pattern','GFN2 validated exact-replica folding and geometry audit'],
            'excluded_sources':['Cancelled DREIDING/DFTB-versus-GFN2 production executor; remains cancelled'],
        },
        'source_files':[{'path':str(p),'sha256':sha(p)} for p in [GROUP,old_poscar,lock_path,OLD/'scripts/run_point.py',BASE/'pe_pnam_candidates_20261020/optimization_summary.json',BASE/'pe_gfn2_first_run_20261020/convergence_relaxation/workflow.py']],
    }
    OUT.parent.mkdir(exist_ok=True, parents=True)
    OUT.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'all_three_H_points_reusable':out['all_three_H_points_reusable'], 'old_reference_comparison':out['old_reference_comparison'], 'methods':[{k:r[k] for k in ['method','reusable_as_final_group001_H','geometry_comparison']} for r in rows]}))

if __name__ == '__main__':
    main()
