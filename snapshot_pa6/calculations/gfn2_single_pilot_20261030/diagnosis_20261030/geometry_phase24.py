"""All 24 approved phases at fixed template0/x=.25/phi=0, diagnostic only."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import spglib
from ase import Atoms
from ase.io import write
from ase.geometry import cellpar_to_cell

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[2]
sys.path.insert(0,str(PROJECT/'scripts'))
import build_ordered_candidates as b
import validate_ordered as v

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def load(path): return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def dump(path,record): path.write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

def main():
    out=HERE/'geometry_phase24';out.mkdir(exist_ok=True)
    paths=[PROJECT/'inputs/allowed_literature_inputs.json',PROJECT/'inputs/independent_template/chain_00.json',
           PROJECT/'inputs/chemistry_template_manifest.json',PROJECT/'branches/O/phase24_first_pass/selection.json',
           Path(b.__file__),Path(v.__file__),Path(__file__)]
    meta,template,manifest,selection=map(load,paths[:4])
    for q in manifest['source_hashes']:
        assert sha(PROJECT/q['path'])==q['sha256']
    assert selection['phase_fractional_c']==(np.arange(24)/48).tolist()
    cell=cellpar_to_cell(meta['cell']['lengths_A']+meta['cell']['angles_deg'])
    protocol={'source_sha256':{str(p):sha(p) for p in paths},'template_index':0,
              'chain_x_fractional':.25,'chain_rotation_about_c_deg':0.,
              'phase_grid':selection['phase_fractional_c'],
              'role':'All 24 entries retained; geometry and H-bond diagnostics only. No H-bond-based filtering, ranking, selection, or quantum calculation.',
              'literature_atomic_coordinates_read':False}
    dump(HERE/'geometry_phase24_protocol.json',protocol)
    rows=[]
    for k,phase in enumerate(selection['phase_fractional_c']):
        data=b.expand(np.asarray(template['positions_unwrapped_A']),template['atom_roles'],template['bonds'],template['symbols'],cell,.25,0.,phase)
        builder=b.validate_geometry(data)
        independent=v.validate(data,analyze_hbonds=True)
        symmetry=spglib.get_symmetry_dataset((cell,np.asarray(data['fractional_positions']),[v.NUMBERS[s] for s in data['symbols']]),symprec=1e-5)
        pure_translations=[(np.asarray(t)%1).tolist() for R,t in zip(symmetry.rotations,symmetry.translations) if np.array_equal(R,np.eye(3,dtype=int))]
        bhalf=any(np.linalg.norm((np.array(t)-[0,.5,0])-np.rint(np.array(t)-[0,.5,0]))<1e-7 for t in pure_translations)
        actual={'operation_count':len(symmetry.rotations),'pure_translations_fractional':pure_translations,
                'contains_pure_b_half_translation':bhalf,
                'rotations':symmetry.rotations.tolist(),'translations':symmetry.translations.tolist(),
                'symprec_A':1e-5}
        data.update(diagnostic_id=f'geometry-phase-{k:02d}',template_id=template['template_id'],
                    branch='O',purpose='post-run diagnostic only',builder_geometry=builder,independent_validation=independent,actual_symmetry=actual)
        file=out/f'geometry_phase_{k:02d}.json';dump(file,data)
        write(out/f'geometry_phase_{k:02d}.extxyz',Atoms(data['symbols'],positions=data['positions'],cell=cell,pbc=True))
        hb=independent['hydrogen_bonds'];sheet=independent['hydrogen_bond_sheet_topology']
        row={'phase_index_24':k,'original_phase_index_48':2*k,'phase_fractional_c':phase,
             'builder_geometry_passed':builder['passed'],'independent_geometry_passed':independent['passed'],
             'spglib_number':builder['spglib_number'],'hydrogen_bond_count_diagnostic':len(hb['bonds']),
             'actual_symmetry_operation_count':actual['operation_count'],'contains_pure_b_half_translation':bhalf,
             'covered_donors_diagnostic':hb['covered_donor_count'],'covered_acceptors_diagnostic':hb['covered_acceptor_count'],
             'two_periodic_ac_sheets_diagnostic':sheet['passed'],'nonbonded_collision_count':len(independent['nonbonded_collisions']),
             'record':str(file),'sha256':sha(file)}
        rows.append(row)
    summary={'protocol':protocol,'entry_count':len(rows),'all_entries_retained':True,
             'entries':rows,'geometry_passing_count':sum(r['builder_geometry_passed'] and r['independent_geometry_passed'] for r in rows),
             'geometry_passing_with_full_HB_coverage_count_diagnostic':sum(r['builder_geometry_passed'] and r['independent_geometry_passed'] and r['covered_donors_diagnostic']==8 and r['covered_acceptors_diagnostic']==8 for r in rows),
             'geometry_passing_with_two_sheet_network_count_diagnostic':sum(r['builder_geometry_passed'] and r['independent_geometry_passed'] and r['two_periodic_ac_sheets_diagnostic'] for r in rows),
             'quantum_calculations_started':0,'new_compute_representative_selected':False,
             'interpretation':'Existence of network in a distinct same-No14 seed does not prove a continuous low-barrier pathway from the relaxed phase0 seed.'}
    dump(HERE/'geometry_phase24_summary.json',summary)
    print(json.dumps(rows,indent=2))

if __name__=='__main__':main()
