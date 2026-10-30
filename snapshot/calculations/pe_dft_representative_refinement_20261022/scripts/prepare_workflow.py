"""Freeze the four authorized PE endpoints and settings; no quantum calculation."""
from pathlib import Path
import datetime, hashlib, json, shutil

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parents[1]
OLD = PROJECT/'calculations/pe_gfn2_dft_calibration_20261021'
STAGE2 = PROJECT/'calculations/pe_gpaw_convergence_20261020'
RELEASE = STAGE2/'releases/20261021T043247_514339Z'

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p, obj):
    p.parent.mkdir(parents=True, exist_ok=True)
    content = json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n'
    if p.exists():
        if p.read_text(encoding='utf-8') != content: raise FileExistsError(p)
    else: p.write_text(content, encoding='utf-8')
def copy(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        if sha(src)!=sha(dst): raise FileExistsError(dst)
    else: shutil.copy2(src,dst)

def main():
    from ase.io import read, write
    import numpy as np
    import spglib
    lock=json.loads((RELEASE/'parameters_lock.json').read_text())
    assert lock['status']=='locked_after_all_required_checks_passed'
    assert sha(RELEASE/'parameters_lock.json')=='d3d0de948aebb52865c5e5a391dbfa4d964ae93cae654fbb1c0d6a2e606b10d1'
    point=STAGE2/'scripts/run_point.py'
    assert sha(point)==lock['code_sha256']['run_point']
    copy(point, ROOT/'scripts/run_point.py')
    copy(STAGE2/'environment/environment_verified.json',ROOT/'environment/environment_verified.json')
    copy(RELEASE/'parameters_lock.json', ROOT/'reference/stage2_parameters_lock.json')
    copy(OLD/'protocol.json',ROOT/'reference/stage3_protocol.json')
    copy(OLD/'samples.json',ROOT/'reference/stage3_samples.json')
    copy(OLD/'audit/dft_transfer_gate.json',ROOT/'reference/stage3_dft_transfer_gate.json')
    for paw in lock['environment']['PAW']:
        assert sha(paw['path'])==paw['sha256'], 'PAW changed'
    protocol={
      'schema_version':1,'batch_id':ROOT.name,
      'scope':'Four GFN2-preoptimized PE representatives; PBE-D3(BJ) fixed-cell Pnam refinement and complete-link geometric deduplication.',
      'authorization':'User requested DFT refinement of all four representatives, merge count, retention of all energy data and logs.',
      'method':'PBE-D3(BJ)', 'ecut_eV':1800,'kmesh':[4,8,12],'gpts':[112,72,44],
      'xc_stencil':4,'nbands':16,'mpi_ranks':4,'max_workers':2,
      'scf_energy_eV_electron':1e-9,'scf_density_electron_electron':1e-9,'scf_eigenstates_eV2_electron':1e-11,
      'cell_lengths_A':[7.432,4.929,2.545], 'cell_fixed':True,'cell_tolerance_A':1e-10,
      'spacegroup':{'number':62,'setting':'Pnam','hall_number':297,'symprec_A':1e-5},
      'optimizer':{'name':'LBFGS','fmax_eV_A':0.002,'maxstep_A':0.05,'maxsteps':200},
      'stress_role':'Record raw and symmetry-projected stress; not a fixed-cell convergence target.',
      'cold_final_verification':True,
      'merge':{'rms_A':0.02,'max_A':0.05,'energy_meV_CH2':1.0,'algorithm':'all six pairs; complete linkage; report nontransitive ambiguous matches'},
      'N_CH2':4,'electronic_temperature_K':0,
      'precision_provenance':'Stage3 H production settings; stage2 immutable release is preserved, not rewritten.',
      'point_runner_sha256':sha(point),'stage2_lock_sha256':sha(RELEASE/'parameters_lock.json'),
      'raw_retention':'Frozen inputs, hashes, every evaluated geometry, E/F/stress, SCF iterations, trajectory, optimizer checkpoints, stdout/stderr, failures, final gpw and cold-start verification.',
      'failure_policy':'Stop new submissions on SCF, geometry, integrity or final-force failure. Keep already-running jobs and all failed outputs. No automatic scientific threshold changes.',
      'claims':'Constrained local stationary structures and tolerance-based duplicate groups; not proof of global minima or distinct polymorphs.',
      'scope_exclusions':['No automatic cross-functional optimization in this queue','No variable-cell optimization','No sigma-based candidate rejection']}
    save(ROOT/'protocol.json',protocol)
    manifest=json.loads((OLD/'samples.json').read_text())
    samples=[];queue=[]
    for source in manifest['samples']:
        if not source['id'].startswith('final_group_'):continue
        sid=source['id'].removeprefix('final_')
        src=Path(source['structure'])
        assert sha(src)==source['structure_sha256']
        dst=ROOT/'structures'/sid/'POSCAR'
        copy(src,dst)
        a=read(dst)
        assert a.get_chemical_formula()=='C4H8' and len(a)==12 and a.pbc.all()
        assert np.max(abs(a.cell.array-np.diag(protocol['cell_lengths_A'])))<1e-12
        ds=spglib.get_symmetry_dataset((a.cell.array,a.get_scaled_positions(),a.numbers),symprec=1e-5,hall_number=297)
        assert ds.number==62 and len(ds.rotations)==8
        sample={'id':sid,'source_sample_id':source['id'],'source_candidate_id':source['source_candidate_id'],
          'source_path':str(src),'source_sha256':sha(src),'structure':str(dst),'structure_sha256':sha(dst),
          'original_group_member_count':source['group_member_count'],'initial_angle_deg':source['initial_angle_deg'],
          'natoms':len(a),'N_CH2':4,'cell_A':a.cell.array.tolist(),'spacegroup_number':int(ds.number)}
        samples.append(sample)
        spec={'method':protocol['method'],'structure':str(dst),'source_geometry_sha256_expected':sha(dst),
              'protocol_sha256_expected':sha(ROOT/'protocol.json'),'point_runner':str(ROOT/'scripts/run_point.py'),
              'point_runner_sha256':sha(point),'sample_id':sid,'point_group':True,'cell_fixed':True,
              'cell_tolerance_A':1e-10,'fmax_eV_A':0.002,'maxsteps':200,'maxstep':0.05,'symprec_A':1e-5}
        for key in ['ecut_eV','kmesh','gpts','xc_stencil','nbands','mpi_ranks','scf_energy_eV_electron','scf_density_electron_electron','scf_eigenstates_eV2_electron']:
            spec[key]=protocol[key]
        save(ROOT/'inputs'/f'{sid}.json',spec)
        queue.append({'id':sid,'input':str(ROOT/'inputs'/f'{sid}.json'),'input_sha256':sha(ROOT/'inputs'/f'{sid}.json')})
    assert len(samples)==4
    save(ROOT/'sample_manifest.json',{'schema_version':1,'samples':samples})
    save(ROOT/'queue.json',{'max_workers':2,'protocol_sha256':sha(ROOT/'protocol.json'),'points':queue})
    save(ROOT/'audit/preparation.json',{'passed':True,'source_hashes_matched':4,'PAW_hashes_verified':2,
        'formula':'C4H8','natoms':12,'Pnam_number':62,'quantum_calculations_started':False})
    for folder in ['runs','tables','control']: (ROOT/folder).mkdir(exist_ok=True)
    print(json.dumps({'status':'prepared','samples':len(samples),'protocol_sha256':sha(ROOT/'protocol.json')}))

if __name__=='__main__':main()
