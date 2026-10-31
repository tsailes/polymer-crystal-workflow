"""Read-only numerical diagnostics of saved PA6 forces; NO quantum calls.

Writes only symmetry_* and sampling_* in this diagnosis folder.
"""
from __future__ import annotations
import hashlib
import itertools
import json
from pathlib import Path
import sys

import numpy as np
from ase.constraints import FixSymmetry
from ase.io import read
from ase.spacegroup.symmetrize import symmetrize_rank1
from ase.units import GPa
import spglib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROJECT = ROOT.parents[1]
sys.path.insert(0, str(PROJECT / 'scripts'))
import validate_ordered as validator

def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def save(name, obj):
    def encode(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        raise TypeError(type(x).__name__)
    (HERE / name).write_text(json.dumps(obj, indent=2, default=encode, allow_nan=False)+'\n', encoding='utf-8')

def stats(x):
    x = np.asarray(x)
    return dict(min=float(x.min()), max=float(x.max()), mean=float(x.mean()), rms=float(np.sqrt(np.mean(x*x))))

record = load(ROOT / 'input/candidate.json')
initial = read(ROOT / 'input/candidate.extxyz')
final = read(ROOT / 'relaxation/final.extxyz')
roles = {r['index']:r['role'] for r in record['atom_roles']}
chains = {r['index']:r['chain_id'] for r in record['atom_roles']}
N = len(final)

# Match exactly the actual ASE constraint constructor, using a disposable copy.
working = initial.copy()
constraint = FixSymmetry(working, symprec=1e-5, adjust_cell=False)
def project(f):
    a = np.asarray(f, dtype=float).copy()
    constraint.adjust_forces(working, a)
    return a

P = np.empty((3*N, 3*N))
for col in range(3*N):
    unit = np.zeros((N,3)); unit.flat[col] = 1.
    P[:, col] = project(unit).ravel()
eig = np.linalg.eigvalsh((P+P.T)/2)
actual_ops = list(zip(constraint.rotations, constraint.translations))
def equivalent_op(a,b):
    return bool(np.array_equal(a[0],b[0]) and np.max(np.abs((a[1]-b[1])-np.rint(a[1]-b[1]))) < 1e-8)
declared_indices=[next(i for i,a in enumerate(actual_ops) if equivalent_op(a,b)) for b in validator.OPS]
def declared_project(f):
    return symmetrize_rank1(working.cell.array, working.cell.reciprocal().T, np.asarray(f,dtype=float),
        [constraint.rotations[i] for i in declared_indices],
        [constraint.translations[i] for i in declared_indices],
        [constraint.symm_map[i] for i in declared_indices])
PD=np.empty_like(P)
for col in range(3*N):
    unit=np.zeros((N,3));unit.flat[col]=1.
    PD[:,col]=declared_project(unit).ravel()
eig_declared=np.linalg.eigvalsh((PD+PD.T)/2)
dataset = spglib.get_symmetry_dataset((working.cell.array, working.get_scaled_positions(), working.numbers), symprec=1e-5)
raw_saved = np.asarray(load(ROOT/'relaxation/final_raw_results.json')['forces_eV_A'])
trajectory = read(ROOT/'relaxation/trajectory.traj', ':')
trace = []
for step, atoms in enumerate(trajectory):
    # Saved SinglePointCalculator only: never invoke a new energy/force engine.
    f = np.asarray(atoms.calc.results['forces'])
    projected = project(f)
    trace.append(dict(step=step, raw_fmax=float(np.linalg.norm(f, axis=1).max()),
                      projected_fmax=float(np.linalg.norm(projected, axis=1).max()),
                      max_removed_force_vector_eV_A=float(np.linalg.norm(f-projected, axis=1).max())))

# Rigid chain displacement basis, after subtracting common translation.
rigid = []
for chain in sorted(set(chains.values())):
    ids = [i for i in range(N) if chains[i]==chain]
    for axis in range(3):
        d = np.zeros((N,3)); d[ids,axis]=1.; d-=d.mean(axis=0)
        rigid.append(d.ravel())
B = np.asarray(rigid).T
symmetry = dict(
    scope='Saved-data and ASE geometry/projection analysis only; zero quantum calculations',
    actual_space_group_number=int(dataset.number), actual_international=str(dataset.international),
    actual_hall_number=int(dataset.hall_number), operation_count=len(actual_ops),
    declared_operation_count=len(validator.OPS),
    declared_operations_are_subset=bool(len(declared_indices)==len(validator.OPS)),
    actual_and_declared_operations_equal=bool(len(actual_ops)==len(validator.OPS) and all(any(equivalent_op(a,b) for b in validator.OPS) for a in actual_ops)),
    actual_operations=[dict(rotation=r, translation=t) for r,t in actual_ops],
    additional_operation_indices=[i for i in range(len(actual_ops)) if i not in declared_indices],
    extra_translation=[0,.5,0],
    extra_translation_interpretation='Accidental b/2 translation maps the starting structure to itself; retaining all detected operations also preserves this smaller periodicity. No.14 alone does not describe which translation subgroup was constrained.',
    constraint_initialization_position_max_A=float(np.abs(initial.positions-working.positions).max()),
    constraint_initialization_cell_max_A=float(np.abs(initial.cell.array-working.cell.array).max()),
    position_dofs_unconstrained=3*N, projected_rank=int(np.count_nonzero(eig>.5)),
    excluded_directions=int(np.count_nonzero(eig<.5)),
    declared_four_operations_projected_rank=int(np.count_nonzero(eig_declared>.5)),
    extra_suppressed_directions_due_to_8_vs_4=int(np.count_nonzero(eig_declared>.5)-np.count_nonzero(eig>.5)),
    declared_projection_idempotence_max=float(np.abs(PD@PD-PD).max()),
    projection_idempotence_max=float(np.abs(P@P-P).max()), projection_asymmetry_max=float(np.abs(P-P.T).max()),
    projection_eigenvalue_extrema=[float(eig.min()),float(eig.max())],
    relative_rigid_chain_translation_dimension=int(np.linalg.matrix_rank(B,tol=1e-8)),
    permitted_rigid_chain_translation_dimension=int(np.linalg.matrix_rank(P@B,tol=1e-8)),
    excluded_rigid_chain_translation_dimension=int(np.linalg.matrix_rank((np.eye(3*N)-P)@B,tol=1e-8)),
    declared_four_operations_permitted_rigid_chain_translation_dimension=int(np.linalg.matrix_rank(PD@B,tol=1e-8)),
    additional_rigid_translation_dimension_released_by_8_to_4=int(np.linalg.matrix_rank((PD-P)@B,tol=1e-8)),
    final_saved_force_max_removed_eV_A=float(np.linalg.norm(raw_saved-project(raw_saved),axis=1).max()),
    trajectory_force_projection=trace,
    max_removed_force_vector_over_trajectory_eV_A=max(t['max_removed_force_vector_eV_A'] for t in trace),
    fixed_cell_max_initial_final_change_A=float(np.abs(initial.cell.array-final.cell.array).max()),
    source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [
        Path(__file__),ROOT/'input/candidate.extxyz',ROOT/'input/candidate.json',ROOT/'relaxation/final.extxyz',
        ROOT/'relaxation/trajectory.traj',ROOT/'relaxation/final_raw_results.json',Path(validator.__file__)]},
    interpretation=[
        'Raw and projected force equality excludes a large instantaneous force discarded by the group projector at these symmetric geometries.',
        'It does not exclude negative curvature, bifurcation or lower minima along symmetry-breaking modes. Such first derivatives can vanish identically by symmetry.',
        'FixSymmetry acts on individual atom displacements; it is not a rigid-chain constraint. Chain geometry is allowed to relax subject to symmetry relations.',
        'The fixed lattice excludes strain relaxation. Retained primitive translations exclude long-period distortions. Neither is tested by raw/projected force agreement.'
    ])
save('symmetry_diagnosis.json',symmetry)

base=load(ROOT/'final_cold_tight/result.json')
F0=np.asarray(base['forces'])
groups={'all':np.arange(N), 'amide_N':np.array([i for i in roles if roles[i]=='amide_N']),
        'carbonyl_O':np.array([i for i in roles if roles[i]=='carbonyl_O']),
        'amide_H':np.array([i for i in roles if roles[i]=='amide_H']),
        'carbonyl_C':np.array([i for i in roles if roles[i]=='carbonyl_C']),
        'methylene_C':np.array([i for i in roles if roles[i].startswith('methylene_C')]),
        'methylene_H':np.array([i for i in roles if roles[i] not in ('amide_N','amide_H','carbonyl_O','carbonyl_C') and final.numbers[i]==1])}

# Each donor's nearest cross-chain acceptor by HO or NO. Enumerating [-2,2]^3
# around each N donor is ample here; check selected vectors against exact MIC.
donors=[]
for i in groups['amide_N']:
    bond=next(b for b in record['bonds'] if (b['i']==i and roles[b['j']]=='amide_H') or (b['j']==i and roles[b['i']]=='amide_H'))
    h=bond['j'] if bond['i']==i else bond['i']
    hs=np.array(bond['shift'])*(1 if bond['i']==i else -1)
    pn=final.positions[i]; ph=final.positions[h]+hs@final.cell.array
    pairs=[]
    for o in groups['carbonyl_O']:
        if chains[i]==chains[o]: continue
        for shift in itertools.product(range(-2,3),repeat=3):
            po=final.positions[o]+np.asarray(shift)@final.cell.array
            ho=po-ph; no=po-pn; nh=ph-pn
            angle=np.degrees(np.arccos(np.clip((-nh)@ho/np.linalg.norm(nh)/np.linalg.norm(ho),-1,1)))
            pairs.append(dict(N=int(i),H=int(h),O=int(o),donor_NH_image_shift=hs,
                              oxygen_image_from_N=list(shift),N_O_A=float(np.linalg.norm(no)),
                              H_O_A=float(np.linalg.norm(ho)),N_H_O_deg=float(angle),
                              H_to_O_unit=ho/np.linalg.norm(ho),N_to_O_unit=no/np.linalg.norm(no)))
    for criterion in ['H_O_A','N_O_A']:
        pair=min(pairs,key=lambda p:p[criterion]); pair=dict(pair,selection=criterion)
        donors.append(pair)

forces={'1x1x1':F0}; points=[]
for label in ['1x1x1','2x1x1','1x2x1','1x1x2']:
    if label=='1x1x1':
        r=base; f=F0; replicas=1; deviation=0.
    else:
        r=load(ROOT/f'final_sampling_{label}/result.json')
        replicas=int(np.prod(list(map(int,label.split('x')))))
        ff=np.asarray(r['forces']).reshape(replicas,N,3)
        f=ff.mean(axis=0); deviation=float(np.linalg.norm(ff-f,axis=2).max())
        forces[label]=f
    delta=f-F0
    role_summary={}
    for name,ids in groups.items():
        norms=np.linalg.norm(f[ids],axis=1); dn=np.linalg.norm(delta[ids],axis=1)
        role_summary[name]=dict(atom_count=len(ids),force_norm_eV_A=stats(norms),delta_force_norm_eV_A=stats(dn),
                              max_force_atom=int(ids[norms.argmax()]),max_delta_atom=int(ids[dn.argmax()]))
    points.append(dict(repeat=label,replicas=replicas,
                       energy_delta_meV_per_PA6_repeat=(r['energy']/replicas-base['energy'])*1000/8,
                       fmax_eV_A=float(np.linalg.norm(f,axis=1).max()),
                       removed_symmetry_force_eV_A=float(np.linalg.norm(f-project(f),axis=1).max()),
                       replica_force_max_deviation_eV_A=deviation,
                       stress_voigt_GPa=np.asarray(r['stress'])/GPa,
                       hydrostatic_pressure_GPa=float(-np.asarray(r['stress'])[:3].mean()/GPa),
                       role_summary=role_summary))

contact_projections=[]
for pair in donors:
    n,h,o=pair['N'],pair['H'],pair['O']; u=np.asarray(pair['H_to_O_unit']); uno=np.asarray(pair['N_to_O_unit'])
    projections={}
    for label,f in forces.items():
        projections[label]=dict(H_toward_O_eV_A=float(f[h]@u), O_toward_H_eV_A=float(-f[o]@u),
            HO_generalized_closing_force_eV_A=float((f[h]-f[o])@u/2),
            HO_delta_generalized_closing_force_vs_base_eV_A=float(((f[h]-f[o])-(F0[h]-F0[o]))@u/2),
            NO_generalized_closing_force_eV_A=float((f[n]-f[o])@uno/2))
    contact_projections.append(dict(**pair,force_projections=projections))

geometry_reference_path=HERE/'geometry_diagnosis.json'
crosscheck=[]
if geometry_reference_path.exists():
    reference=load(geometry_reference_path)['frames']['final']['donors']
    for pair in contact_projections:
        donor=next(d for d in reference if d['N']==pair['N'])
        other=donor['nearest_by_HO' if pair['selection']=='H_O_A' else 'nearest_by_NO']
        crosscheck.append(dict(N=pair['N'],selection=pair['selection'],
             same_indices=bool(pair['H']==other['H'] and pair['O']==other['O']),
             same_oxygen_shift=bool(np.array_equal(pair['oxygen_image_from_N'],other['oxygen_image_from_N'])),
             max_distance_difference_A=max(abs(pair[k]-other[k]) for k in ['N_O_A','H_O_A']),
             angle_difference_deg=abs(pair['N_H_O_deg']-other['N_H_O_deg'])))

sampling=dict(scope='Same final geometry, fixed cell, accuracy .01; saved single points only; zero new quantum calls',
              reference='final_cold_tight/result.json', points=points,
              contact_force_projection_definition={
                'u':'unit vector H to selected periodic O image',
                'HO_generalized_closing_force':'(F_H-F_O) dot u / 2, for H += q*u/2 and O -= q*u/2; q>0 reduces H-O separation',
                'sign':'Positive favors approach along this artificial pair coordinate; negative favors separation.',
                'limitation':'Current residual-force projection, NOT a measured hydrogen-bond attraction, isolated pair interaction, curvature or restoring force. Moving H alone distorts N-H; this is a diagnostic projection only.'},
              contacts=contact_projections,
              independent_geometry_agent_crosscheck=crosscheck,
              source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),
                  ROOT/'final_cold_tight/result.json',*[ROOT/f'final_sampling_{s}/result.json' for s in ['2x1x1','1x2x1','1x1x2']]]},
              convergence_limitations=['Doubling one direction at a time is sensitivity, not combined-grid convergence.',
                  'Gamma replicates retain identical primitive geometry; no longer-period relaxation was allowed.',
                  'Coarse .05 eV/A target cannot be replaced by tighter production force accuracy without reoptimization.'])
save('sampling_diagnosis.json',sampling)
print(json.dumps({'symmetry':{k:v for k,v in symmetry.items() if k not in ['trajectory_force_projection','actual_operations','interpretation']},
                  'sampling_points':points},default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x,indent=2))
