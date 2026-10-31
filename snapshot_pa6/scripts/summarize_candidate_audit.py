"""Summarize geometry-only audit records; never modify candidates or run energies."""
from pathlib import Path
from collections import defaultdict, Counter
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'branches/O/audit_20261029'

def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    path = OUT / 'hbond/grid_records.json'
    records = read(path)
    phase = read(OUT / 'phase/phase_symmetry_summary.json')
    maps = read(OUT / 'phase/x_mirror_equivalences.json')
    assert len(records) == 1440
    assert phase['x_mirror_passed'] == phase['x_mirror_tests'] == 864
    assert all(m['equivalent'] and m['role_map_passed'] and
               m['periodic_covalent_graph_passed'] for m in maps)
    maximum = max(m['max_residual_A'] for m in maps)
    assert maximum < 1e-7
    groups = defaultdict(list)
    for record in records:
        p = record['parameters']
        key = (p['template'], round(min(p['x'], .5-p['x']), 3),
               p['rotation_deg'], round(p['phase_c']*96))
        groups[key].append(record)
    flags = ['original_full_HB_coverage',
             'original_builder_geometry_without_HB_passed',
             'independent_geometry_without_HB_passed']
    result = []
    for key, members in sorted(groups.items()):
        for flag in flags:
            assert len({m[flag] for m in members}) == 1, (key, flag)
        representative = min(members, key=lambda m: m['parameters']['x'])
        result.append({
            'group_key': list(key), 'representative_audit_id': representative['audit_id'],
            'parameters': representative['parameters'],
            'members': [m['audit_id'] for m in members],
            **{flag: representative[flag] for flag in flags},
            'donor_coverage': representative['donors_covered']})
    summary = {
        'scope': 'Audit index only. Known x <-> 0.5-x numerical equivalence of current templates; not proven globally minimal.',
        'no_candidates_modified': True, 'energy_or_preoptimization_runs': 0,
        'input_entries': len(records), 'known_relation_groups': len(result),
        'maximum_mapped_residual_A': maximum,
        'equivalence_tolerance_A': 1e-7,
        'counts_after_known_relation': {
            flag: sum(g[flag] for g in result) for flag in flags},
        'groups_missing_complete_HB_but_original_geometry_passed': sum(
            g['original_builder_geometry_without_HB_passed'] and not g['original_full_HB_coverage']
            for g in result),
        'phase_step_A': 17.46/96, 'relative_opposite_chain_phase_step_A': 2*17.46/96,
        'source_hashes': {'records': sha(path), 'script': sha(Path(__file__)),
                          'mapping_evidence': sha(OUT / 'phase/x_mirror_equivalences.json')},
        'limitations': [
            'No proof of future relaxation outcomes.',
            'Do not transfer x folding to templates lacking the internal screw relation.',
            'Not a reduction of 48 distinct phases at fixed other parameters.',
            'Original candidate and rejected-design provenance are preserved.']}
    assert len(result) == 864
    assert summary['counts_after_known_relation']['original_full_HB_coverage'] == 87
    (OUT / 'known_equivalence_groups.json').write_text(
        json.dumps({'summary': summary, 'groups': result}, indent=2)+'\n', encoding='utf-8')
    (OUT / 'combined_summary.json').write_text(
        json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
