"""Independent analytic closure record; no refined coordinates or quantum calls.

Writes only the conformation audit JSON and Markdown beside this file.
"""
from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import gemmi
import numpy as np
import spglib
from ase.geometry import get_dihedrals

HERE = Path(__file__).resolve().parent
CELL = np.diag([7.432, 4.929, 2.545])
CC = 1.54
CH = 1.09


def main():
    translation = CELL[2]
    transverse = math.sqrt(CC**2 - (CELL[2, 2] / 2)**2)
    theta = math.radians(35)
    e = np.array([math.cos(theta), math.sin(theta), 0.0])
    a = 0.5 * transverse * e + 0.25 * translation
    b = -0.5 * transverse * e + 0.75 * translation
    sequence = np.array([a, b, a + translation, b + translation, a + 2 * translation])
    bonds = np.diff(sequence, axis=0)
    dihedrals = [
        float(get_dihedrals(bonds[k][None], bonds[k+1][None], bonds[k+2][None])[0])
        for k in (0, 1)
    ]
    assert np.allclose(dihedrals, 180.0, atol=1e-10)
    assert np.allclose(np.linalg.norm(bonds, axis=1), CC, atol=1e-12)
    assert np.allclose(sequence[2:] - sequence[:-2], translation, atol=1e-12)

    states = {"t": 180.0, "g+": 60.0, "g-": -60.0}
    records = []
    for labels in itertools.product(states, repeat=2):
        desired = [states[label] for label in labels]
        accepted = labels == ("t", "t")
        records.append({
            "torsion_labels": list(labels),
            "requested_torsions_deg": desired,
            "periodic_closure_compatible": accepted,
            "decision": "accepted_analytic_with_geometric_check" if accepted
                        else "analytic closure rejection",
            "rejection_basis": None if accepted else
                "Two-carbon pure-translation periodicity requires both nondegenerate "
                "consecutive dihedrals to equal 180 degrees modulo 360; at least one "
                "requested gauche dihedral contradicts this necessary condition.",
            "numerical_construction_attempted": accepted,
            "construction_failure_claimed": False,
            "geometric_torsions_deg": dihedrals if accepted else None,
            "quantum_calculation_performed": False,
        })

    gemmi_ops = list(gemmi.SpaceGroup("P n a m").operations())
    database = spglib.get_symmetry_from_database(297)

    def key(r, t):
        return (tuple(np.rint(r).astype(int).ravel()),
                tuple(np.round(np.asarray(t) % 1, 12)))

    expected = {key(r, t) for r, t in
                zip(database["rotations"], database["translations"])}
    actual = {key(np.asarray(op.rot) / op.DEN, np.asarray(op.tran) / op.DEN)
              for op in gemmi_ops}
    assert expected == actual
    archive = {
        "scope": "All nine ordered pairs over {t,g+,g-} for one nondegenerate "
                 "two-carbon chain under pure translation T=(0,0,c). This is not "
                 "an enumeration of longer-period, helical, defective or disordered chains.",
        "experimental_refined_coordinates_read": False,
        "quantum_calculations_performed": False,
        "inputs": {"cell_A": CELL.tolist(), "CC_A": CC, "CH_A": CH,
                   "spacegroup": "Pnam", "hall_number": 297, "hall_symbol": "-P 2c 2n",
                   "number_of_chains": 2, "CH2_per_cell": 4},
        "torsion_convention": {"t_deg": 180.0, "g_plus_deg": 60.0, "g_minus_deg": -60.0,
                              "values_equivalent_modulo_360": True},
        "enumerated_count": len(records),
        "closure_compatible_count": sum(r["periodic_closure_compatible"] for r in records),
        "analytic_rejection_count": sum(not r["periodic_closure_compatible"] for r in records),
        "records": records,
        "proof": [
            "Two-carbon pure-translation periodicity makes consecutive carbon positions "
            "A, B, A+T, B+T, A+2T, ... . Set v=B-A and w=T-v.",
            "Consecutive bond vectors are v,w,v,w,... . The normals for the first "
            "dihedral are v cross w and w cross v = -(v cross w); the next "
            "dihedral has the same opposite-normal relation. If v cross w != 0, "
            "both dihedrals are trans, 180 degrees modulo 360.",
            "Equal C-C bond lengths give |v|^2=|T-v|^2, hence v dot T=|T|^2/2. "
            "For |T|=c and bond length l, each bond has axial projection c/2 and "
            "transverse magnitude sqrt(l^2-c^2/4). Here l>c/2, so the chain is "
            "noncollinear and the dihedrals are defined.",
            "A gauche request is rejected by this analytic necessary closure "
            "condition. No failed numerical builder or optimization is asserted.",
        ],
        "geometric_check": {
            "constructed_from": "Chemical bond length, cell translation, arbitrary "
                                "setting angle 35 degrees; no coordinate template",
            "carbon_positions_A": sequence.tolist(),
            "torsions_deg": dihedrals,
            "CC_lengths_A": np.linalg.norm(bonds, axis=1).tolist(),
            "two_step_translation_error_A": float(
                np.linalg.norm(sequence[2:] - sequence[:-2] - translation, axis=1).max()),
            "CCC_deg": float(np.degrees(2 * np.arcsin(CELL[2, 2] / (2 * CC)))),
            "zigzag_half_width_A": transverse / 2,
        },
        "symmetry": {
            "gemmi_spglib_operation_sets_identical": expected == actual,
            "operations_xyz": [op.triplet() for op in gemmi_ops],
            "fixed_origin_convention": "Pnam Hall 297, -P 2c 2n",
            "chain_centers_fractional": [
                {"first": [x, y], "second": [(x + 0.5) % 1, (0.5 - y) % 1]}
                for x, y in itertools.product((0.0, 0.5), repeat=2)
            ],
            "carbon_phase_fractional": [0.25, 0.75],
            "mirror_argument": "The operation (x,y,1/2-z) fixes the transverse "
                "coordinates. The two carbons on a nonzero-width zigzag have distinct "
                "transverse positions. Thus each carbon lies on a mirror plane: "
                "z=1/4 or 3/4 modulo 1. This imposes no condition on its x,y coordinates "
                "and therefore no condition on the in-plane setting angle.",
            "axis_argument": "The operation G:(x,y,z)->(x+1/2,1/2-y,z) exchanges the "
                "two chain axes, so Baxis=(X+1/2,1/2-Y). The screw S:(x,y,z)->"
                "(-x,-y,z+1/2) cannot exchange these two axes, because its y equation "
                "would require -Y=1/2-Y modulo 1. It therefore fixes each axis, "
                "requiring 2X=2Y=0 modulo 1.",
            "atomic_sites": "Each initial C/H independent atom is on a 4c mirror "
                "site (x,y,1/4); x and y are free. Fixed chain axes must not be "
                "misdescribed as fixed atomic coordinates.",
            "setting_angle": "Continuous angle of the whole trans zigzag plane "
                "around c; the second chain is symmetry related. It is not a "
                "C-C-C-C torsion.",
            "origin_preserving_global_translations": [
                list(x) for x in itertools.product((0.0, 0.5), repeat=3)
            ],
            "general_origin_change_rule": "For f'=f+delta, keep the same physical "
                "structure by transforming each translation to t'=t+(I-R)delta "
                "modulo integers. An arbitrary origin change is a gauge choice, "
                "not an independent packing or axial-slide degree of freedom.",
            "independent_relative_chain_slide_allowed_in_this_model": False,
            "independent_relative_chain_rotation_allowed_in_this_model": False,
            "linear_target_group_displacement_dimension": 6,
            "rigid_initial_geometry_continuous_parameter_count": 1,
        },
        "generation_reconciliation": {
            "setting_angles_deg": list(range(0, 180, 5)),
            "centers_count": 4, "phases_count": 2, "grid_count": 288,
            "excluded_endpoint_180": "Same physical chain-plane family as 0; "
                                      "verified separately by periodic atom matching.",
            "expected_geometric_survivors": 248,
            "expected_distinct_initial_structures": 16,
            "interpretation": "These counts concern the declared rigid-chain grid "
                              "after analytic torsion closure; not an exhaustive "
                              "search of continuous molecular crystal structures.",
        },
    }
    destination = HERE / "conformation_enumeration.json"
    destination.write_text(json.dumps(archive, ensure_ascii=False, indent=2,
                                     allow_nan=False) + "\n", encoding="utf-8")
    note = """# PE 两碳周期构象与 Pnam 自由度独立审查

本档案只使用 C–C/C–H 化学几何、指定晶胞和公开空间群操作，不读取实验精修原子坐标，也不运行量化计算。

## 九个扭转组合的闭合结论

按 t=180°、g+=60°、g−=−60° 枚举两个连续二面角，共 9 个有序组合。只有 (t,t) 与两碳纯平移周期相容。其余 8 个明确记为 `analytic closure rejection`，表示解析闭合约束排除，不表示数值构建、优化或量化计算失败。

证明：两碳周期使连续原子必为 A、B、A+T、B+T、A+2T，其中 T=(0,0,c)。相邻键向量交替为 v=B−A 和 w=T−v。二面角两侧平面的法向量分别为 v×w 与 w×v=−v×w；非共线时两连续二面角均为 180°。等长键再给出 |v|²=|T−v|²，因此 v·T=|T|²/2。给定 c=2.545 Å、C–C=1.54 Å，存在非零横向分量，二面角定义良好。脚本由化学长度生成 5 个连续碳，实际验证两个二面角和两步平移；不为已解析排除的 gauche 组合虚构失败日志。

该证明不排除更长周期、螺旋、缺陷或无序 PE 的 gauche 构象；它仅限定本次两碳纯平移周期。

## 镜面、链轴与平移

采用 Pnam #62、Hall 297（−P 2c 2n）的固定原点。镜面操作 (x,y,1/2−z) 要求这些非退化平面链的原子处于 z=1/4 或 3/4。镜面不限制 x、y，所以绕 c 旋转整个锯齿平面的 setting angle 可以连续变化；setting angle 与主链二面角是不同变量。

令第一条链轴为 (X,Y)，操作 (x+1/2,1/2−y,z) 产生第二条链轴 (X+1/2,1/2−Y)。二重螺旋 (−x,−y,z+1/2) 不可能交换这两轴，因为这将要求 −Y=1/2−Y（模 1）；因此它固定每条链轴，要求 X、Y 各为 0 或 1/2。两条链不具有独立的相对连续滑移或转角自由度。

链轴位置固定不等于原子坐标固定：C、H、H 三个独立原子分别属于 4c (x,y,1/4) 镜面位点，具有总计 6 个线性坐标自由度。初始刚性化学链将这些自由度关联后，连续扫描变量为一个 setting angle。约束优化可以在这 6 个坐标自由度内改变键长和键角。

全晶体任意整体平移属于原点规范变化。若保持当前 Hall 操作不变，允许的整体平移为三方向的整数/半整数格点组合；若采用任意原点偏移 δ，则必须同步变换操作平移部分 t′=t+(I−R)δ。整体原点变化不应当作新的独立堆积。

## 与生成计数的关系

主生成网格取 0–175°、步长 5°，每角枚举 4 个链轴表示和 2 个碳相位，合计 36×4×2=288 条。180° 是 0° 的周期坐标等价端点。1–2、1–3 周期图邻接从非键排斥筛查中排除；其余原子按 0.85 倍范德华半径和及链间 C···C≥3.4 Å 检查。已有独立几何枚举和主程序均得到 248 个几何通过记录、16 个去重初始结构。θ=0° 可有比目标群更高的对称性；判断包含指定 Pnam 子群即可，不能因空间群识别出超群就误拒绝。

这些计数仅对已声明的有限网格完备。它们不是所有连续构象或所有可能 PE 晶体的穷举结论。
"""
    (HERE / "conformation_enumeration.md").write_text(note, encoding="utf-8")
    print(json.dumps({"file": str(destination), "enumerated": len(records),
                      "accepted": 1, "analytic_rejections": 8,
                      "torsions_deg": dihedrals, "no_quantum_calculations": True}))


if __name__ == "__main__":
    main()
