# 原始候选到 DFT 最终分组的对应关系

只有四个代表结构实际进行本轮 DFT 精修。其余编号经已有 GFN2 分组继承代表的最终组号，不表示已对全部十五个候选分别做过 DFT 优化。

|原始候选|初始角度 (°)|GFN 组|本轮 DFT 输入代表|最终组|是否亲自完成 DFT 精修|
|---|---:|---|---|---|---|
|cand_0001|0||||否；GFN 不稳定控制点|
|cand_0009|5|group_001|cand_0025|dft_group_002|否；继承代表映射|
|cand_0017|10|group_001|cand_0025|dft_group_002|否；继承代表映射|
|cand_0025|15|group_001|cand_0025|dft_group_002|是|
|cand_0033|20|group_001|cand_0025|dft_group_002|否；继承代表映射|
|cand_0041|25|group_003|cand_0041|dft_group_001|是|
|cand_0049|30|group_001|cand_0025|dft_group_002|否；继承代表映射|
|cand_0057|35|group_001|cand_0025|dft_group_002|否；继承代表映射|
|cand_0065|40|group_001|cand_0025|dft_group_002|否；继承代表映射|
|cand_0073|45|group_001|cand_0025|dft_group_002|否；继承代表映射|
|cand_0081|50|group_004|cand_0081|dft_group_001|是|
|cand_0089|55|group_001|cand_0025|dft_group_002|否；继承代表映射|
|cand_0097|60|group_001|cand_0025|dft_group_002|否；继承代表映射|
|cand_0105|65|group_001|cand_0025|dft_group_002|否；继承代表映射|
|cand_0113|70|group_002|cand_0113|dft_group_001|是|
|cand_0121|75|group_001|cand_0025|dft_group_002|否；继承代表映射|

最终 dft_group_001 = 原 g2、g3、g4，代表为原 g4（cand_0081）。最终 dft_group_002 = 原 g1，代表为 cand_0025。两组是固定晶胞 Pnam 约束和既定几何/能量容差下的分组，不是两种已确认晶型。

cand_0001（0°）保留为历史 GFN 不稳定驻点记录，未进入本轮 DFT 精修。完整源文件 SHA256 见 tables/candidate_to_final_group_provenance.json。
