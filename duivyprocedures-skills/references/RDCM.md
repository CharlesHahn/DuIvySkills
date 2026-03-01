# RDCM

残基距离接触矩阵分析。可基于 RDCM 计算 RMSD、RMSF、DCCM、PCA、聚类等衍生分析。

## YAML 配置

```yaml
- RDCM:
    type_select: center
    atom_selection: protein
    frames_output_step: -1
    calc_RMSD: no
    RMSD_Matrix_step: -1
    calc_RMSF: no
    calc_DCCM: no
    Pearson_Observe: ""
    calc_PCA: no
    clustering_step: -1
    calc_contact: yes
    contact_cutoff: 1.5
    calc_encounter: yes
    encounter_low_cutoff: 0.8
    encounter_high_cutoff: 1.0
    calc_encounter_DCCM: no
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `type_select` | 否 | 距离类型：`center`/`res_com`（质心）、`min`/`res_min`（最小距离）、`atom`（原子坐标）、`res_cog`（几何中心）、`res_coc`（电荷中心）。**注意：`min` 计算非常慢** |
| `atom_selection` | 是 | 原子选择，MDAnalysis 语法。`center`/`min` 时请包含残基所有原子 |
| `frames_output_step` | 否 | 输出矩阵的帧步长，`-1` 不输出 |
| `calc_RMSD` | 否 | 是否基于 RDCM 计算 RMSD |
| `RMSD_Matrix_step` | 否 | RMSD 矩阵计算步长，`-1` 不计算。步长太小时耗时会明显增加 |
| `calc_RMSF` | 否 | 是否基于 RDCM 计算 RMSF |
| `calc_DCCM` | 否 | 是否基于 RDCM 计算 DCCM |
| `Pearson_Observe` | 否 | 自定义变量（xvg 文件）计算 Pearson 相关系数。变量维度需与帧数一致 |
| `calc_PCA` | 否 | 是否基于 RDCM 计算 PCA |
| `clustering_step` | 否 | 残基和帧聚类步长，`-1` 不聚类 |
| `calc_contact` | 否 | 是否计算 contact |
| `contact_cutoff` | 否 | contact 距离阈值（nm） |
| `calc_encounter` | 否 | 是否计算 encounter |
| `encounter_low_cutoff` | 否 | encounter 形成阈值（nm） |
| `encounter_high_cutoff` | 否 | encounter 断裂阈值（nm） |
| `calc_encounter_DCCM` | 否 | 是否基于 encounter 矩阵计算 DCCM |

## 帧选择参数

```yaml
frame_start: 1000
frame_end: 5001
frame_step: 10
```

**注意**：蛋白质较大且帧数较多时，建议设置帧选择参数减少计算量，否则可能内存不足。

## 输出

**基础输出**：
- RDCM 初始帧、结束帧、中间帧矩阵（.xpm + .csv + .png）
- 帧间差值矩阵
- RDCM 平均值和标准偏差图

**衍生分析输出**（如设置）：
- RMSD 曲线、RMSD 矩阵
- RMSF 曲线
- DCCM 矩阵
- Pearson 相关系数矩阵及 p_value 矩阵
- PCA 散点图
- 残基聚类和帧聚类树状图

**Contact/Encounter 输出**：
- Contact/Encounter 占有率矩阵
- 局部接触时间曲线
- Encounter 形成时间、平均时间长度、次数矩阵
- C50/C70/C90 无量纲数（屏显和 log）

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/RDCM.html
