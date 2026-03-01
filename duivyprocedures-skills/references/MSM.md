# MSM

马尔可夫状态模型构建。**此模块为 demo 性质，需谨慎用于研究。**

## 构建 MSM 的一般步骤

1. 对研究对象和 CV 进行充分采样
2. 选择合适的 features（坐标、二面角、二级结构等）
3. 对 features 进行 tICA 降维
4. 对降维数据进行聚类
5. 分析不同 lag 参数下的 implied timescales
6. 选取合适 lag 参数构建 MSM
7. 使用 ck_test 检验 MSM 有效性
8. 使用 PCCA 分析 meta states
9. 绘制状态转移网络图

## YAML 配置

```yaml
- MSM:
    atom_selection: protein and name CA
    byType: atom
    target: coordinates
    coordinate_fit: no
    dimension: 3
    lag4tica: 20
    number_clusters: 100
    transition_count_mode: "effective"
    lag4ITS_range: [1, 10, 20, 30, 50, 100, 200]
    lag4MSM: 20
    meta_number: 3
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `atom_selection` | 是 | 原子选择，MDAnalysis 语法。二面角分析需包含骨架原子 |
| `byType` | 否 | 计算方式：`atom`、`res_com`、`res_cog`、`res_coc`（仅 `coordinates` 有效） |
| `target` | 否 | 目标：`coordinates` 或 `dihedrals` |
| `coordinate_fit` | 否 | 是否对坐标进行点云 fit |
| `dimension` | 否 | tICA 降维后维度 |
| `lag4tica` | 否 | tICA 的 lag 参数 |
| `number_clusters` | 否 | 聚类子状态数量 |
| `transition_count_mode` | 否 | 转移矩阵计算方法：`sliding`、`sample`、`effective`、`sliding-effective` |
| `lag4ITS_range` | 否 | ITS 计算的 lag 参数列表 |
| `lag4MSM` | 否 | 构建 MSM 的 lag 参数 |
| `meta_number` | 否 | meta states 数量 |

## 二面角分析注意

二面角具有周期性，此模块将角度转换为 sin/cos 值再分析。请对照文献确认计算过程是否合适。

## 输出

- tICA 散点分布图
- 聚类结果图
- implied timescales 图
- timescales 图
- ck_test 图（预测线与测量线应重合）
- PCCA 结果（屏显：transition matrix, stationary probability, MFPT 等）
- 状态转移网络图

## 有效性检验

- `state_fraction` 和 `count_fraction` 应为 1.0
- ck_test 图上预测线与测量线应重合
- ITS 图上彩色线应平稳

**参数不合适时会报错终止，需调整参数重试。**

## 参考资料

- http://www.emma-project.org/latest/tutorials/notebooks/00-pentapeptide-showcase.html
- https://deeptime-ml.github.io/trunk/index.html

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/MSM.html
