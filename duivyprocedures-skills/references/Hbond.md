# Hbond

纯 Python 模块，氢键分析。使用 MDAnalysis 原子选择语法。

## YAML 配置

```yaml
- Hbond:
    donor_group: protein
    acceptor_group: protein
    update_selection: no
    d_h_cutoff: 0.12
    d_a_cutoff: 0.30
    d_h_a_angle_cutoff: 150
    only_calc_number: no
    top2show: 10
    calc_lifetime: no
    tau_max: 20
    window_step: 1
    intermittency: 0
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `donor_group` | 是 | 供体组，MDAnalysis 语法 |
| `acceptor_group` | 是 | 受体组，MDAnalysis 语法 |
| `update_selection` | 否 | 每帧刷新选择。涉及 `around` 等动态选择时设为 `yes`，否则 `no` 以提高速度 |
| `d_h_cutoff` | 否 | 供体-氢距离阈值（nm），默认 0.12。用于识别供体-氢原子对 |
| `d_a_cutoff` | 否 | 供体-受体距离阈值（nm），默认 0.30。小于此值形成氢键 |
| `d_h_a_angle_cutoff` | 否 | 供体-氢-受体角度阈值（度），默认 150。大于此值形成氢键 |
| `only_calc_number` | 否 | 仅计算氢键数量。预估氢键数特别大时（如蛋白-水）设为 `yes` |
| `top2show` | 否 | 展示前 N 个高占有率氢键，默认 10 |
| `calc_lifetime` | 否 | 计算生命周期（仅对 top2show 个氢键）。动态选择时不可用 |
| `tau_max` | 否 | 生命周期最大时间（帧），窗口越大积分越准确 |
| `window_step` | 否 | 生命周期窗口平移步长（帧） |
| `intermittency` | 否 | 允许帧间隔，默认 0 表示必须连续发生 |

## 输出

- 氢键数量图（hbnum.png）
- 氢键占有率图（hbmap.png）
- 高占有率氢键详情图（hbmap_top10.png）
- CSV 文件：占有率、出现帧数、平均距离、距离标准误、平均角度、角度标准误
- 生命周期自相关函数图及积分（如计算）

## 氢键命名格式

`残基名残基号原子(原子号)@氢原子...受体原子`

如：`GLY35N(312)@GLY35H(313)...ASP33OD2(301)`

## 生命周期计算注意

自相关函数未降到 0 时，需调大 `tau_max` 以获得准确积分。生命周期通过 simpson 积分计算，准确度一般。

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/Hbond.html
