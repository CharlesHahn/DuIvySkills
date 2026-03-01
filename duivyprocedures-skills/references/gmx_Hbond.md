# gmx_Hbond

依赖 GROMACS 计算氢键。

## YAML 配置

```yaml
- gmx_Hbond:
    group1: Protein
    group2: Ligands
    top2show: 6
    only_calc_number: no
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `group1` | 是 | 第一组原子名（**不能以数字开头**） |
| `group2` | 是 | 第二组原子名 |
| `top2show` | 否 | 展示前 N 个占有率最高的氢键，默认 6 |
| `only_calc_number` | 否 | 仅计算数量，不计算距离角度，默认 `no` |
| `mkdir` | 否 | 输出目录名 |
| `gmx_parm` | 否 | 可用参数：`a`, `r`, `da`, `b`, `e`, `dt` |

## 输出

- 氢键数量图（hbnum.png）
- 氢键占有率图（hbmap.png）
- 最高占有率氢键的距离/角度图
- CSV 文件：占有率、平均距离、平均角度

## 注意

组名不能以数字开头，如 `1ZIN` 应改为 `ZIN1` 或 `Ligands`。

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/gmx_Hbond.html