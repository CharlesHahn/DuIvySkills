# gmx_DSSP

依赖 GROMACS 计算蛋白质二级结构。

## YAML 配置

```yaml
- gmx_DSSP:
    group: Protein
    gmx_parm:
      tu: ns
      dt: 0.5
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `group` | 是 | 蛋白质组名，需包含骨架原子 |
| `mkdir` | 否 | 输出目录名 |
| `gmx_parm` | 否 | 传递给 `gmx do_dssp` 的参数 |

## 输出

- 二级结构热图
- 二级结构含量时间变化图
- 各残基二级结构分布图

## DSSP 安装

- GROMACS 2022 及以下：需安装 DSSP 3.0 并设置环境变量
- GROMACS 2023 及以上：无需额外安装

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/gmx_DSSP.html