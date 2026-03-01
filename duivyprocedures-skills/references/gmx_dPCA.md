# gmx_dPCA

依赖 GROMACS 进行蛋白质骨架二面角主成分分析。

## YAML 配置

```yaml
- gmx_dPCA:
    group: Protein
    fast_mode: no
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `group` | 是 | 蛋白质组名，需包含骨架原子 |
| `fast_mode` | 否 | 快速模式，`yes`/`no`。设为 `yes` 可跳过耗时统计计算 |
| `mkdir` | 否 | 输出目录名 |
| `gmx_parm` | 否 | 传递给 `gmx anaeig` 的参数 |

## 输出

- 前 3 主成分两两散点图
- 主成分占比图
- 主成分数据文件（可用于 `gmx_FEL`）
- 余弦含量

## 注意

此模块使用了一些技巧，不能保证所有情况下都能成功执行。建议仔细检查结果。

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/gmx_dPCA.html