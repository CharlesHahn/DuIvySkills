# RMSF

纯 Python 模块，计算均方根涨落。使用 MDAnalysis 原子选择语法。

## YAML 配置

```yaml
- RMSF:
    calc_group: protein
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `calc_group` | 是 | 计算组，MDAnalysis 语法 |
| `mkdir` | 否 | 输出目录名 |
| `frame_start` | 否 | 起始帧 |
| `frame_end` | 否 | 结束帧 |
| `frame_step` | 否 | 帧步长 |

## 输出

- RMSF 数据文件（.xvg）
- 残基级 RMSF 图
- 原子级 RMSF 图

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/RMSF.html