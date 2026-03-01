# Density

计算质量和电荷密度沿 X、Y、Z 方向的分布。

## YAML 配置

```yaml
- Density:
    atom_selection: protein
```

## 参数

| 参数 | 必需 | 说明 |
|------|------|------|
| `atom_selection` | 是 | 原子选择，MDAnalysis 语法 |
| `mkdir` | 否 | 输出目录名 |
| `frame_start` | 否 | 起始帧 |
| `frame_end` | 否 | 结束帧 |
| `frame_step` | 否 | 帧步长 |

## 输出

- 质量密度分布（X/Y/Z 方向）
- 电荷密度分布（X/Y/Z 方向）

## 官方文档

https://duivyprocedures-docs.readthedocs.io/en/latest/Density.html