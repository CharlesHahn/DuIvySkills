# DuIvyTools 安装和配置

## 安装方法

DuIvyTools (DIT) 可以通过以下两种方式安装：

### 方法1：从源码安装

```bash
git clone https://github.com/CharlesHahn/DuIvyTools.git
cd DuIvyTools
pip install -e .
```

### 方法2：通过 pip 安装（推荐）

```bash
pip install DuIvyTools -i https://pypi.tuna.tsinghua.edu.cn/simple
```

或者使用官方 PyPI：

```bash
pip install DuIvyTools
```

## 验证安装

安装完成后，在命令行中输入：

```bash
dit
```

如果显示 DuIvyTools 的帮助信息和所有可用命令，说明安装成功。

## 绘图引擎配置

DuIvyTools 支持四种绘图引擎：

1. **matplotlib**（默认）：无需额外配置
2. **plotly**：需要安装 `pip install plotly`
3. **gnuplot**：需要自行安装 gnuplot 并确保在命令行中可以调用
4. **plotext**：需要安装 `pip install plotext`（仅用于在终端中绘制简单图形）

## 绘图样式自定义

### matplotlib 样式

在工作目录中创建 `dit_mplstyle.mplstyle` 文件，自定义 matplotlib 绘图样式。DIT 会自动加载该文件。

示例样式参数：

```python
axes.labelsize:     12
axes.linewidth:     1
xtick.labelsize:    12
ytick.labelsize:    12
lines.linewidth:    2
legend.fontsize:    12
font.family:        Arial
image.cmap:         coolwarm
figure.dpi:         100
savefig.dpi:        300
```

### plotly 样式

在工作目录中创建 `DIT_plotly.json` 文件，自定义 plotly 配置。

### gnuplot 样式

在工作目录中创建 `DIT_gnuplot.gp` 文件，自定义 gnuplot 配置。

## 系统要求

- Python 3.6+
- matplotlib（默认绘图引擎）
- numpy
- pandas

## 更新 DuIvyTools

```bash
pip install --upgrade DuIvyTools
```

## 卸载

```bash
pip uninstall DuIvyTools
```

## 引用

如果在研究中使用了 DuIvyTools，请通过以下 DOI 引用：

https://doi.org/10.5281/zenodo.6339993

## 获取帮助

- 查看所有命令：`dit`
- 查看全局参数：`dit -h` 或 `dit --help`
- 查看特定命令帮助：`dit <command> -h`
- 例如：`dit xvg_show -h`

## 相关资源

- GitHub 仓库：https://github.com/CharlesHahn/DuIvyTools
- 在线文档：https://duivytools.readthedocs.io/
- PyPI 页面：https://pypi.org/project/DuIvyTools/