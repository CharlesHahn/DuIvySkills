# DuIvyTools 命令参考文档

## 目录

- [全局参数](#全局参数)
- [XVG 命令](#xvg-命令)
- [XPM 命令](#xpm-命令)
- [NDX 命令](#ndx-命令)
- [其他命令](#其他命令)

---

## 全局参数

所有 DuIvyTools 命令都支持以下参数：

### 输入输出参数

- `-f INPUT [INPUT ...], --input INPUT [INPUT ...]`：指定输入文件（一个或多个）
- `-o OUTPUT, --output OUTPUT`：指定输出文件名

### 数据选择参数

- `-c COLUMNS [COLUMNS ...], --columns COLUMNS [COLUMNS ...]`：选择数据列索引（从0开始）
  - 格式：`1-7,10` 表示第1到6列和第10列
  - 多组文件：`1-7,10 0,1,4` 第一组文件选第1-6列和第10列，第二组选第0、1、4列
  - 注意：列索引从0开始，`-` 表示左闭右开区间
- `-b BEGIN, --begin BEGIN`：指定起始行索引（包含）
- `-e END, --end END`：指定结束行索引（不包含）
- `-dt DT, --dt DT`：指定步长，默认为1

### 绘图标签参数

- `-l LEGENDS [LEGENDS ...], --legends LEGENDS [LEGENDS ...]`：指定图例标签
  - 不显示图例：`-l "" "" ""`
  - 支持 LaTeX 语法：`-l "$nm^2$" "$\Delta G_{energy}$"`
- `-x XLABEL, --xlabel XLABEL`：指定 X 轴标签
- `-y YLABEL, --ylabel YLABEL`：指定 Y 轴标签
- `-z ZLABEL, --zlabel ZLABEL`：指定 Z 轴标签（用于颜色、3D图等）
- `-t TITLE, --title TITLE`：指定图片标题
  - 不显示标题：`-t ""`
  - 支持 LaTeX 语法

### 数据范围参数

- `-xmin XMIN, --xmin XMIN`：指定 X 轴最小值
- `-xmax XMAX, --xmax XMAX`：指定 X 轴最大值
- `-ymin YMIN, --ymin YMIN`：指定 Y 轴最小值
- `-ymax YMAX, --ymax YMAX`：指定 Y 轴最大值
- `-zmin ZMIN, --zmin ZMIN`：指定 Z 轴最小值（颜色值）
- `-zmax ZMAX, --zmax ZMAX`：指定 Z 轴最大值（颜色值）

### 数据变换参数

**⚠️ 警告**：使用数据变换参数时请务必谨慎：
- 确保变换适用于数据类型（时间、能量、距离等）
- 变换后必须更新对应的轴标签
- 确保变换因子正确，不会误导科学解释
- 默认情况下应保持原始单位，仅在用户明确要求时才转换

- `-xs XSHRINK, --xshrink XSHRINK`：X 轴数据缩放因子，默认1.0
- `-ys YSHRINK, --yshrink YSHRINK`：Y 轴数据缩放因子，默认1.0
- `-zs ZSHRINK, --zshrink ZSHRINK`：Z 轴数据缩放因子，默认1.0
- `-xp XPLUS, --xplus XPLUS`：X 轴数据偏移量，默认0.0
- `-yp YPLUS, --yplus YPLUS`：Y 轴数据偏移量，默认0.0
- `-zp ZPLUS, --zplus ZPLUS`：Z 轴数据偏移量，默认0.0

### 数据精度参数

- `--x_precision X_PRECISION`：X 轴数据显示精度（小数位数）
- `--y_precision Y_PRECISION`：Y 轴数据显示精度（小数位数）
- `--z_precision Z_PRECISION`：Z 轴数据显示精度（小数位数）

### 刻度参数

- `--x_numticks X_NUMTICKS`：X 轴刻度数量
- `--y_numticks Y_NUMTICKS`：Y 轴刻度数量
- `--z_numticks Z_NUMTICKS`：Z 轴刻度数量（仅 matplotlib 有效）

### 统计分析参数

- `-smv [{,CI,origin}], --showMV [{,CI,origin}]`：显示滑动平均值
  - 无参数：显示滑动平均值，原始数据作为背景
  - `origin`：同无参数
  - `CI`：显示置信区间作为背景
- `-ws WINDOWSIZE, --windowsize WINDOWSIZE`：滑动平均窗口大小，默认50
- `-cf CONFIDENCE, --confidence CONFIDENCE`：置信区间可信度，默认0.95

### 绘图控制参数

- `-ns, --noshow`：不显示图形（用于批量处理或生成脚本）
- `--alpha ALPHA`：图形元素透明度
- `-csv CSV, --csv CSV`：将数据导出到 CSV 文件
- `-eg {matplotlib,plotext,plotly,gnuplot}, --engine {matplotlib,plotext,plotly,gnuplot}`：指定绘图引擎
  - `matplotlib`：默认引擎，功能最全面
  - `plotly`：交互式图形，支持3D和轮廓图
  - `gnuplot`：生成 gnuplot 脚本，需自行安装 gnuplot
  - `plotext`：在终端中绘制简单图形
- `-cmap COLORMAP, --colormap COLORMAP`：指定颜色映射（matplotlib 和 plotly 有效）
- `--colorbar_location {None,left,top,bottom,right}`：颜色条位置（仅 matplotlib 有效）
- `--legend_location {inside,outside}`：图例位置（matplotlib 和 gnuplot 有效）
- `--legend_ncol LEGEND_NCOL`：图例列数，默认1（matplotlib 有效）

### 插值参数

- `-ip INTERPOLATION, --interpolation INTERPOLATION`：插值方法
  - matplotlib imshow：支持多种插值方法（如 bilinear, bicubic, nearest 等）
  - 其他引擎：`linear`, `cubic`, `quintic`
- `-ipf INTERPOLATION_FOLD, --interpolation_fold INTERPOLATION_FOLD`：插值倍数，默认10

### 其他参数

- `-m {withoutScatter,pcolormesh,3d,contour,AllAtoms,pdf,cdf}, --mode {withoutScatter,pcolormesh,3d,contour,AllAtoms,pdf,cdf}`：模式选择
  - `withoutScatter`：不显示散点图（xvg_box_compare）
  - `pcolormesh`：使用 pcolormesh 绘制（xpm_show）
  - `3d`：3D 模式（xpm_show）
  - `contour`：等高线模式（xpm_show）
  - `AllAtoms`：在所有原子中寻找中心（find_center）
  - `pdf`：概率密度函数（xvg_show_distribution）
  - `cdf`：累积密度函数（xvg_show_distribution）
- `-al ADDITIONAL_LIST [ADDITIONAL_LIST ...], --additional_list ADDITIONAL_LIST [ADDITIONAL_LIST ...]`：附加参数列表（命令特定）

---

## XVG 命令

### xvg_show

**功能**：快速显示一个或多个 XVG 文件中的所有数据。

**说明**：
- 默认第0列为 X 值，第1列及之后为 Y 值
- 自动解析图例、X/Y 轴标签
- 从 v0.6.0 起支持无 legend 数据列的可视化

**参数**：
- `-f`：输入 XVG 文件（可多个）
- `-c`：选择要显示的列（可选）
- `-b`, `-e`, `-dt`：数据范围选择
- 其他标准绘图参数

**示例**：
```bash
# 显示单个文件
dit xvg_show -f rmsd.xvg

# 显示多个文件
dit xvg_show -f rmsd.xvg gyrate.xvg

# 选择特定列并设置标签
dit xvg_show -f rmsd.xvg -c 1 -x "Time (ns)" -y "RMSD (nm)" -t "RMSD Analysis"
```

---

### xvg_compare

**功能**：比较多个 XVG 文件中的数据，支持灵活的列选择和滑动平均。

**说明**：
- 比 xvg_show 更灵活，推荐使用
- 每个输入文件可以选择不同的列
- 支持滑动平均和置信区间显示

**参数**：
- `-f`：输入 XVG 文件（可多个）
- `-c`：列索引选择（每组文件对应一个索引组）
- `-l`：图例标签（数量需与选择的列数一致）
- `-smv`：显示滑动平均值
- `-ws`：滑动平均窗口大小
- `-cf`：置信区间可信度
- `-csv`：导出数据到 CSV
- 其他标准绘图参数

**示例**：
```bash
# 比较两个文件的不同列
dit xvg_compare -f rmsd.xvg gyrate.xvg -c 1 1,2

# 设置图例和轴标签（注意：使用原始单位 ps）
dit xvg_compare -f rmsd.xvg gyrate.xvg -c 1 1,2 -l RMSD gyrate gyrate_X -x "Time(ps)" -y "(nm)" -t hhh

# 显示滑动平均值
dit xvg_compare -f energy.xvg -c 1,3 -l "LJ(SR)" "Coulomb(SR)" -smv

# 导出数据到 CSV
dit xvg_compare -f energy.xvg -c 1,3 -l "LJ(SR)" "Coulomb(SR)" -ns -csv data.csv

# 使用 plotly 引擎（单位转换示例 - 仅作技术演示）
dit xvg_compare -f energy.xvg -c 1,3 -l "LJ(SR)" "Coulomb(SR)" -xs 0.001 -x "Time(ns)" -smv -eg plotly

# 使用 gnuplot 引擎（单位转换示例 - 仅作技术演示）
dit xvg_compare -f energy.xvg -c 1,3 -l "LJ(SR)" "Coulomb(SR)" -xs 0.001 -x "Time(ns)" -smv -eg gnuplot
```

---

### xvg_ave

**功能**：计算 XVG 文件中每一列数据的平均值、标准偏差和标准误差。

**说明**：
- 用于计算平衡时期数据的平均值
- `-b` 指定的起始行包含在计算中
- `-e` 指定的结束行不包含在计算中

**参数**：
- `-f`：输入 XVG 文件
- `-b`：起始行索引（包含）
- `-e`：结束行索引（不包含）

**示例**：
```bash
# 计算平均值
dit xvg_ave -f rmsd.xvg -b 1000 -e 2001

# 输出示例：
# >>>>>>>>>                    rmsd.xvg                    <<<<<<<<<<<<<<
# ----------------------------------------------------------------------------
# |                  |     Average      |     Std.Dev      |     Std.Err      |
# ----------------------------------------------------------------------------
# |    Time (ps)     |   15000.000000   |   2891.081113    |    91.378334     |
# ----------------------------------------------------------------------------
# |    RMSD (nm)     |     0.388980     |     0.038187     |     0.001207     |
# ----------------------------------------------------------------------------
```

---

### xvg_show_distribution

**功能**：显示数据的分布情况，支持直方图、PDF 和 CDF。

**说明**：
- 默认显示数据列的分布（直方图）
- `-m pdf`：显示概率密度函数（Kernel Density Estimation）
- `-m cdf`：显示累积概率密度函数

**参数**：
- `-f`：输入 XVG 文件
- `-c`：选择要显示的列
- `-m`：显示模式（distribution, pdf, cdf）
- 其他标准绘图参数

**示例**：
```bash
# 显示分布
dit xvg_show_distribution -f gyrate.xvg -c 1,2

# 显示 PDF（概率密度函数）
dit xvg_show_distribution -f gyrate.xvg -c 1,2 -m pdf -eg plotly

# 显示 CDF（累积密度函数）
dit xvg_show_distribution -f gyrate.xvg -c 1,2 -m cdf -eg gnuplot
```

---

### xvg_show_stack

**功能**：绘制堆积折线图。

**说明**：
- 适用于绘制蛋白质二级结构含量变化等场景
- 自动将选中的列绘制为堆积折线图

**参数**：
- `-f`：输入 XVG 文件
- `-c`：选择要堆积的列
- 其他标准绘图参数

**示例**：
```bash
# 绘制二级结构含量堆积图（单位转换示例 - 仅作技术演示）
dit xvg_show_stack -f dssp_sc.xvg -c 2-7 -xs 0.001 -x "Time (ns)"
```

---

### xvg_show_scatter

**功能**：绘制散点图，支持两列或三列数据（第三列用于着色）。

**说明**：
- 两列：X-Y 散点图
- 三列：X-Y 散点图，第三列用于颜色映射

**参数**：
- `-f`：输入 XVG 文件
- `-c`：选择列（2列或3列）
- 其他标准绘图参数

**示例**：
```bash
# 两列散点图
dit xvg_show_scatter -f gyrate.xvg -c 1,2

# 三列散点图（用时间着色 - 单位转换示例，仅作技术演示）
dit xvg_show_scatter -f gyrate.xvg -c 1,2,0 -zs 0.001 -z "Time(ns)" -eg plotly --x_precision 2 --y_precision 2
```

---

### xvg_box_compare

**功能**：以小提琴图和散点图的形式比较数据列。

**说明**：
- 小提琴图显示数据分布
- 散点图显示原始数据点
- 可以使用 `-m withoutScatter` 隐藏散点图

**参数**：
- `-f`：输入 XVG 文件
- `-c`：选择要比较的列
- `-l`：图例标签
- `-z`：第三维度标签（用于着色）
- `-m withoutScatter`：不显示散点图
- 其他标准绘图参数

**示例**：
```bash
# 显示小提琴图和散点图（单位转换示例 - 仅作技术演示）
dit xvg_box_compare -f gyrate.xvg -c 1,2,3,4 -l Gyrate Gx Gy Gz -z "Time(ns)" -zs 0.001

# 仅显示小提琴图（单位转换示例 - 仅作技术演示）
dit xvg_box_compare -f gyrate.xvg -c 1,2,3,4 -l Gyrate Gx Gy Gz -z "Time(ns)" -zs 0.001 -m withoutScatter

# 使用 plotly 引擎（单位转换示例 - 仅作技术演示）
dit xvg_box_compare -f gyrate.xvg -c 1,2,3,4 -l Gyrate Gx Gy Gz -z "Time(ns)" -zs 0.001 -eg plotly

# 使用 gnuplot 引擎并设置 Y 轴范围（单位转换示例 - 仅作技术演示）
dit xvg_box_compare -f gyrate.xvg -c 1,2,3,4 -l Gyrate Gx Gy Gz -z "Time(ns)" -zs 0.001 -eg gnuplot -ymin 2
```

---

### xvg_combine

**功能**：从多个 XVG 文件中读取数据并组合成一个新的 XVG 文件。

**说明**：
- 与 xvg_compare 类似的输入方式
- 将选中的数据输出到新的 XVG 文件

**参数**：
- `-f`：输入 XVG 文件（可多个）
- `-c`：列索引选择（每组文件对应一个索引组）
- `-l`：输出数据的标签
- `-x`：X 轴标签
- `-o`：输出文件名

**示例**：
```bash
# 组合两个文件的数据
dit xvg_combine -f RMSD.xvg Gyrate.xvg -c 0,1 1 -l RMSD Gyrate -x "Time(ps)"

# 组合并输出到新文件
dit xvg_combine -f f1.xvg f2.xvg -c 1,2 2,3 -o res.xvg
```

---

### xvg_ave_bar

**功能**：计算多组数据的平均值和误差，并绘制柱状图。

**说明**：
- 适用于比较多个平行模拟的结果
- 自动计算每个体系的平均值和误差
- 支持将数据导出到 CSV

**参数**：
- `-f`：输入 XVG 文件（逗号分隔表示平行实验）
- `-c`：选择要统计的列
- `-l`：每个体系的标签
- `-al`：X 轴标签
- `-csv`：导出数据到 CSV
- `-y`：Y 轴标签
- 其他标准绘图参数

**示例**：
```bash
# 比较三个配体体系的氢键数量
dit xvg_ave_bar -f bar_0_0.xvg,bar_0_1.xvg bar_1_0.xvg,bar_1_1.xvg -c 1,2 -l MD_0 MD_1 -al Hbond Pair -csv hhh.csv -y Number
```

---

### xvg_rama

**功能**：将 Ramachandran 数据（phi 和 psi 二面角）绘制成拉式图。

**说明**：
- 读取 `gmx rama` 命令生成的 XVG 文件
- 借鉴自 PyRAMA 项目

**参数**：
- `-f`：输入 XVG 文件（包含 phi 和 psi 数据）
- 其他标准绘图参数

**示例**：
```bash
# 绘制拉式图
dit xvg_rama -f rama.xvg
```

---

### xvg_energy_compute

**功能**：计算蛋白质和配体之间的能量。

**参数**：
- `-f`：输入 XVG 文件
- 其他标准参数

---

## XPM 命令

### xpm_show

**功能**：可视化 XPM 矩阵文件。

**说明**：
- 支持 4 种绘图引擎：matplotlib, plotly, gnuplot, plotext
- 支持 4 种绘图模式：imshow, pcolormesh, 3d, contour
- matplotlib 支持 imshow, pcolormesh, 3d, contour
- plotly 和 gnuplot 支持 pcolormesh, 3d, contour
- plotext 仅支持绘制小尺寸灰度图
- Discrete 类型 XPM：matplotlib 的 imshow、plotly 和 gnuplot 的 pcolormesh 使用原始颜色
- Continuous 类型 XPM：所有引擎使用 colormap 着色
- 支持插值和图像切割

**参数**：
- `-f`：输入 XPM 文件
- `-m`：绘图模式（imshow, pcolormesh, 3d, contour）
- `-eg`：绘图引擎（matplotlib, plotly, gnuplot, plotext）
- `-cmap`：颜色映射
- `-ip`：插值方法
- `-ipf`：插值倍数
- `-xmin`, `-xmax`, `-ymin`, `-ymax`：图像切割（像素索引）
- `-zmin`, `-zmax`：颜色值范围
- 其他标准绘图参数

**示例**：
```bash
# 基本显示
dit xpm_show -f DSSP.xpm
dit xpm_show -f fel.xpm

# 不显示图形，保存到文件
dit xpm_show -f DSSP.xpm -ns -o dssp.png

# 使用 pcolormesh 模式并插值
dit xpm_show -f FEL.xpm -m pcolormesh -ip linear -ipf 5 -cmap Greys_r

# 3D 模式
dit xpm_show -f FEL.xpm -m 3d -x PC1 -y PC2 -z Energy -t FEL --alpha 0.5
dit xpm_show -f FEL.xpm -m 3d --x_precision 1 --y_precision 2 --z_precision 0 -cmap summer --colorbar_location bottom

# 等高线模式
dit xpm_show -f FEL.xpm -m contour -cmap jet
dit xpm_show -f FEL.xpm -m contour -cmap jet -zmin 0 -zmax 20

# 图像切割（显示部分区域）
dit xpm_show -f DSSP.xpm -xmin 1000 -xmax 2001 -ymin 50 -ymax 101

# 使用 plotly 引擎
dit xpm_show -f FEL.xpm -eg plotly -m 3d -cmap spectral
dit xpm_show -f FEL.xpm -eg plotly -m contour

# 使用 gnuplot 引擎
dit xpm_show -f DSSP.xpm -eg gnuplot --legend_location outside
dit xpm_show -f FEL.xpm -eg gnuplot -m 3d -ip cubic
dit xpm_show -f FEL.xpm -eg gnuplot -m contour -ns -o contour.png

# 设置刻度数量
dit xpm_show -f dccm.xpm --x_numticks 5 --y_numticks 5 --z_numticks 5 -zmin -1
```

---

### xpm2csv

**功能**：将 XPM 文件转换为 CSV 格式。

**说明**：
- 输出三列数据：x, y, v（横坐标、纵坐标、值）

**参数**：
- `-f`：输入 XPM 文件
- `-o`：输出 CSV 文件名

**示例**：
```bash
dit xpm2csv -f test.xpm -o test.csv
```

---

### xpm2dat

**功能**：将 XPM 文件转换为 M×N 的 DAT 格式。

**说明**：
- 第一行注释：X 轴标题、Y 轴标题、矩阵标题
- 第二行：X 轴数据
- 第三行：Y 轴数据（从下到上）
- 后续：M×N 数据矩阵

**参数**：
- `-f`：输入 XPM 文件
- `-o`：输出 DAT 文件名

**示例**：
```bash
dit xpm2dat -f test.xpm -o test.dat
```

---

### xpm_diff

**功能**：计算两个相同尺寸 XPM 文件的差值。

**说明**：
- 适用于比较不同模拟的残基接触矩阵差异
- 适用于比较不同 DSSP 的差异

**参数**：
- `-f`：输入两个 XPM 文件
- `-o`：输出 XPM 文件名

**示例**：
```bash
dit xpm_diff -f DCCM0.xpm DCCM1.xpm -o DCCM0-1.xpm
```

---

### xpm_merge

**功能**：将两个相同尺寸的 XPM 文件沿对角线一半一半拼接。

**说明**：
- 适用于沿对角线对称的 XPM 矩阵图
- 节省篇幅

**参数**：
- `-f`：输入两个 XPM 文件
- `-o`：输出 XPM 文件名

**示例**：
```bash
dit xpm_merge -f DCCM0.xpm DCCM1.xpm -o DCCM0-1.xpm
```

---

## NDX 命令

### ndx_add

**功能**：向 NDX 索引文件添加新的索引组。

**说明**：
- 使用 `-c` 和 `-al` 参数添加新组
- 索引从1开始（GROMACS 约定）

**参数**：
- `-f`：输入 NDX 文件
- `-o`：输出 NDX 文件名
- `-al`：组名（可多个）
- `-c`：原子索引（从1开始）

**示例**：
```bash
# 添加单个组
dit ndx_add -f index.ndx -o test.ndx -al lig -c 1-10

# 添加多个组
dit ndx_add -al lig mol -c 1-10-3,11-21 21-42
```

---

### ndx_split

**功能**：将一个索引组均匀切分成多个组。

**参数**：
- `-f`：输入 NDX 文件
- `-al`：组名和切分数量
- `-o`：输出 NDX 文件名（可选）

**示例**：
```bash
# 按组索引切分
dit ndx_split -f index.ndx -al 1 2

# 按组名切分
dit ndx_split -f index.ndx -al Protein 2

# 切分并输出到新文件
dit ndx_split -f index.ndx -al Protein 2 -o test.ndx
```

---

### ndx_show

**功能**：显示 NDX 文件中所有索引组的名称。

**参数**：
- `-f`：输入 NDX 文件

**示例**：
```bash
dit ndx_show -f test.ndx
```

---

### ndx_rm_dup

**功能**：删除 NDX 文件中所有重复的索引组。

**说明**：
- 重复是指名称和索引都相同

**参数**：
- `-f`：输入 NDX 文件
- `-o`：输出 NDX 文件名

**示例**：
```bash
dit ndx_rm_dup -f test.ndx -o res.ndx
```

---

## 其他命令

### mdp_gen

**功能**：生成 GROMACS MDP 文件模板。

**说明**：
- 提供常见生物体系模拟的 MDP 模板
- **重要**：生成的模板不一定适合您的体系，请根据实际情况调整参数
- 仅用于快速生成基础模板，避免手动复制

**参数**：
- `-o`：输出 MDP 文件名（可选）

**示例**：
```bash
# 生成默认模板
dit mdp_gen

# 生成指定文件名的模板
dit mdp_gen -o nvt.mdp
```

---

### show_style

**功能**：显示或生成不同绘图引擎的格式控制文件。

**说明**：
- 默认显示 DIT 使用的默认格式控制文件
- 可以通过 `-eg` 指定绘图引擎
- 生成的样式文件放在工作目录会自动被 DIT 加载

**参数**：
- `-eg`：绘图引擎（matplotlib, plotly, gnuplot）
- `-o`：输出文件名

**示例**：
```bash
# 显示默认 matplotlib 样式
dit show_style

# 显示 plotly 样式
dit show_style -eg plotly

# 显示 gnuplot 样式
dit show_style -eg gnuplot

# 生成 plotly 样式文件
dit show_style -eg plotly -o DIT_plotly.json
```

---

### find_center

**功能**：寻找 GRO 文件中原子组的几何中心。

**说明**：
- 可以指定索引文件和索引组
- 不指定索引时计算所有原子的几何中心
- `-m AllAtoms`：在所有原子中寻找指定原子组的几何中心

**参数**：
- `-f`：输入 GRO 文件
- `-m AllAtoms`：在所有原子中寻找中心
- （可选）索引文件和索引组

**示例**：
```bash
# 计算所有原子的几何中心
dit find_center -f test.gro

# 使用索引文件
dit find_center -f test.gro index.ndx

# 在所有原子中寻找指定组的几何中心
dit find_center -f test.gro index.ndx -m AllAtoms
```

---

### dccm_ascii

**功能**：将 GROMACS `covar` 命令的 ASCII 格式协方差矩阵转换为动态互相关矩阵（DCCM）XPM 文件。

**参数**：
- `-f`：输入协方差矩阵 DAT 文件
- `-o`：输出 DCCM XPM 文件

**示例**：
```bash
dit dccm_ascii -f covar.dat -o dccm.xpm
```

---

### dssp

**功能**：处理 GROMACS 2023 的 `dssp` 命令生成的 DAT 文件，转换为 DSSP XPM 和 sc.xvg 文件。

**说明**：
- 将 GROMACS 2023 格式转换为旧版本格式
- 生成 XPM 矩阵图和二级结构含量 XVG 文件

**参数**：
- `-f`：输入 DSSP DAT 文件
- `-o`：输出 DSSP XPM 文件
- `-c`：残基索引（用于时间序列生成）
- `-b`, `-e`, `-dt`：时间序列参数
- `-x`：X 轴标签

**示例**：
```bash
# 基本转换
dit dssp -f dssp.dat -o dssp.xpm

# 生成时间序列
dit dssp -f dssp.dat -c 1-42,1-42,1-42 -b 1000 -e 2001 -dt 10 -x "Time (ps)"
```

---

## 重要提示

### 索引规则

- **列索引**：从 0 开始
- **行索引**：从 0 开始
- **原子索引（NDX）**：从 1 开始（GROMACS 约定）
- **范围区间**：`-` 表示左闭右开区间，如 `1-7` 表示第1到第6列

### 数据范围

- `-b` 指定的起始行**包含**在计算中
- `-e` 指定的结束行**不包含**在计算中
- 如果要计算到末尾，只指定 `-b` 即可

### 绘图引擎选择

- **matplotlib**：默认引擎，功能最全面，支持所有模式和参数
- **plotly**：交互式图形，适合展示和分享，支持 3D 和轮廓图
- **gnuplot**：适合批量处理和高性能绘图，需自行安装 gnuplot
- **plotext**：适合在终端中快速查看简单图形

### 颜色映射

- matplotlib 和 plotly 支持自定义 colormap
- 常用 colormap：coolwarm, jet, viridis, plasma, Blues, Reds, Greys 等
- 可以随便写一个 colormap 名称，错误信息会列出所有可用的 colormap

### LaTeX 支持

- 图例、轴标签、标题支持 LaTeX 语法
- 例如：`-l "$\Delta G_{energy}$"`, `-x "$nm^2$"`

### 获取帮助

- 查看所有命令：`dit`
- 查看全局参数：`dit -h`
- 查看特定命令帮助：`dit <command> -h`