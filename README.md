# Sound source tracking项目代码
 - 编程语言：MATLAB
 - 编程环境：MATLAB r2017a for mac
# 文件目录
## 根目录下文件
|文件名称|说明|
|:---|:---|
|CostFunction.m|用于PSO程序中计算粒子的适应度|
|figure_generator.m|生成论文配图|
|framing.m|对音频信号进行分帧|
|gcc_phat_w.m|使用以phat作为参数的GCC程序计算TDOA|
|langevin.m|郎之万运动模型|
|main_w_npso.m|改进后的PSOPF算法主函数|
|main_w_pso.m|未改进的PSOPF算法主函数|
|main_w.m|基本粒子滤波算法主函数|
|multinomialR.m|多项式随机重采样|
|particle_weight_generator.m|计算粒子权重|
|raw.wav|音频信号|
|rir_example.m|使用Image模型产生RIR|
|rir_generator.cpp|Image模型实现|
|tdoaT_generator.m|计算真实TDOA|
|test_w.m|用于生成空间中各处的权重|
|wavread.m|wav文件读取函数|
## 根目录下文件夹
|文件夹名称|说明|
|:---|:---|
|jpg|存储粒子在各时刻的位置图|
|mat|存储main_w.m、main_w_pso.m、main_w_npso.m三个主函数生成的数据，以及figure_generator.m生成的论文配图|
 - jpg文件夹和mat文件夹中又依据三种不同的运动轨迹为依据创建了semicircle、straightLine、triangle三种文件夹；
 - jpg文件夹中的这三种文件夹中又依据使用算法的不同创建了ps、psopf、npsopf三种文件夹，它们中便是存储`粒子在各时刻的位置图`的位置；
 - mat文件夹中的这三种文件夹中便是存储`main_w.m、main_w_pso.m、main_w_npso.m三个主函数生成的数据`的位置，mat文件夹中的jpg文件夹是存储`figure_generator.m生成的论文配图`的位置。
# 引用其他项目
 - gcc_phat_w.m来自[LucasVale22/ufrj-programas](https://github.com/LucasVale22/ufrj-programas/blob/ce9f7da61604d233e6edbd879647755c2caae90d/IC-PADS/Testes%20Gravacao%202.0/doa_gcc0603.m)，使用中进行了更改；
 - rir_example.m与rir_generator.cpp来自[ehabets/RIR-Generator](https://github.com/ehabets/RIR-Generator)；
 - 另有多个文件是由书籍[《粒子滤波原理及应用――MATLAB仿真》](https://item.jd.com/12178346.html)中的代码以及网上代码更改而来。
# 使用
## 编译运行
 1. 编译`rir_generator.cpp`
```
mex rir_generator.cpp
```
 2. 逐个运行`main_w.m`、`main_w_pso.m`、`main_w_npso.m`；
 3. 在`/jpg/semicircle`或`/jpg/straightLine`或`/jpg/triangle`目录下查看`粒子在各时刻的位置图`；
 4. 运行`figure_generator.m`；
 5. 在`/mat/jpg`目录下查看`figure_generator.m生成的论文配图`。
## 参数设置
 - 可以通过在`main_w.m`、`main_w_pso.m`、`main_w_npso.m`三个文件的第29行以及`figure_generator.m`的第5行设置`track`变量的值来改变说话人运行的轨迹。其中，`track = 1`表示圆形轨迹；`track = 2`表示三角轨迹；`track = <其他>`表示直线轨迹。
# 开源协议
MIT License
