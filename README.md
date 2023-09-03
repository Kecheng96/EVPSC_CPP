The cpp version of Elastic Visco-plastic Self-Consistent model
# Elastic visco-plastic self-consistent model
## 1. 晶体变形运动学
描述单晶体大变形，即描述单晶体从初始构型（参考构型）变形到当前构型（变形后构型）时的几何变化和应力状态变化。

我们定义：
$\boldsymbol X$ ：物质点坐标，即在晶粒参考构型中的坐标
$\boldsymbol{x(X)}$ ：物质点在晶粒当前构型的坐标
$\boldsymbol {u=x-X}$ ：物质点的位移

根据有限变形理论，晶粒的变形可以通过变形梯度张量 $\boldsymbol  F$ 及速度梯度张量 $\boldsymbol l$ ，定义为：
$$\boldsymbol {l} = \frac{\partial \boldsymbol {\dot u}}{\partial \boldsymbol x} = \frac{\partial \boldsymbol v}{\partial \boldsymbol x} \tag{1-1}$$

$$\boldsymbol{ F= \frac{\partial  x}{\partial X}= \frac{\partial  u}{\partial  X} + I }\tag{1-2}$$	

根据它们的偏微分关系，有：
$$\boldsymbol {\dot F} = \boldsymbol {l \cdot F}  \tag{1-3}$$

如图1所示，变形梯度张量 $\boldsymbol {F}$ 可以分解为：
$$\boldsymbol {F} = \boldsymbol {F}^e \cdot \boldsymbol {F}^p \tag{1-4}$$

晶体在外力的作用下， 会发生晶格畸变，同时由于晶粒边界的约束和变形协调的要求，发生刚体转动， $\boldsymbol {F}^e$ 即表示由晶格畸变和刚体转动所产生的变形梯度， $\boldsymbol {F}^p$ 则表示晶体由于滑移/孪晶系统产生的均匀剪切产生的变形梯度。

![Fig1](/PNGs/F=FeFp.png)
图 1. 晶体弹塑性变形几何学

根据公式(1-3)和公式(1-4)，速度梯度张量可以写成：

$$\begin{align}\boldsymbol l &= \boldsymbol {\dot F} \cdot \boldsymbol {F}^{-1}\\
&= (\boldsymbol {\dot F}^e \cdot \boldsymbol { F}^p+\boldsymbol { F}^e \cdot \boldsymbol {\dot F}^p)\cdot(\boldsymbol { F}^e \cdot \boldsymbol {F}^p)^{-1}\\
&=\boldsymbol {\dot F}^e \cdot (\boldsymbol { F}^e)^{-1} + \boldsymbol { F}^e \cdot \boldsymbol {\dot F}^p \cdot (\boldsymbol {F}^p)^{-1} \cdot  (\boldsymbol {F}^e)^{-1}\\
\end{align}\tag{1-5}$$

令 $\boldsymbol {l}^e=\boldsymbol {\dot F}^e \cdot (\boldsymbol { F}^e)^{-1}$ , $\boldsymbol {l}^p= \boldsymbol { F}^e \cdot \boldsymbol {\dot F}^p \cdot (\boldsymbol {F}^p)^{-1} \cdot  (\boldsymbol {F}^e)^{-1}$ 则速度梯度分解成弹性和塑性部分：
$$\boldsymbol {l} = \boldsymbol {l}^e+\boldsymbol {l}^p \tag{1-6}$$

速度梯度张量可以分解成对称张量（应变率张量）$\boldsymbol d$ 和反对称张量（旋率张量）$\boldsymbol w$:
$$\boldsymbol {l = d + w}\tag{1-7}$$

$\boldsymbol {l}$ 的弹性部分和塑性部分也可以进行类似的分解：
$$\boldsymbol {l}^e = \boldsymbol {d}^e + \boldsymbol {w}^e\tag{1-8}$$
$$\boldsymbol {l}^p = \boldsymbol {d}^p + \boldsymbol {w}^p\tag{1-9}$$

单晶体材料的塑性变形由滑移或孪生引起，在图 1所示的中间构型中，晶格矢量不发生变化，记第 $\alpha$ 个滑移/孪晶系统的变形方向为 $\boldsymbol {s}^\alpha_0$、变形法向为 $\boldsymbol {n}^\alpha_0$，则发生的塑性变形：
$$\boldsymbol {\dot F}^p \cdot (\boldsymbol {F}^p)^{-1} = \sum_\alpha \dot\gamma^\alpha\boldsymbol {s}^\alpha_0\cdot(\boldsymbol {n}^\alpha_0)^T\tag{1-10}$$

其中，$\dot\gamma^\alpha$ 为变形系 $\alpha$ 引起的剪切应变率。

在发生晶格畸变后，晶格矢量将发生拉伸和转动，他们会保持正交关系，但一般不再是单位矢量。晶格畸变后，变形方向和法向分别为 $\boldsymbol {s}^\alpha$，$\boldsymbol {n}^\alpha$:
$$\boldsymbol {s}^\alpha=\boldsymbol{F}^e\cdot\boldsymbol {s}^\alpha_0\tag{1-11a}$$
$$\boldsymbol {n}^\alpha=(\boldsymbol{F}^e)^{-T}\cdot\boldsymbol {n}^\alpha_0\tag{1-11a}$$

那么塑性速度梯度张量 $\boldsymbol {l}^p$:
$$\boldsymbol {l}^p=\sum_\alpha \dot\gamma^\alpha \boldsymbol {s}^\alpha\cdot(\boldsymbol {n}^\alpha)^T\tag{1-12}$$

根据公式(1-12)，可以得到应变率张量 $\boldsymbol d$ 和旋率张量 $\boldsymbol w$ 的弹性和塑性部分
$$\boldsymbol{d}^e= \frac{1}{2}[\boldsymbol{\dot F}^e\cdot(\boldsymbol{F}^e)^{-1}+(\boldsymbol{F}^e)^{-T}\cdot(\boldsymbol{\dot F}^e)^T]\tag{1-13}$$
$$\boldsymbol{w}^e= \frac{1}{2}[\boldsymbol{\dot F}^e\cdot(\boldsymbol{F}^e)^{-1}-(\boldsymbol{F}^e)^{-T}\cdot(\boldsymbol{\dot F}^e)^T]\tag{1-14}$$
$$\boldsymbol {d}^p=\sum_\alpha \dot\gamma^\alpha \boldsymbol {P}^\alpha\tag{1-15}$$
$$\boldsymbol {w}^p=\sum_\alpha \dot\gamma^\alpha \boldsymbol {R}^\alpha\tag{1-16}$$

其中 $\boldsymbol {P}^\alpha=\frac{1}{2}[\boldsymbol {s}^\alpha\cdot(\boldsymbol {n}^\alpha)^T+\boldsymbol {n}^\alpha\cdot(\boldsymbol {s}^\alpha)^T]$（称为Schmid tensor施密特张量）,$\boldsymbol {P}^\alpha=\frac{1}{2}[\boldsymbol {s}^\alpha\cdot(\boldsymbol {n}^\alpha)^T+\boldsymbol {n}^\alpha\cdot(\boldsymbol {s}^\alpha)^T]$

## 2. 单晶本构关系
设晶体的弹性性质不受滑移/孪晶变形的影响，则单晶的本构方程为:
$$\boldsymbol\sigma^{\nabla*}+\boldsymbol\sigma\ tr(\boldsymbol{d}^e)=\boldsymbol L:\boldsymbol{d}^e\tag{2-1}$$
