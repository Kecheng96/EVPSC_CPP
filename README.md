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
$$
\begin{align}
\boldsymbol {l}&=\boldsymbol {\dot F} \cdot \boldsymbol{F}^{-1}\\
&=(\boldsymbol{\dot F}^e \cdot \boldsymbol {F}^p + \boldsymbol{F}^e \cdot \boldsymbol {\dot F}^p)\cdot(\boldsymbol {F}^e \cdot \boldsymbol {F}^p)^{-1}\\
&= \boldsymbol{\dot F}^e \cdot (\boldsymbol{F}^e)^{-1} +
\boldsymbol{ F}^e \cdot \boldsymbol{F}^p\cdot(\boldsymbol{F}^p)^{-1}\cdot(\boldsymbol{F}^e)^{-1}
\end{align}
$$

令
