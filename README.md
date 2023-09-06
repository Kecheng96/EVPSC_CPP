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

其中， $\dot\gamma^\alpha$ 为变形系 $\alpha$ 引起的剪切应变率。

在发生晶格畸变后，晶格矢量将发生拉伸和转动，他们会保持正交关系，但一般不再是单位矢量。晶格畸变后，变形方向和法向分别为 $\boldsymbol {s}^\alpha$ ， $\boldsymbol {n}^\alpha$ :
$$\boldsymbol {s}^\alpha=\boldsymbol{F}^e\cdot\boldsymbol {s}^\alpha_0\tag{1-11a}$$
$$\boldsymbol {n}^\alpha=(\boldsymbol{F}^e)^{-T}\cdot\boldsymbol {n}^\alpha_0\tag{1-11a}$$

那么塑性速度梯度张量 $\boldsymbol {l}^p$:
$$\boldsymbol {l}^p=\sum_\alpha \dot\gamma^\alpha \boldsymbol {s}^\alpha\cdot(\boldsymbol {n}^\alpha)^T\tag{1-12}$$

根据公式(1-12)，可以得到应变率张量 $\boldsymbol d$ 和旋率张量 $\boldsymbol w$ 的弹性和塑性部分
$$\boldsymbol{d}^e= \frac{1}{2}[\boldsymbol{\dot F}^e\cdot(\boldsymbol{F}^e)^{-1}+(\boldsymbol{F}^e)^{-T}\cdot(\boldsymbol{\dot F}^e)^T]\tag{1-13}$$
$$\boldsymbol{w}^e= \frac{1}{2}[\boldsymbol{\dot F}^e\cdot(\boldsymbol{F}^e)^{-1}-(\boldsymbol{F}^e)^{-T}\cdot(\boldsymbol{\dot F}^e)^T]\tag{1-14}$$
$$\boldsymbol {d}^p=\sum_\alpha \dot\gamma^\alpha \boldsymbol {P}^\alpha\tag{1-15}$$
$$\boldsymbol {w}^p=\sum_\alpha \dot\gamma^\alpha \boldsymbol {R}^\alpha\tag{1-16}$$

其中 $\boldsymbol {P}^\alpha=\frac{1}{2}[\boldsymbol {s}^\alpha\cdot(\boldsymbol {n}^\alpha)^T+\boldsymbol {n}^\alpha\cdot(\boldsymbol {s}^\alpha)^T]$（称为Schmid tensor施密特张量）,
$\boldsymbol {R}^\alpha=\frac{1}{2}[\boldsymbol {s}^\alpha\cdot(\boldsymbol {n}^\alpha)^T-\boldsymbol {n}^\alpha\cdot(\boldsymbol {s}^\alpha)^T]$

## 2. 单晶本构关系
设晶体的弹性性质不受滑移/孪晶变形的影响，则单晶的本构方程为:
$$\boldsymbol\sigma^{\nabla*}+\boldsymbol\sigma\ tr(\boldsymbol{d}^e)=\boldsymbol L:\boldsymbol{d}^e\tag{2-1}$$

其中， $\boldsymbol L$ 为四阶弹性模量张量， $\boldsymbol\sigma^{\nabla*}$ 为基于中间构型的Cauchy应力张量的客观率（Jaumann率）:
$$\boldsymbol\sigma^{\nabla*}=\boldsymbol {\dot\sigma}-\boldsymbol {w}^e\cdot\boldsymbol\sigma+\boldsymbol\sigma\cdot\boldsymbol {w}^e\tag{2-2}$$

而基于初始构型的Cauchy应力张量的Jaumann率为：

$$\begin{align}\boldsymbol{\sigma}^\nabla &= \boldsymbol{\dot\sigma}-\boldsymbol w \cdot \boldsymbol \sigma + \boldsymbol \sigma \cdot \boldsymbol w \\
&= (\boldsymbol{\sigma}^{\nabla*}+\boldsymbol w^e\cdot\boldsymbol\sigma-\boldsymbol \sigma\cdot\boldsymbol w^e) - \boldsymbol w \cdot \boldsymbol \sigma + \boldsymbol \sigma\cdot \boldsymbol w\\
&=\boldsymbol\sigma^{\nabla*}-(\boldsymbol {w}-\boldsymbol {w}^e)\cdot\boldsymbol\sigma+\boldsymbol\sigma\cdot(\boldsymbol {w}-\boldsymbol {w}^e)\\
&=\boldsymbol\sigma^{\nabla*}-\boldsymbol {w}^p\cdot\boldsymbol\sigma+\boldsymbol\sigma\cdot\boldsymbol {w}^p\end{align}\tag{2-3}$$

将式(2-3)代入单晶本构方程(2-1)，得到在参考构型下的单晶本构方程：
$$\boldsymbol{\sigma}^\nabla+\boldsymbol {w}^p\cdot\boldsymbol\sigma-\boldsymbol\sigma\cdot\boldsymbol {w}^p + \boldsymbol\sigma\ tr(\boldsymbol{d}-\boldsymbol{d}^p)=\boldsymbol L : (\boldsymbol{d}-\boldsymbol{d}^p)\tag{2-4}$$

整理得到:
$$\boldsymbol{\sigma}^\nabla= \boldsymbol {L'}:(\boldsymbol{d}-\boldsymbol{d}^p)+\boldsymbol\sigma^0\tag{2-5}$$

其中, $L'_{ijkl}=L_{ijkl}-\sigma_{ij}\delta_{kl}$ , $\sigma^0_{ij}=w^p_{ik}\sigma_{kj}-\sigma_{ik}w^p_{kj}$ , 进而得到应变率张量与客观应力率的关系:
$$\boldsymbol d = \boldsymbol M^e:\boldsymbol\sigma^\nabla+\boldsymbol d^p + \boldsymbol w^0\tag{2-6}$$

其中 $\boldsymbol M^e=(\boldsymbol L')^{-1}$ 为弹性柔度张量， $\boldsymbol w^0=\boldsymbol M^e:\boldsymbol\sigma^0$ 

根据式(1-14a)，单晶的本构关系中还需要明确 $\alpha$ 滑移/孪生系的剪切应变率 $\dot\gamma^\alpha$ ，而对率相关材料， $\dot\gamma^\alpha$ 取决于变形系的分解剪切应力 $\tau^\alpha$ 、临界剪切应力 $\tau_{cr}^\alpha$ 以及率相关系数 $m$ 等:
$$\dot\gamma^\alpha=\dot\gamma^\alpha(\tau^\alpha,\tau^\alpha_{cr},m,...)\tag{2-7}$$
其中，变形系的分解剪切应力 $\tau^\alpha=\boldsymbol P^\alpha: \boldsymbol\sigma'$ ,  $\boldsymbol\sigma'$ 为应力偏张量。可以看出，这样定义的分解剪切应力 $\tau^\alpha$ 与剪切应变率 $\dot\gamma^\alpha$ 是功共轭的；而临界剪切应力 $\tau^\alpha_{cr}$ 反映了变形系的硬化/软化行为，其变化率 $\dot\tau^\alpha_{cr}$ 与当前临界剪切应力 $\tau^\alpha_{cr}$ 、其他变形系 $\beta$ 的累计剪切应变 $\gamma^\beta$ 和剪切应变率 $\dot\gamma^\beta$ 、以及孪晶系 $\kappa$ 的孪晶体积分数 $f^\kappa$ 相关:
$$\dot\tau^\alpha_{cr}=\dot\tau^\alpha_{cr}(\tau^\alpha_{cr},\gamma^\beta,\dot\gamma^\beta,f^\kappa,...)\tag{2-8}$$

其中，孪晶系 $\kappa$ 的孪晶体积分数 $f^\kappa$ 的变化率为:
$$\dot{f}^\kappa=\frac{\gamma^\kappa}{\gamma^{tw}}\tag{2-9}$$
$\gamma^{tw}$ 为孪晶系的特征剪切应变值，为常数。

根据 Tomé等人的研究，对滑移系:
$$\dot\gamma^\alpha=\dot\gamma_0|\frac{\tau^\alpha}{\tau^\alpha_{cr}}|^\frac{1}{m}sgn(\tau^\alpha)\tag{2-10}$$

考虑到孪晶变形的极性，对孪晶系:

$$\dot\gamma^\alpha=\begin{cases}\dot\gamma_0|\frac{\tau^\alpha}{\tau^\alpha_{cr}}|^\frac{1}{m} & , \tau^\alpha \gt 0 \\
1 & , \tau^\alpha \le 0\end{cases}\tag{2-11}$$

$\dot\gamma_0$ 为参考剪切应变率， $sgn$ 为符号函数。

由此建立了塑性应变率 $\boldsymbol d^p$ 与应力张量 $\boldsymbol\sigma$ 的联系，注意到此时单晶体本构关系(2-6)成为一个非线性方程，可以通过不同的方法将该方程进行准线性化，
本构关系可以进一步表示为:
$$\boldsymbol d = \boldsymbol M^e:\boldsymbol\sigma^\nabla+\boldsymbol M^{vp}:\boldsymbol\sigma'+\boldsymbol d^0 \tag{2-12}$$
其中, $\boldsymbol M^{vp}$ 为粘塑性模量， $\boldsymbol d^0$ 为使该准线性方程成立的反推项，并且有 $\boldsymbol d^e=\boldsymbol M^e:\boldsymbol\sigma^\nabla$, $\boldsymbol d^p=\boldsymbol M^{vp}:\boldsymbol\sigma^{'}+\boldsymbol d^0$ 

## 3. 多晶体自洽模型
