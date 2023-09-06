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

$$\boldsymbol{F}= \frac{\partial  \boldsymbol x}{\partial \boldsymbol X}= \frac{\partial \boldsymbol u}{\partial \boldsymbol X} + \boldsymbol I \tag{1-2}$$	

根据它们的偏微分关系，有：
$$\boldsymbol {\dot F} = \boldsymbol {l \cdot F}  \tag{1-3}$$

如图1所示，变形梯度张量 $\boldsymbol {F}$ 可以分解为：
$$\boldsymbol {F} = \boldsymbol {F}^e \cdot \boldsymbol {F}^p \tag{1-4}$$

晶体在外力的作用下， 会发生晶格畸变，同时由于晶粒边界的约束和变形协调的要求，发生刚体转动， $\boldsymbol {F}^e$ 即表示由晶格畸变和刚体转动所产生的变形梯度， $\boldsymbol {F}^p$ 则表示晶体由于滑移/孪晶系统产生的均匀剪切产生的变形梯度。

<div align=center><img src="/PNGs/F=FeFp.png" width = "50%" /></div>
<div align=center>图 1. 晶体弹塑性变形几何学</div>


根据公式(1-3)和公式(1-4)，速度梯度张量可以写成：

$$\begin{align}\boldsymbol l &= \boldsymbol {\dot F} \cdot \boldsymbol {F}^{-1}\\
&= (\boldsymbol {\dot F}^e \cdot \boldsymbol { F}^p+\boldsymbol { F}^e \cdot \boldsymbol {\dot F}^p)\cdot(\boldsymbol { F}^e \cdot \boldsymbol {F}^p)^{-1}\\
&=\boldsymbol {\dot F}^e \cdot (\boldsymbol { F}^e)^{-1} + \boldsymbol { F}^e \cdot \boldsymbol {\dot F}^p \cdot (\boldsymbol {F}^p)^{-1} \cdot  (\boldsymbol {F}^e)^{-1}\\
\end{align}\tag{1-5}$$

令 $\boldsymbol {l}^e=\boldsymbol {\dot F}^e \cdot (\boldsymbol { F}^e)^{-1}$ , $\boldsymbol {l}^p= \boldsymbol { F}^e \cdot \boldsymbol {\dot F}^p \cdot (\boldsymbol {F}^p)^{-1} \cdot  (\boldsymbol {F}^e)^{-1}$ 则速度梯度分解成弹性和塑性部分：
$$\boldsymbol {l} = \boldsymbol {l}^e+\boldsymbol {l}^p \tag{1-6}$$

速度梯度张量可以分解成对称张量（应变率张量） $\boldsymbol d$ 和反对称张量（旋率张量） $\boldsymbol w$ :
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

其中, $L_{ijkl}'=L_{ijkl}-\sigma_{ij}\delta_{kl}$ , $\sigma^0_{ij}=w^p_{ik}\sigma_{kj}-\sigma_{ik}w^p_{kj}$ , 进而得到应变率张量与客观应力率的关系:
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

多晶体材料可以看作是许多单晶的集合，基于单晶体本构模型和以各个晶体的取向，并结合应力应变协调方程可得到多晶体材料在宏观载荷下的力学响应，结合有限元理论即可实现（Crystal-plasticity finite element method, CPFEM方法）。然而，由于实际多晶体材料中包含的晶粒数量较多，通过CPFEM来计算宏观力学响应往往意味着较大的计算开销。多晶体自洽模型则是一种能在较小的计算开销下，得到多晶体材料的宏微观力学响应的均质化假设模型。


对多晶体材料，宏观应变率张量 $\boldsymbol D$ 、旋率张量 $\boldsymbol W$ 及 Cauchy 应力张量 $\boldsymbol\Sigma$ 可看作其包含的所有单晶体的对应张量的体积平均:
$$\boldsymbol D = \langle\boldsymbol d\rangle=\frac{1}{V}\int\boldsymbol d\ dV\tag{3-1a}$$
$$\boldsymbol W = \langle\boldsymbol w\rangle=\frac{1}{V}\int\boldsymbol w\ dV\tag{3-1b}$$
$$\boldsymbol\Sigma = \langle\boldsymbol\sigma\rangle=\frac{1}{V}\int\boldsymbol\sigma\ dV\tag{3-1c}$$

$V$ 为多晶体的体积，算符〈⋯〉表示求体积平均。均匀化处理之后可以得到多晶体的本构方程：
$$\boldsymbol D = \overline{\boldsymbol M}^e:\boldsymbol\Sigma^\nabla+\overline{\boldsymbol M}^{vp}:\boldsymbol\Sigma'+\boldsymbol D^0 \tag{3-2}$$

其中， $\overline{\boldsymbol M}^e$ , $\overline{\boldsymbol M}^{vp}$ 和 $\boldsymbol D^0$ 分别为宏观弹性模量张量，宏观粘塑性模量张量和反推项。在实际计算中，仅仅只有宏观的边界条件是已知的，而宏观的模量和各个晶粒的模量以及塑性应变都是未知的，需要利用多晶体聚合体和晶粒之间的联系来迭代求解。

### 3.1 椭球体夹杂问题
上文提到的求解宏微观模量的方法，是通过将多晶体看作无限大均匀介质，从而使得需要求解的晶粒与晶粒之间的应力状态差异问题，转换为了根据本征应变求解本征应力的问题。现在假设均匀介质中存在一个局部区域 $\Omega$ ,若假设 $\Omega$ 不受到周围介质的约束，在由于某种物理或化学的原因，则会产生不产生应力场的均匀的局部应变 $\boldsymbol\varepsilon^+$，该局部应变就被称为本征应变。本征应变是一个广义的概念，它可以是热应变、相变应变和残余应力等，实际上区域 $\Omega$ 在周围介质的约束下，无法自由发生变形，从而在局部区域 $\Omega$ 内外产生应力场，产生本征应变的区域 $\Omega$ 称为夹杂。

本征应变问题可以分解为三个问题的叠加
	
（i）将局部区域 $\Omega$ 剥离，让 $\Omega$ 自由产生本征应变 $\varepsilon_{ij}^+$ ，此时区域 $\Omega$ 内没有产生应力，而剩余区域 $\Omega^-$ 也不产生应变；

（ii）当本征应变在区域 $\Omega$ 内均匀分布，则可以通过在夹杂边界 $S$ 附加虚拟面力 $p_i^+$ ，从而使区域 $\Omega$ 产生弹性应变 $-\varepsilon_{ij}^+$ ，恢复取出时的形状，这样产生的弹性应力场为：


$$\sigma_{ij}^+=-C_{ijkl}\varepsilon_{ij}^+\tag{3-3}$$

式中 $C_{ijkl}$ 为材料的弹性刚度， $σ_{ij}^+$ 即为对应的本征应力，虚拟面力 $p_i^+$ 为:

$$p_i^+=-\sigma^{+}_{ij}n_j\tag{3-4}$$

式中 $n_j$ 为边界的外法向。至此，局部区域 $\Omega$ 已经恢复成原来的形状，只是在边界上存在虚拟面力;

（iii）将局部区域 $\Omega$ 放回均匀介质中，在夹杂边界上施加 $-p_i^+$ ，同时让整个均匀介质一起变形，释放虚拟载荷。求解此时无限大介质内的位移解即为夹杂问题的位移解。

Eshebly (1957) 证明，当介质为线弹性，夹杂体形状为椭球体，而且本征应变 $\boldsymbol\varepsilon^+$ 为常应变（应变大小在夹杂体内不随位置改变），最终求解得到的夹杂内的实际应变 $\boldsymbol\varepsilon$ 也是常应变，二者之间满足:

$$\varepsilon_{ij}^+=S_{ijkl}\varepsilon_{ij}^+\tag{3-5}$$

$S_{ijkl}$ 称为Eshebly张量，它仅与介质的弹性性质和椭球体的形状与取向有关。 $S_{ijkl}$ 关于 $i$ 和 $j$ ， $k$ 和 $l$ 对称，但一般关于 $(i,j)$ 与 $(k,l)$ 不对称，故一般不具有 Voigt 对称性。

### 3.2 粘塑性介质中粘塑性夹杂问题
将多晶体视为无限大粘塑性介质，而某一晶粒则为夹杂体。根据单晶体塑性应变率 $\boldsymbol d^p=\boldsymbol M^{vp}:\boldsymbol\sigma'+\boldsymbol d^0$ 和多晶体塑性应变率表达式 $\boldsymbol D^p=\overline{\boldsymbol M}^{vp}:\boldsymbol\Sigma'+ \boldsymbol D^0$,将单晶体的塑性应变率通过宏观粘塑性张量整理成:
$$\boldsymbol d^p=\overline{\boldsymbol M}^{vp}:\boldsymbol\Sigma'+\boldsymbol d^0 + \boldsymbol d^+\tag{3-6}$$

这样， $\boldsymbol d^+=(\boldsymbol M^{vp}-\overline{\boldsymbol M}^{vp}):\boldsymbol\sigma'+(\boldsymbol d^0-\boldsymbol D^0)$ 则是此时的本征应变率，考虑到粘塑性刚度张量 $\overline{\boldsymbol L}^{vp}=(\overline{\boldsymbol M}^{vp})^{-1}$ ，并记 $\boldsymbol{\widetilde\sigma}'=\boldsymbol\sigma'-\boldsymbol\Sigma'$ ， $\boldsymbol{\widetilde d}^p=\boldsymbol d^p-\boldsymbol D^p$ 式（3-6）可以改写成：
$$\boldsymbol{\widetilde\sigma}'=\overline{\boldsymbol L}^{vp}:(\widetilde{\boldsymbol d}^p-\boldsymbol d^+)\tag{3-7a}$$

记材料点的坐标为 $\boldsymbol x$ ,并将张量形式展开：

$$\widetilde\sigma_{ij}'(\boldsymbol x)=\overline{L_{ijkl}^{vp}}:(\widetilde d_{kl}^p(\boldsymbol x)-d_{kl}^+(\boldsymbol x)) \tag{3-7b}$$

平衡方程为:

$$\sigma_{ij,j}(\boldsymbol x)=(\widetilde\sigma_{ij}(\boldsymbol x)+\Sigma_{ij}(\boldsymbol x)),j=\widetilde{\sigma}_{ij,j}(\boldsymbol x)=0\tag{3-8}$$

