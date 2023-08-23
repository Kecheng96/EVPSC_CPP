The cpp version of Elastic Visco-plastic Self-Consistent model
- # Elastic visco-plastic self-consistent model
- ## 1. 晶体变形运动学
	- 描述单晶体大变形，即描述单晶体从初始构型（参考构型）变形到当前构型（变形后构型）时的几何变化和应力状态变化。我们定义：
	- $\mathbf X$ ：物质点在晶粒初始构型中的坐标
	  $\mathbf {x(X)}$ ：物质点在晶粒当前构型的坐标  
	  $\mathbf {u=x-X}$ ：物质点的位移  
	- 晶粒的变形可以通过变形梯度张量 $\mathbf F$ 及速度梯度张量 $\mathbf l$ 来描述：
	- $\mathbf {F} = \frac{\partial \mathbf x}{\partial \mathbf X}$
	- $\mathbf {l} = \frac{\partial \mathbf {\dot u}}{\partial \mathbf x}$
