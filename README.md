# CO2-Carbon-Dioxide-Property-Calculation
High-resolution CO2 (Carbon Dioxide) Properties Calculation

_CO2PC程序说明：

1.	本程序名为_CO2PC.exe（CO2 Properties Calculation)，版本为1.0beta。技术参照：
	CO2热物性计算：R. Span, W. Wagner, A New Equation of State for Carbon Dioxide Covering the Fluid Region from the Triple –Point Temperature to 1100K at Pressures up to 800MPa, J. Phys. Chem. Ref. Data, Vol 25, No.6, 1996
	CO2热导率计算：Marcia L Huber, Marc J Assel, R. A. Perkins, Reference Correlation of the Thermal Conductivity of Carbon Dioxide from the Triple Point to 1100K and up to 200MPa. Journal of Physical and Chemical References Date. Feb. 2016
	CO2粘度计算：William A. Wakeham, The Viscosity of Carbon Dioxide, Journal of Physical and Chemical Reference Data, Jan, 1998.

2.	程序功能是根据给定压力（p）、温度（t）、密度（r）、比焓（h）、比熵（s）5个参数的两个，计算出一系列热物性参数和输运参数。

3.	程序输入文件：_input.dat。格式为每行三个数据，第一个为整数，表示计算类型（给定参数）：1-给定pt，2-给定ph，3-给定ps，4-给定hs，5-给定tr，6-给定pr，7-给定th，8-给定ts。第二、三个参数为浮点数（可以使用C语言认可的科学计数法，如103=1.03e2。请见范例文件）。请务必注意参数顺序，严格按照先后顺序。各参数之间用逗号隔开。所有行之后不得再有任何空格、回车等。
给定的参数采用国际单位制：
压力：Pa
温度：K
密度：kg/m3
比焓：J/kg
比熵：J/(kgK)

4.	程序输出文件：_properties.dat，格式见范例文件。
Vsound：声速
V_FRAC：两相时的气态质量比
VISC：动力粘度
THCOND：导热系数（热导率）

