一、使用方法：确保Matlab包含当前路径和 .\functions路径，用trajectory生成轨迹数据存盘；ins读数据进行捷联惯导解算；federal_rk直接解微分方程仿真车载组合导航系统；federal 进行车载组合导航综合仿真；plotfigure 为以上个程序运行结果之后的作图程序。

二、程序包括：
	trajectory 模拟轨迹仿真
	ins 捷联惯导解算
	federal_rk 车载组合导航仿真（直接解微分方程）
	federal 车载组合导航综合仿真
	plotfigure 以上个程序运行结果之后的作图程序
	functions下的子程序
		Asym 求反对称阵
		Att2Mat 姿态角转化为姿态矩阵
		Att2Quat 姿态角转化为姿态四元数
		Rv2Mat 旋转矢量转化为变换矩阵
		Rv2Quat 旋转矢量转化为变换四元数
		QuatMul 四元数相乘
		Quat2Mat 姿态四元数转化为姿态矩阵 
		Quat2Att 姿态四元数转化为姿态角
		Mat2Att 姿态矩阵转化为姿态角
		Mat2Quat 姿态矩阵转化为姿态四元数   
		sins 捷联惯导解算
		getf 求系统矩阵Ft
		geth 求观测矩阵Ht
		Kfilter 卡尔曼滤波
