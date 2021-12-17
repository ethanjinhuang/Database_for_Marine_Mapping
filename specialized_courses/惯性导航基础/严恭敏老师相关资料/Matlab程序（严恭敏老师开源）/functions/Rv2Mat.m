function m = Rv2Mat(rv)  
%旋转矢量转化为变换矩阵
    m = [1, 0, 0; 0, 1, 0; 0, 0, 1] - Asym(rv);
