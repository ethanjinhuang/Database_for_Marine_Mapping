global Re e wie g0 ppm ug deg min sec hur dph PKG PKA  %全局变量
Re = 6378160;               %地球半径
e = 1/298.3;                %椭圆度
wie = 7.2921151467e-5;      %自转角速率
g0 = 9.7803267714;          %重力加速度
ppm = 1.0e-6;               %百万分之一
ug = 1.0e-6*g0;             %毫重力加速度
deg = pi/180;               %角度
min = deg/60;               %角分
sec = min/60;               %角秒
hur = 3600;                 %小时
dph = deg/hur;              %度每小时
PKG = 0.932*sec;            %陀螺脉冲当量
PKA = 1/2500*g0;            %加速度计脉冲当量

%马尔可夫过程参数
%陀螺漂移：          eb_=-1/tao(1:3)*eb+web;       
%里程仪刻度系数误差：dKD_=-1/tao(4)+wdKD;
%GPS：         速度：dvnS_=-1/tao(5:7)*dvnS+wdvnS;
%              位置：dposS_=-1/tao(8:10)*dposS+wdposS; 
tao = [3600; 3600; 3600;   10e4;   5; 5; 5; 10; 10; 10];
Rv = [(0.03*dph)^2; (0.03*dph)^2; (0.03*dph)^2;  0.0001^2;  0.1^2; 0.1^2; 0.1^2; (20/Re)^2; (20/Re)^2; 50^2];
%捷联算法、里程仪、GPS更新周期
Tm = 0.02; TD = 0.1; TS = 1;
%马尔可夫过程离散参数
e_tao = exp([-1./tao(1:3)*Tm; -1./tao(4)*TD; -1./tao(5:10)*TS]);
sqw = sqrt( Rv.*(ones(10,1)-exp(-2*[1./tao(1:3)*Tm; 1./tao(4)*TD; 1./tao(5:10)*TS])) );
%系统状态初值设置
fi = [.3 ; .3; 3]*min;     dvn = [0.01; 0.01; 0.01];           dpos = [20/Re; 20/Re; 10];	
dKG = [50; 49; 48]*ppm;    eb = [0.01; 0.01; 0.01]*dph;  
dKA = [100; 90; 80]*ppm;	db = [100; 90; 80]*ug;
dposD = [20/Re; 20/Re; 10];	dKD = 0.001;        
dvnS = [0.01; 0.01; 0.01];     dposS = [20/Re; 20/Re; 50];	
x0 = [fi; dvn; dpos; dKG; eb; dKA; db; dposD; dKD; dvnS; dposS];
%系统方程白噪声均方差
wfi = [.01; 0.01; 0.01]*min; wdvn = [0.001; 0.001; 0.001];    wdpos = [1/Re; 1/Re; 1];
wdKG = [0; 0; 0];            web = sqrt( 2*Rv(1:3).*(1./tao(1:3)) );
wdKA = [0; 0; 0];            wdb = [0; 0; 0];
wdposD = [1/Re; 1/Re; 1];    wdKD = sqrt(2*Rv(4)*1/tao(4));
wdvnS = sqrt( 2*Rv(5:7).*(1./tao(5:7)) );    wdposS = sqrt( 2*Rv(8:10).*(1./tao(8:10)) );
w = [wfi; wdvn; wdpos; wdKG; web; wdKA; wdb; wdposD; wdKD; wdvnS; wdposS];
%观测方程白噪声均方差
v = [0.01; 0.01; 0.01; 9/Re; 9/Re; 1;     0.01; 0.01; 0.01; 9/Re; 9/Re; 1];
%卡尔曼滤波时间
TKF = 1;     
% KF1:  X1_=Ft1*X1+w1        E(w1)=0, Cov(w1)=Qt1
%       Z1 =Ht1*X1+v1        E(v1)=0, Cov(v1)=Rk1
w1=w(1:25); v1=v(1:6); x10=x0(1:25);
Qt1=diag(w1.^2);     Rk1=diag(v1.^2);    Pk1=diag(x10.^2);    Xk1=zeros(25,1);
% KF2:  X2_=Ft2*X2+w2        E(w2)=0, Cov(w2)=Qt2
%       Z2 =Ht2*X2+v2        E(v2)=0, Cov(v2)=Rk2
w2=[w(1:21);w(26:31)]; v2=v(7:12); x20=[x0(1:21);x0(26:31)];
Qt2=diag(w2.^2);     Rk2=diag(v2.^2);    Pk2=diag(x20.^2);    Xk2=zeros(27,1);
%安装误差
dG = [0, 0, 0; 0, 0, 0; 0, 0, 0]*ppm;    dA = [0, 0, 0; 0, 0, 0; 0, 0, 0]*ppm;
KG = diag(ones(3,1)+dKG)+dG; KA = diag(ones(3,1)+dKA)+dA;

