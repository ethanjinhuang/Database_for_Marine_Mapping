function Alignment
%卡尔曼滤波初始精对准（降维），缺乏机动效果不好！
global Re e wie g0 ppm ug deg min sec hur dph  tao %全局变量
glb;
Tm  = 0.02;    
TKF = 1;     

%陀螺漂移：          eb_=-1/tao(1:3)*eb+web;       
tao = [360; 360; 360];
Rv = [(0.01*dph)^2; (0.01*dph)^2; (0.01*dph)^2];
% 马尔可夫过程离散参数
e_tao = exp(-1./tao(1:3)*Tm);
sqw = sqrt( Rv.*(ones(3,1)-exp(-2*1./tao*Tm)) );
%系统状态初值设置
fi = [2 ; 3; 7]*min;     dvn = [.1; .1; .1];           dpos = [20/Re; 20/Re; 10];	
dKG = [100; 99; 98]*ppm;    eb = [0.01; 0.01; 0.01]*dph;  
dKA = [100; 99; 98]*ppm;	 db = [100; 99; 98]*ug;
x0 = [fi; dvn; dpos; dKG; eb; dKA; db];
%系统方程白噪声均方差
wfi = [.01; 0.01; 0.01]*min; wdvn = [0.01; 0.01; 0.01];    wdpos = [5/Re; 5/Re; 5];
wdKG = [0; 0; 0];            web = sqrt( 2*Rv.*(1./tao) );
wdKA = [0; 0; 0];            wdb = [0; 0; 0];
w = [wfi; wdvn; wdpos; wdKG; web; wdKA; wdb];
%观测方程白噪声均方差
v = [0.01; 0.01; 0.01; 5/Re; 5/Re; 5];
% KF1:  X1_=Ft1*X1+w1        E(w1)=0, Cov(w1)=Qt1
%       Z1 =Ht1*X1+v1        E(v1)=0, Cov(v1)=Rk1
Qt=diag(w(1:9).^2);     Rk=diag(v.^2);    Pk=diag(x0(1:9).^2);    Xk=zeros(9,1);

    % 轨迹设置
    att = pi/180*[0; 0; 0];	 %[pitch roll azimuth]  
	vb = [0; 0; 0];          %[0 vby 0]       
	pos = [34*deg+14.76289014*min; 108*deg+54.57983*min; 380];  %[latitude longitude height]
    %rk(1:3) att(pitch roll azimuth);   rk(2:6) vn(vnE vnN vnU);   rk(7:9) pos(lti lgi hgt)
    %rk(10:12) wm;  rk(13:15) vm;         
    rk = [att; Att2Mat(att)*vb; pos; zeros(6,1)];
    wat=zeros(6,3);
    rk = track(rk, wat, Tm/2);    
    dwvm=rk(10:15);
    
    % 捷联算法设置
    qnb = Att2Quat(att);    vnm = Att2Mat(att)*vb;    posm = pos;
    qnb = QuatMul(Rv2Quat(fi), qnb);  vnm = vnm + dvn;  posm = posm + dpos;  %初始误差

    %fout = fopen('e:/ygm/vehicle/kf4.bin','wb');
    m=1;
    for k=2 :2: 1000*100  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1惯导更新
        ns = randn(3,2);
        eb = e_tao.* eb  + sqw.*ns(:,1);
        dwvm1 = (ones(6,1)+[dKG;dKA]).*dwvm + [eb; db]*Tm/2;
        eb = e_tao.* eb  + sqw.*ns(:,2);
        dwvm2 = (ones(6,1)+[dKG;dKA]).*dwvm + [eb; db]*Tm/2;        %刻度系数、漂移
        [qnb, vnm, posm] = sins(qnb, vnm, posm, dwvm1, dwvm2, Tm);
        if mod(k,100)==0   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2卡尔曼滤波
            Cnb = Quat2Mat(qnb);
            [Ft, Ht]= getfh(Cnb, vnm, posm, dwvm(1:3)/Tm/2+eb, dwvm(4:6)/Tm/2+db);
            Zk=[vnm; posm-pos];
            [Xk, Pk] = kfilter(Ft, Xk, Qt, Ht, Zk, Rk, Pk, TKF, 9);
            E = Quat2Mat(qnb)*Att2Mat(att)'; fi = [-E(2,3);E(1,3);-E(1,2)];
            %fwrite(fout, [Xk; x], 'real*8');
            x(m,:) = [fi; vnm; posm-pos; dKG; eb; dKA; db]';
            y(m,:) = Xk';
            m=m+1;
            if mod(m,10)==0
                step=m/10,
            end
%            if mod(k,400*100) == 0
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3反馈校正
%                qnb = QuatMul(Rv2Quat(Xk(1:3)),qnb);
 %               vnm = vnm-Xk(4:6);
  %              posm = posm-Xk(7:9);
   %             Pk=diag(x0(1:9).^2);    Xk=zeros(9,1);
   %        end %end 100*100
        end
    end
    %fclose(fout)
figure;
xk=x;
subplot(3,3,1);  plot(1/min*xk(:,1:3));        ylabel('fx fy fz(min)');       
subplot(3,3,2);  plot(xk(:,4:6));              ylabel('dVnx dVny dVnz(m/s)');
subplot(3,3,3);  plot(Re*xk(:,7:8));           ylabel('dLti dLgi(m)');      
subplot(3,3,4);  plot(xk(:,9));                ylabel('dH(m)');
subplot(3,3,5);  plot(1/ppm*xk(:,10:12));      ylabel('dKGx dKGy dKGz(ppm)'); 
subplot(3,3,6);  plot(1/dph*xk(:,13:15));      ylabel('ebx eby ebz(deg/h)');    
subplot(3,3,7);  plot(1/ppm*xk(:,16:18));      ylabel('dKAx dKAy dKAz(ppm)');
subplot(3,3,8);  plot(1/ug*xk(:,19:21));       ylabel('dbx dby dbz(ug)');
figure;
xk=y;
subplot(3,3,1);  plot(1/min*xk(:,1:3));        ylabel('fx fy fz(min)');       
subplot(3,3,2);  plot(xk(:,4:6));              ylabel('dVnx dVny dVnz(m/s)');
subplot(3,3,3);  plot(Re*xk(:,7:8));           ylabel('dLti dLgi(m)');      
subplot(3,3,4);  plot(xk(:,9));                ylabel('dH(m)');
figure;
xk=x(:,1:9)-y;
subplot(3,3,1);  plot(1/min*xk(:,1:3));        ylabel('fx fy fz(min)');       
subplot(3,3,2);  plot(xk(:,4:6));              ylabel('dVnx dVny dVnz(m/s)');
subplot(3,3,3);  plot(Re*xk(:,7:8));           ylabel('dLti dLgi(m)');      
subplot(3,3,4);  plot(xk(:,9));                ylabel('dH(m)');

function [Ft, Ht] = getfh(Cnb, vn, pos, wb, fb)
%求系统矩阵Ft
global Re e g0 wie tao Rv
	sl = sin(pos(1)); cl = cos(pos(1)); tl = sl/cl; secl = 1/cl; secl2 = secl^2; sl2 = sl^2;
    RM = Re*(1-2*e+3*e*sl2); RN = Re*(1+e*sl2); 
    f_RMh = 1/(RM + pos(3)); f_RNh = 1/(RN + pos(3)); f_RMh2 = f_RMh^2; f_RNh2 = f_RNh^2;
    wnie = wie * [0; cl; sl];    wnen = [-vn(2)*f_RMh; vn(1)*f_RNh; vn(1)*f_RNh*tl];
    wnin = wnie + wnen;    
    %%%
    M1 = [0, 0, 0; -wie*sl, 0, 0; wie*cl, 0, 0];
    M2 = [0, -f_RMh, 0; f_RNh, 0, 0; f_RNh*tl, 0, 0];
    M3 = [0, 0, vn(2)*f_RMh2; 0, 0, -vn(1)*f_RNh2; vn(1)*secl2*f_RNh, 0, -vn(1)*tl*f_RNh2];
    M13 = M1+M3;
    M4 = Asym(vn)*M2 - Asym(2*wnie+wnen);                                                                         
    M5 = Asym(vn)*(2*M1+M3);
    M6 = [0, f_RMh, 0; secl*f_RNh, 0, 0; 0, 0, 1];
    M7 = [0, 0, -vn(2)*f_RMh2; vn(1)*secl*tl*f_RNh, 0, vn(1)*secl*f_RNh2; 0, 0, 0];
    aG = diag(-1./tao);
    S1 = Asym(wnin); S2 = Cnb*diag(wb); S3 = Asym(Cnb*fb); S4 = Cnb*diag(fb); 

    o3 = zeros(3,3); I3=eye(3);
    %%%%%  fi    dvn   dpos       
    Ft = [ -S1   M2    M13     ;
           S3    M4    M5         ;
           o3    M6    M7    ];
                    
    Ht = [ o3    I3    o3 ;
           o3    o3    I3 ];



        