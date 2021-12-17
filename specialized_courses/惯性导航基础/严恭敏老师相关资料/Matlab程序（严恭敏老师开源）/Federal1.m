function federal
%惯导里程仪卡尔曼滤波
global Re e wie g0 ppm ug deg min sec hur dph  tao %全局变量
%马尔可夫过程参数
%陀螺漂移：          eb_=-1/tao(1:3)*eb+web;       
%里程仪刻度系数误差：dKD_=-1/tao(4)+wdKD;
tao = [3600; 3600; 3600;   100];
Rv = [(0.01*dph)^2; (0.01*dph)^2; (0.01*dph)^2;  0.003^2];
%系统状态初值设置
fi = [0.7 ; 0.7; 5]*min;     dvn = [.01; .01; .01];           dpos = [20/Re; 20/Re; 10];	
dKG = [100; 99; 98]*ppm;    eb = [0.01; 0.01; 0.01]*dph;  
dKA = [100; 99; 98]*ppm;	   db = [100; 99; 98]*ug;
dposD = [20/Re; 20/Re; 10];	dKD = sqrt(Rv(4));        
x0 = [fi; dvn; dpos; dKG; eb; dKA; db; dposD; dKD];
%系统方程白噪声均方差
wfi = [0.1; 0.1; 0.1]*min; wdvn = [0.01; 0.01; 0.01];    wdpos = [5/Re; 5/Re; 5];
wdKG = [0; 0; 0];            web = sqrt( 2*Rv(1:3).*(1./tao(1:3)) );
wdKA = [0; 0; 0];            wdb = [0; 0; 0];
wdposD = [5/Re; 5/Re; 5];    wdKD = sqrt(2*Rv(4)*1/tao(4));
w = [wfi; wdvn; wdpos; wdKG; web; wdKA; wdb; wdposD; wdKD];
%观测方程白噪声均方差
v = [0.01; 0.01; 0.01; 5/Re; 5/Re; 5];
%卡尔曼滤波时间
TKF = 1;     
% KF1:  X1_=Ft1*X1+w1        E(w1)=0, Cov(w1)=Qt1
%       Z1 =Ht1*X1+v1        E(v1)=0, Cov(v1)=Rk1
Qt=diag(w.^2);     Rk=diag(v.^2);    Pk=diag(x0.^2);    Xk=zeros(25,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 初始设置
    dG = [0, 0, 0; 0, 0, 0; 0, 0, 0]*ppm;    dA = [0, 0, 0; 0, 0, 0; 0, 0, 0]*ppm;
    KG = diag(ones(3,1)+dKG)+dG; KA = diag(ones(3,1)+dKA)+dA;
    dvnD = [0; 0; 0];
    % 捷联算法设置
    Tm  = 0.02;    
    fid = fopen('e:/ygm/vehicle/trace2_.bin','r');   %打开轨迹数据文件,读数据进行初始化
%    fseek(fid, 21*8*1000*100, 'bof');
    [data, n] = fread(fid, 21, 'real*8');  %data [att, vn, pos, wm, vm, wb, fb] 
    qnb = Att2Quat(data(1:3));    vnm = data(4:6);    posm = data(7:9);
    qnb = QuatMul(Rv2Quat(fi),qnb);  vnm = vnm + dvn;  posm = posm + dpos; 
%    E = Quat2Mat(qnb)*Att2Mat(data(1:3))'; fii = [-E(2,3);E(1,3);-E(1,2)];    qnb = QuatMul(Rv2Quat(-fii),qnb);
    % 里程仪算法设置
    TD = 0.1;
    vnD = data(4:6); posD = data(7:9);
    vnD = vnD + dvnD; posD = posD + dposD;    %初始误差
    % 马尔可夫过程离散参数
    e_tao = exp([-1./tao(1:3)*Tm; -1./tao(4)*TD]);
    sqw = sqrt( Rv.*(ones(4,1)-exp(-2*[1./tao(1:3)*Tm; 1./tao(4)*TD])) );
    
    fout = fopen('e:/ygm/vehicle/kf4.bin','wb');
    wvm2 = data(10:15);
%    for k=2:2:1830*100  
    for k=2 :2: 800*100  
        %读数据,获得角增量比力增量
        [data, n] = fread(fid, 21+21, 'real*8');  
        wvm0 = wvm2; wvm1 = data(10:15); wvm2 = data((21+10):(21+15));
        dwvm1 = wvm1-wvm0; dwvm2 = wvm2-wvm1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1惯导更新
        ns = randn(3,2);
        eb = e_tao(1:3).* eb  + sqw(1:3).*ns(:,1);
        dwvm1 = [KG*dwvm1(1:3); KA*dwvm1(4:6)] + [eb; db]*Tm/2;
        eb = e_tao(1:3).* eb  + sqw(1:3).*ns(:,2);
        dwvm2 = [KG*dwvm1(1:3); KA*dwvm1(4:6)] + [eb; db]*Tm/2;        %刻度系数、漂移
        [qnb, vnm, posm] = sins(qnb, vnm, posm, dwvm1, dwvm2, Tm);
        if mod(k, 10)==0
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2里程仪更新       
            dKD = e_tao(4)*dKD + sqw(4)*randn(1,1);
            vD = (1+dKD) * sqrt(sum(data(4:6).^2));
            Cnb = Quat2Mat(qnb);
            vnD = Cnb * [0; vD; 0];
            sl = sin(posD(1)); cl = cos(posD(1)); 
            RM = Re*(1-2*e+3*e*sl^2); RN = Re*(1+e*sl^2); RMhD = RM + posD(3); RNhD = RN + posD(3);
            posD = posD + TD*[vnD(2)/RMhD; vnD(1)/(RNhD*cl); vnD(3)];
            if mod(k,100)==0      
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3卡尔曼滤波
                [Ft, Ht] = getfh(Cnb, vnm, posm, data(16:18)+eb, data(19:21)+db, vnD, posD);
                Zk = [vnm-vnD; posm-posD];
                [Xk, Pk] = kfilter(Ft, Xk, Qt, Ht, Zk, Rk, Pk, TKF, 25);
                %%%%%  fi    dvn   dpos   dKG   eb    dKA   db  22  dposD dKD 
                E = Quat2Mat(qnb)*Att2Mat(data(1:3))'; fi = [-E(2,3);E(1,3);-E(1,2)];
                x = [fi; vnm-data(4:6); posm-data(7:9); dKG; eb; dKA; db; vnD-data(4:6); posD-data(7:9); dKD];
                fwrite(fout, [x; Zk; Xk], 'real*8');
                if mod(k,200*100) == 0
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4反馈校正 
                    qnb = QuatMul(Rv2Quat(Xk(1:3)),qnb);
                    vnm = vnm-Xk(4:6);
                    posm = posm-Xk(7:9);
                    posD = posD-Xk(22:24);
                    Pk=diag(x0.^2);    Xk=zeros(25,1);
                end %end 100*100
                if mod(k,10*100) == 0
                    step=k/100,    %进度、时间显示
                end 
            end %end 100
        end %end 10
    end    %end for
fclose(fid);
fclose(fout);

function [Ft, Ht] = getfh(Cnb, vn, pos, wb, fb, vnD, posD)
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
    aG = diag(-1./tao(1:3));
    S1 = Asym(wnin); S2 = Cnb*diag(wb); S3 = Asym(Cnb*fb); S4 = Cnb*diag(fb); 
    %%%
    slD = sin(posD(1)); clD = cos(posD(1)); tlD = slD/clD; seclD = 1/clD; slD2 = slD^2;
    RMD = Re*(1-2*e+3*e*slD2); RND = Re*(1+e*slD2); 
    f_RMhD = 1/(RM+posD(3)); f_RNhD = 1/(RN+posD(3)); f_RMhD2 = f_RMhD^2; f_RNhD2 = f_RNhD^2;
    MD6 = [0, f_RMhD, 0; seclD*f_RNhD, 0, 0; 0, 0, 1];
    MD7 = [0, 0, -vnD(2)*f_RMhD2; vnD(1)*seclD*tlD*f_RNhD, 0, vnD(1)*seclD*f_RNhD2; 0, 0, 0];
    SD1 = MD6*Asym(vnD); SD2 = MD6*vnD; aD = -1/tao(4);
        
    o13 = [0, 0, 0]; o31 = [0; 0; 0]; o3 = [0, 0, 0; 0, 0, 0; 0, 0, 0];  I3=[1,0,0; 0,1,0; 0,0,1];

    %%%%%  fi    dvn   dpos   dKG   eb    dKA   db  22  dposD dKD  
    Ft = [ -S1   M2    M13    -S2   -Cnb  o3    o3      o3    o31      ;
           S3    M4    M5     o3    o3    S4    Cnb     o3    o31      ;
           o3    M6    M7     o3    o3    o3    o3      o3    o31      ;
           o3    o3    o3     o3    o3    o3    o3      o3    o31      ;
           o3    o3    o3     o3    aG    o3    o3      o3    o31      ;
           o3    o3    o3     o3    o3    o3    o3      o3    o31      ;
           o3    o3    o3     o3    o3    o3    o3      o3    o31      ; 
       
           SD1   o3    o3     o3    o3    o3    o3      MD7   SD2      ; 
           o13   o13   o13    o13   o13   o13   o13     o13   aD        ]; 
           
    SD1 = Asym(vnD);
    %%%%%  fi    dvn   dpos   dKG   eb    dKA   db %22  dposD dKD %26
    Ht = [ -SD1  I3    o3     o3    o3    o3    o3      o3    -vnD     ;
           o3    o3    I3     o3    o3    o3    o3      -I3   o31       ]; 