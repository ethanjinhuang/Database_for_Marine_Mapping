function FBKfilter
%惯导-里程仪组合反馈校正卡尔曼滤波
    glb;
    % 捷联算法设置
    Tm  = 0.02;    
    fid = fopen('e:/ygm/vehicle/trace2_.bin','r');   %打开轨迹数据文件,读数据进行初始化
    [data, n] = fread(fid, 21, 'real*8');  %data [att, vn, pos, wm, vm] 
    qnb = Att2Quat(data(1:3));    vnm = data(4:6);    posm = data(7:9);
    qnb = QuatMul(Rv2Quat(-fi),qnb);  vnm = vnm + dvn;  posm = posm + dpos;  %初始误差
    % 里程仪算法设置
    TD = 0.1;
    vnD = data(4:6); posD = data(7:9);   
    dvnD = [0; 0; 0];
    vnD = vnD + dvnD; posD = posD + dposD;    %初始误差
    % GPS设置
    TS = 1.0;
    % 马尔可夫过程离散参数
    e_tao = exp([-1./tao(1:3)*Tm; -1./tao(4)*TD; -1./tao(5:10)*TS]);
    sqw = sqrt( Rv.*(ones(10,1)-exp(-2*[1./tao(1:3)*Tm; 1./tao(4)*TD; 1./tao(5:10)*TS])) );
    
    fout = fopen('e:/ygm/vehicle/kf5.bin','wb');
    wvm2 = data(10:15);
    for k=2:2:999*100  
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
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3卡尔曼滤波
                wb=data(10:12)/0.01; fb=data(13:15)/0.01; 
                [Ft, Ht] = getfh(Cnb, vnm, posm, wb+eb, fb+db, vnD, posD, tao);
                Ft1=Ft(1:25,:); Ft1(:,26:31)=[]; Ht1=Ht(1:6,:); Ht1(:,26:31)=[];
                Zk1=[vnm-vnD; posm-posD];
                if mod(k,300*100)==0        %反馈
                    [Xk1, Pk1] = kfilter(Ft1, Xk1, Qt1, Ht1, Zk1, Rk1, Pk1, TKF, 25);%卡尔曼滤波
                    qnb = QuatMul(Rv2Quat(Xk1(1:3)),qnb);
                    vnm = vnm-Xk1(4:6);
                    posm = posm-Xk1(7:9);
                    dKG = dKG-Xk1(10:12);
                    eb = eb-Xk1(13:15);
                    dKA = dKA-Xk1(16:18);
                    db = db-Xk1(19:21);
                    dposD = dposD-Xk1(22:24);
                    dKD = dKD-Xk1(25);
                    %Pk1=diag(x10.^2);    Xk1=zeros(25,1);
                elseif mod(k,301*100)==0        
                    [Xk1, Pk1] = kfilter1(Ft1, Xk1, Qt1, Ht1, Zk1, Rk1, Pk1, TKF, 25);%校正滤波
                else
                    [Xk1, Pk1] = kfilter(Ft1, Xk1, Qt1, Ht1, Zk1, Rk1, Pk1, TKF, 25);%卡尔曼滤波
                end
                
                E = Quat2Mat(qnb)*Att2Mat(data(1:3))'; fi = -[-E(2,3);E(1,3);-E(1,2)];
                err = [fi; vnm-data(4:6); posm-data(7:9); dKG; eb; dKA; db; posD-data(7:9); dKD];
                fwrite(fout, [err; Xk1], 'real*8');

                if mod(k,10*100) == 0
                    step=k/100,    %进度、时间显示
                end 
            end  %end 100
        end %10
    end   %end for
fclose(fid);
fclose(fout);

function [Xk, Pk] = kfilter1(Ft, Xk_1, Qt, Hk, Zk, Rk, Pk_1, Tkf, n)
%卡尔曼滤波
    In = eye(n);
    Fikk_1=In +Tkf*Ft +Tkf^2/2*Ft^2 +Tkf^3/6*Ft^3 +Tkf^4/24*Ft^4 +Tkf^5/120*Ft^5; 
    M1=Qt; M2=Ft*M1+(Ft*M1)'; M3=Ft*M2+(Ft*M2)'; M4=Ft*M3+(Ft*M3)'; M5=Ft*M4+(Ft*M4)';
    Qk=M1*Tkf +M2*Tkf^2/2 +M3*Tkf^3/6 +M4*Tkf^4/24 +M5*Tkf^5/120;
    
    Pkk_1=Fikk_1*Pk_1*Fikk_1'+Qk;                    %Qk
    Kk=Pkk_1*Hk'*(Hk*Pkk_1*Hk'+Rk)^-1;               %Rk
    Pk=(In-Kk*Hk)*Pkk_1;
    
    Xk=Kk*Zk;
%    Xkk_1=Fikk_1*Xk_1;    
%    Xk=Xkk_1+Kk*(Zk-Hk*Xkk_1);                       %Zk


