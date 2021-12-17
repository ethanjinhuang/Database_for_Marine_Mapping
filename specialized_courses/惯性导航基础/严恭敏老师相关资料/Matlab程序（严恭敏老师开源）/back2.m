function federal_rk
%状态方程正向用微分方程法解算
%global Re e wie g0 ppm ug deg min sec hur dph PKG PKA  %全局变量
glb;

    tt = 0.02;               %龙格-库塔步长
    Gt = w/sqrt(tt/2);       %在数值解法中,系统噪声离散化
    fid = fopen('e:/ygm/vehicle/trace4.bin','r');   %打开轨迹数据文件
    fout = fopen('e:/ygm/vehicle/federal_rk4.bin','wb');
    [data2, n] = fread(fid, 21, 'real*8');  %data [att, vn, pos, wm, vm, wb, fb] 
    ww2 = zeros(31,1);
    for k = 2:2:500*100
        [data, n] = fread(fid, 21+21, 'real*8');  %data [att, vn, pos, wm, vm, wb, fb] 
        ww = randn(31,2);
        %          Cnb                 vn         pos        wb           fb           vnD  posD
        data0=data2;        data1 = data(1:21); data2 = data(22:42);
        ww0 = ww2;          ww1=ww(:,1);  ww2=ww(:,2);      
        k1 = getdx(x0,         Att2Mat(data0(1:3)), data0(4:6), data0(7:9), data0(16:18), data0(19:21), data0(4:6), data0(7:9), tao) + Gt.*ww0;
        k2 = getdx(x0+tt/2*k1, Att2Mat(data1(1:3)), data1(4:6), data1(7:9), data1(16:18), data1(19:21), data1(4:6), data1(7:9), tao) + Gt.*ww1;
        k3 = getdx(x0+tt/2*k2, Att2Mat(data1(1:3)), data1(4:6), data1(7:9), data1(16:18), data1(19:21), data1(4:6), data1(7:9), tao) + Gt.*ww1;
        k4 = getdx(x0+  tt*k3, Att2Mat(data2(1:3)), data2(4:6), data2(7:9), data2(16:18), data2(19:21), data2(4:6), data2(7:9), tao) + Gt.*ww2;
        x0 = x0 + tt/6*(k1+2*k2+2*k3+k4);

        if mod(k,100) == 0
            Ht = geth(data2(4:6));         zk = Ht*x0 + v.*randn(12,1); 
            Ft = getf(Att2Mat(data2(1:3)), data2(4:6), data2(7:9), data2(16:18), data2(19:21), data2(4:6), data2(7:9), tao);
            % kf1
            Ft1=Ft(1:25,:); Ft1(:,26:31)=[]; Zk1=zk(1:6); Ht1=Ht(1:6,:); Ht1(:,26:31)=[];
            [Xk1, Pk1] = kfilter(Ft1, Xk1, Qt1, Ht1, Zk1, Rk1, Pk1, TKF, 25);
            % kf2
            Ft2=Ft; Ft2(22:25,:)=[]; Ft2(:,22:25)=[]; Zk2=zk(7:12); Ht2=Ht(7:12,:); Ht2(:,22:25)=[];
            [Xk2, Pk2] = kfilter(Ft2, Xk2, Qt2, Ht2, Zk2, Rk2, Pk2, TKF, 27);
            
            fwrite(fout, [x0;zk;Xk1;Xk2], 'real*8');

            if mod(k,10*100) == 0
              step=k/100,    %进度、时间显示
            end
        end 
    end   
fclose(fid);
fclose(fout);

