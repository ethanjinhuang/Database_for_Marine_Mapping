function VerSINS
%直接解误差微分方程 和 捷联算法误差之间的差别(还可验证失准角书写方式)
%global Re e wie g0 ppm ug deg min sec hur dph PKG PKA  %全局变量
glb;

    Tm = 0.02;               %龙格-库塔步长
    Gt = w/sqrt(Tm/2);       %在数值解法中,系统噪声离散化
    fid = fopen('e:/ygm/vehicle/trace4.bin','r');   %打开轨迹数据文件
    fout = fopen('e:/ygm/vehicle/VerSINS.bin','wb');
    [data2, n] = fread(fid, 21, 'real*8');  %data [att, vn, pos, wm, vm, wb, fb] 
    qnb = Att2Quat(data2(1:3));    vnm = data2(4:6);    posm = data2(7:9);                  wvm2 = data2(10:15);
    qnb = QuatMul(Rv2Quat(-fi), qnb);  vnm = vnm + dvn;  posm = posm + dpos;  %初始误差
    ww2 = zeros(31,1);

    for k = 2:2:900*100
        [data, n] = fread(fid, 21+21, 'real*8');
        ww = randn(31,2);
        %          Cnb                 vn         pos        wb           fb           vnD  posD
        data0=data2;        data1 = data(1:21); data2 = data(22:42);
        ww0 = ww2;          ww1=ww(:,1);  ww2=ww(:,2);      
        k1 = getdx(x0,         Att2Mat(data0(1:3)), data0(4:6), data0(7:9), data0(16:18), data0(19:21), data0(4:6), data0(7:9), tao) + Gt.*ww0;
        k2 = getdx(x0+Tm/2*k1, Att2Mat(data1(1:3)), data1(4:6), data1(7:9), data1(16:18), data1(19:21), data1(4:6), data1(7:9), tao) + Gt.*ww1;
        k3 = getdx(x0+Tm/2*k2, Att2Mat(data1(1:3)), data1(4:6), data1(7:9), data1(16:18), data1(19:21), data1(4:6), data1(7:9), tao) + Gt.*ww1;
        k4 = getdx(x0+  Tm*k3, Att2Mat(data2(1:3)), data2(4:6), data2(7:9), data2(16:18), data2(19:21), data2(4:6), data2(7:9), tao) + Gt.*ww2;
        x0 = x0 + Tm/6*(k1+2*k2+2*k3+k4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1惯导更新
        wvm0 = wvm2; wvm1 = data(10:15); wvm2 = data((21+10):(21+15));
        dwvm1 = wvm1-wvm0; dwvm2 = wvm2-wvm1;
%        eb = e_tao(1:3).* eb  + sqw(1:3).*ww1(13:15); 
        eb = x0(13:15);     %%%%%  fi    dvn   dpos   dKG   eb    dKA   db  22  dposD dKD 
        dwvm1 = [KG*dwvm1(1:3); KA*dwvm1(4:6)] + [x0(13:15); db]*Tm/2;
%        eb = e_tao(1:3).* eb  + sqw(1:3).*ww2(13:15);
        dwvm2 = [KG*dwvm1(1:3); KA*dwvm1(4:6)] + [x0(13:15); db]*Tm/2;        %刻度系数、漂移
        [qnb, vnm, posm] = sins(qnb, vnm, posm, dwvm1, dwvm2, Tm);
        if mod(k,100)==0
            E = Quat2Mat(qnb)*Att2Mat(data(1:3))'; fi = -[-E(2,3);E(1,3);-E(1,2)];
            x1 = [fi; vnm-data(4:6); posm-data(7:9); dKG; eb; dKA; db];
            fwrite(fout, [x0(1:21); x1], 'real*8');
            if mod(k,1000) == 0
                step=k/100,    %进度、时间显示
            end %end 1000
        end
    end   
fclose(fid);
fclose(fout);

