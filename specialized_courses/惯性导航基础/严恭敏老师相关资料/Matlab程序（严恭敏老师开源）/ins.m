function ins
global Re e g0 wie 
Re = 6378160;   e = 1/298.3;   wie = 7.2921151467e-5;   g0 = 9.7803267714;
Tm  = 0.02;                     %姿态更新周期20ms

    fid = fopen('e:/ygm/vehicle/trace3_.bin','r');   %打开轨迹数据文件,读数据进行初始化
%    status = fseek(fid,15*8*201*100,0);
    [data, n] = fread(fid, 15, 'real*8');         %data=[att, vn, pos, am, vm, wb,fb];
    qnb = Att2Quat(data(1:3)); vnm = data(4:6); posm = data(7:9);

    fout = fopen('e:/ygm/vehicle/trsins1.bin','wb');  
    wvm2 = data(10:15);
    for k=2:2:500*100  
        %读数据,获得角增量比力增量
        [data, n] = fread(fid, 15+15, 'real*8');  
        wvm0 = wvm2;    wvm1 = data(10:15);     wvm2 = data((15+10):(15+15));
        %捷联解算
        [qnb, vnm, posm] = sins(qnb, vnm, posm, wvm1-wvm0, wvm2-wvm1, Tm);
        if mod(k,10*100)==0
            fwrite(fout, [data(1:9); pi/180*Quat2Att(qnb); vnm; posm], 'real*8');
            if mod(k,1000)==0
                step = k/1000,                        %进度显示
            end
        end
    end
fclose(fid);
fclose(fout);
    
