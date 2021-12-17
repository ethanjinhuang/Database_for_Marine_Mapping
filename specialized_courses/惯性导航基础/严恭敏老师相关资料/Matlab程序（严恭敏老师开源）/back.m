function federal_rk1
%验证 Ft 展开式对错
global Re e wie g0 ppm ug deg min sec hur dph   %全局变量

glb;
th = 0.02;               %龙格-库塔步长

    Gt = w/sqrt(th/2);       %在数值解法中,系统噪声离散化
    fid = fopen('e:/ygm/vehicle/trace4.bin','r');   %打开轨迹数据文件
    fout = fopen('e:/ygm/vehicle/federal_rk4.bin','wb');  
    [data2, n] = fread(fid, 21, 'real*8');  %data [att, vn, pos, wm, vm, wb, fb] 
    ww2 = randn(31,1);
    x1 = x0;
    Ft2 = getf(Att2Mat(data2(1:3)), data2(4:6), data2(7:9), data2(16:18), data2(19:21), data2(4:6), data2(7:9), tao);
    for k = 2:2:400*100
        [data, n] = fread(fid, 21+21, 'real*8'); 
        data0 = data2;      data1 = data(1:21);  data2 = data(22:42);
        ww = randn(31,2);
        ww0 = ww2;          ww1=ww(:,1);  ww2=ww(:,2);      
        %                      Cnb                  vn          pos         wb            fb            vnD         posD        
        k1 = getdx(x0,         Att2Mat(data0(1:3)), data0(4:6), data0(7:9), data0(16:18), data0(19:21), data0(4:6), data0(7:9), tao) + Gt.*ww0;
        k2 = getdx(x0+th/2*k1, Att2Mat(data1(1:3)), data1(4:6), data1(7:9), data1(16:18), data1(19:21), data1(4:6), data1(7:9), tao) + Gt.*ww1;
        k3 = getdx(x0+th/2*k2, Att2Mat(data1(1:3)), data1(4:6), data1(7:9), data1(16:18), data1(19:21), data1(4:6), data1(7:9), tao) + Gt.*ww1;
        k4 = getdx(x0+  th*k3, Att2Mat(data2(1:3)), data2(4:6), data2(7:9), data2(16:18), data2(19:21), data2(4:6), data2(7:9), tao) + Gt.*ww2;
        x0 = x0 + th/6*(k1+2*k2+2*k3+k4);
        
        Ft0 = Ft2;  
        Ft1 = getf(Att2Mat(data1(1:3)), data1(4:6), data1(7:9), data1(16:18), data1(19:21), data1(4:6), data1(7:9), tao);
        Ft2 = getf(Att2Mat(data2(1:3)), data2(4:6), data2(7:9), data2(16:18), data2(19:21), data2(4:6), data2(7:9), tao);
        m1 = Ft0* x1          + Gt.*ww0;
        m2 = Ft1*(x1+th/2*m1) + Gt.*ww1;
        m3 = Ft1*(x1+th/2*m2) + Gt.*ww1;
        m4 = Ft2*(x1+     m3) + Gt.*ww2;
        x1 = x1 + th/6*(m1+2*m2+2*m3+m4);

        if mod(k,100) == 0
            fwrite(fout, [x0;x1], 'real*8');
            if mod(k,10*100) == 0
              step=k/100,    %进度、时间显示
            end
        end 
    end   
fclose(fid);
fclose(fout);

function dx = getdx(x, Cnb, vn, pos, wb, fb, vnD, posD, tao)
global Re e wie  

%%%%%  fi    dvn   dpos   dKG   eb    dKA   db  22  dposD dKD  26  dvnS  dposS
    fi=x(1:3);  dvn=x(4:6); dpos=x(7:9);  dKG=x(10:12);  eb=x(13:15);  dKA=x(16:18); db=x(19:21);
    dposD=x(22:24); dKD=x(25); dvnS=x(26:28); dposS=x(29:31);
    
    %%
    sl = sin(pos(1)); cl = cos(pos(1)); tl = sl/cl; secl = 1/cl; secl2 = secl^2; sl2 = sl^2;
    RM = Re*(1-2*e+3*e*sl2); RN = Re*(1+e*sl2); 
    f_RMh = 1/(RM + pos(3)); f_RNh = 1/(RN + pos(3)); f_RMh2 = f_RMh^2; f_RNh2 = f_RNh^2;
    wnie = wie * [0; cl; sl];    wnen = [-vn(2)*f_RMh; vn(1)*f_RNh; vn(1)*f_RNh*tl];
    wnin = wnie + wnen;    

    dwnib = Cnb*(diag(wb)*dKG+eb);
    dfn = Cnb*(diag(fb)*dKA+db);
    dwnie = [0; -wie*sl*dpos(1); wie*cl*dpos(1)];
    dwnen = [-f_RMh*dvn(2)+f_RMh2*vn(2)*dpos(3); f_RNh*dvn(1)-f_RNh2*vn(1)*dpos(3); f_RNh*tl*dvn(1)-f_RNh2*vn(1)*tl*dpos(3)+f_RNh*vn(1)*secl2*dpos(3)];
    dwnin = dwnie+dwnen;
    
    dfi = -cross(wnin,fi)+dwnin-dwnib; 
    ddvn = cross(Cnb*fb,fi)-cross(2*wnie+wnen,dvn)+cross(vn,2*dwnie+dwnen)+dfn;
    ddLti = f_RMh*dvn(2)-f_RMh2*vn(2)*dpos(3);
    ddLgi = f_RNh*secl*dvn(1)+f_RNh*vn(1)*secl*tl*dpos(1)+f_RNh2*vn(1)*secl*dpos(3);
    ddHgt = dvn(3);
    ddKG = [0; 0; 0];
    deb = -1./tao(1:3).*eb;% + w(1:3);
    ddKA = [0; 0; 0];
    ddb = [0; 0; 0];
    
    %%
    slD = sin(posD(1)); clD = cos(posD(1)); tlD = slD/clD; seclD = 1/clD; slD2 = slD^2;
    RMD = Re*(1-2*e+3*e*slD2); RND = Re*(1+e*slD2); 
    f_RMhD = 1/(RM+posD(3)); f_RNhD = 1/(RN+posD(3)); f_RMhD2 = f_RMhD^2; f_RNhD2 = f_RNhD^2;
 
    dvnD = cross(vnD,fi)+vnD*dKD;
    
    ddLtiD = f_RMhD*dvnD(2)-f_RMhD2*vnD(2)*dposD(3);
    ddLgiD = f_RNhD*seclD*dvnD(1)+f_RNhD*vnD(1)*seclD*tlD*dposD(1)+f_RNhD2*vnD(3)*seclD*dposD(3);
    ddHgtD = dvnD(3);
    ddKD = -1/tao(4)*dKD;% + w(4);
    
    %%
    ddvnS = -1./tao(5:7).*dvnS;% + w(5:10);
    ddposS = -1./tao(8:10).*dposS;% + w(5:10);
    
    dx = [dfi; ddvn; ddLti;ddLgi;ddHgt; ddKG; deb; ddKA; ddb;    ddLtiD;ddLgiD;ddHgtD; ddKD;   ddvnS; ddposS] ;
    
    
