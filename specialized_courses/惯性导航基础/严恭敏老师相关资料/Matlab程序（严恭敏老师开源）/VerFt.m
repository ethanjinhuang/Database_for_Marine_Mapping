function VerFt
%验证 Ft 展开式对错
global Re e wie g0 ppm ug deg min sec hur dph   %全局变量

glb;

    for k = 1:1:20
        x=10e10*randn(31,1); 
        Cnb=randn(3,3);  vn=randn(3,1);  pos=randn(3,1); wb=randn(3,1); fb=randn(3,1);  vnD=randn(3,1); posD=randn(3,1);
        %%%%%  fi    dvn   dpos   dKG   eb    dKA   db  22  dposD dKD  26  dvnS  dposS
     %  Ft = [ -S1   M2    M13    -S2   -Cnb  o3    o3      o3    o31      o3    o3;
        dx0 = getdx(x, Cnb, vn, pos, wb, fb, vnD, posD, tao);   %Ft 不展开
        Ft = getf(Cnb, vn, pos, wb, fb, vnD, posD, tao);        %Ft 展开
        dx1 = Ft*x;  
        err(k,:)=(dx1-dx0)'./(abs(dx0)+10e-13*ones(31,1))';        %计算(运算顺序先后不同)存在误差，若相对误差10e-10，则在误差范围之内(?)
    end   
figure;
subplot(4,4,1);  plot(1/min*err(:,1:3));        ylabel('fx fy fz(min)');       
subplot(4,4,2);  plot(err(:,4:6));              ylabel('dVnx dVny dVnz(m/s)');
subplot(4,4,3);  plot(Re*err(:,7:8));        ylabel('dLti dLgi(m)');      
subplot(4,4,4);  plot(err(:,9));                ylabel('dH(m)');
subplot(4,4,5);  plot(1/ppm*err(:,10:12));      ylabel('dKGx dKGy dKGz(ppm)'); 
subplot(4,4,6);  plot(1/dph*err(:,13:15));      ylabel('ebx eby ebz(deg/h)');    
subplot(4,4,7);  plot(1/ppm*err(:,16:18));      ylabel('dKAx dKAy dKAz(ppm)');
subplot(4,4,8);  plot(1/ug*err(:,19:21));       ylabel('dbx dby dbz(ug)');
subplot(4,4,9);  plot(Re*err(:,22:23));      ylabel('dLtiD dLgiD(m)');
subplot(4,4,10); plot(err(:,24));               ylabel('dHD(m)');
subplot(4,4,11); plot(err(:,25));               ylabel('dKD');            
subplot(4,4,13); plot(err(:,26:28));            ylabel('dVnxS dVnyS dVnzS(m/s)');
subplot(4,4,14); plot(Re*err(:,29:30));      ylabel('dLtiS dLgiS(m)');
subplot(4,4,15); plot(err(:,31));               ylabel('dHS(m)');           

