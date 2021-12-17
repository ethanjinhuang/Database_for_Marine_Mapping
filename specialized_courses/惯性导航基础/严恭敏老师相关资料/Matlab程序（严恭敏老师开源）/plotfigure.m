function plottype
global Re e wie g0 ppm ug deg min sec hur dph   %全局变量

fid = fopen('e:/ygm/vehicle/kf4.bin','r');
for k=1:1:1800
    [data, n] = fread(fid, 28+6+25, 'real*8');   %fwrite(fout, [xk;zk1;Xk1;Xk2], 'real*8');
    xk(k,:) = data(1:28)';
    zk1(k,:) = data(29:34)';   
    xk1(k,:) = data(35:59)';
end
figure;%解算值 
subplot(4,4,1);  plot(1/min*xk(:,1:3));        ylabel('fx fy fz(min)');       
subplot(4,4,2);  plot(xk(:,4:6));              ylabel('dVnx dVny dVnz(m/s)');
subplot(4,4,3);  plot(Re*xk(:,7:8));        ylabel('dLti dLgi(m)');      
subplot(4,4,4);  plot(xk(:,9));                ylabel('dH(m)');
subplot(4,4,5);  plot(1/ppm*xk(:,10:12));      ylabel('dKGx dKGy dKGz(ppm)'); 
subplot(4,4,6);  plot(1/dph*xk(:,13:15));      ylabel('ebx eby ebz(deg/h)');    
subplot(4,4,7);  plot(1/ppm*xk(:,16:18));      ylabel('dKAx dKAy dKAz(ppm)');
subplot(4,4,8);  plot(1/ug*xk(:,19:21));       ylabel('dbx dby dbz(ug)');
subplot(4,4,9);  plot(xk(:,22:24));      ylabel('dVnDx dVnDy dVnDz(m/s)');
subplot(4,4,10);  plot(Re*xk(:,25:26));      ylabel('dLtiD dLgiD(m)');
subplot(4,4,11); plot(xk(:,27));               ylabel('dHD(m)');
subplot(4,4,12); plot(xk(:,28));               ylabel('dKD');            
figure;%观测值 
subplot(2,3,1);  plot(zk1(:,1:3));              ylabel('z11 z12 z13(m/2)');       
subplot(2,3,2);  plot(Re*zk1(:,4:5));        ylabel('z14 z15(m)');
subplot(2,3,3);  plot(zk1(:,6));                ylabel('z16(m)');
figure;%滤波值误差1
xk(:,22:24)=[];
%xk1(:,1:25) = xk1(:,1:25) - xk(:,1:25);
subplot(4,4,1);  plot(1/min*xk1(:,1:3));        ylabel('fx fy fz(min)');       
subplot(4,4,2);  plot(xk1(:,4:6));              ylabel('dVnx dVny dVnz(m/s)');
subplot(4,4,3);  plot(Re*xk1(:,7:8));        ylabel('dLti dLgi(m)');      
subplot(4,4,4);  plot(xk1(:,9));                ylabel('dH(m)');
subplot(4,4,5);  plot(1/ppm*xk1(:,10:12));      ylabel('dKGx dKGy dKGz(ppm)'); 
subplot(4,4,6);  plot(1/dph*xk1(:,13:15));      ylabel('ebx eby ebz(deg/h)');    
subplot(4,4,7);  plot(1/ppm*xk1(:,16:18));      ylabel('dKAx dKAy dKAz(ppm)');
subplot(4,4,8);  plot(1/ug*xk1(:,19:21));       ylabel('dbx dby dbz(ug)');
subplot(4,4,9);  plot(Re*xk1(:,22:23));      ylabel('dLtiD dLgiD(m)');
subplot(4,4,10); plot(xk1(:,24));               ylabel('dHD(m)');
subplot(4,4,11); plot(xk1(:,25));               ylabel('dKD');      
return


fid = fopen('e:/ygm/vehicle/kf6.bin','r');
for k=1:1:1820
    [data, n] = fread(fid, 28+6+12, 'real*8');   %fwrite(fout, [xk;zk1;Xk1;Xk2], 'real*8');
    xk(k,:) = data(1:28)';
    zk1(k,:) = data(29:34)';   
    xk1(k,:) = data(35:46)';
end
figure;%解算值 
subplot(4,4,1);  plot(1/min*xk(:,1:3));        ylabel('fx fy fz(min)');       
subplot(4,4,2);  plot(xk(:,4:6));              ylabel('dVnx dVny dVnz(m/s)');
subplot(4,4,3);  plot(Re*xk(:,7:8));        ylabel('dLti dLgi(m)');      
subplot(4,4,4);  plot(xk(:,9));                ylabel('dH(m)');
subplot(4,4,5);  plot(1/ppm*xk(:,10:12));      ylabel('dKGx dKGy dKGz(ppm)'); 
subplot(4,4,6);  plot(1/dph*xk(:,13:15));      ylabel('ebx eby ebz(deg/h)');    
subplot(4,4,7);  plot(1/ppm*xk(:,16:18));      ylabel('dKAx dKAy dKAz(ppm)');
subplot(4,4,8);  plot(1/ug*xk(:,19:21));       ylabel('dbx dby dbz(ug)');
subplot(4,4,9);  plot(xk(:,22:24));      ylabel('dVnDx dVnDy dVnDz(m/s)');
subplot(4,4,10);  plot(Re*xk(:,25:26));      ylabel('dLtiD dLgiD(m)');
subplot(4,4,11); plot(xk(:,27));               ylabel('dHD(m)');
subplot(4,4,12); plot(xk(:,28));               ylabel('dKD');            
figure;%观测值 
subplot(2,3,1);  plot(zk1(:,1:3));              ylabel('z11 z12 z13(m/2)');       
subplot(2,3,2);  plot(Re*zk1(:,4:5));        ylabel('z14 z15(m)');
subplot(2,3,3);  plot(zk1(:,6));                ylabel('z16(m)');
figure;%滤波值误差1
xk(:,28)=[]; xk(:,10:24)=[];
xk1 = xk1 - xk;
subplot(4,4,1);  plot(1/min*xk1(:,1:3));        ylabel('fx fy fz(min)');       
subplot(4,4,2);  plot(xk1(:,4:6));              ylabel('dVnx dVny dVnz(m/s)');
subplot(4,4,3);  plot(Re*xk1(:,7:8));        ylabel('dLti dLgi(m)');      
subplot(4,4,4);  plot(xk1(:,9));                ylabel('dH(m)');

subplot(4,4,9);  plot(Re*xk1(:,10:11));      ylabel('dLtiD dLgiD(m)');
subplot(4,4,10); plot(xk1(:,12));               ylabel('dHD(m)');
return;

