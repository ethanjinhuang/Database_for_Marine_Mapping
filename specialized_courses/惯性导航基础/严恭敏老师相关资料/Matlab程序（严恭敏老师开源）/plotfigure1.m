function fff
global Re e wie g0 ppm ug deg min sec hur dph   %全局变量

fid = fopen('e:/ygm/vehicle/kf4.bin','r');
%fid = fopen('e:/ygm/vehicle/kf3.bin','r');
for k=1:1:980
    [data, n] = fread(fid, 31+12+25+27, 'real*8');   %fwrite(fout, [xk;zk;Xk1;Xk2], 'real*8');
    xk(k,:) = data(1:31)';
    zk(k,:) = data(32:43)';   
    xk1(k,:) = data(44:68)';
    xk2(k,:) = data(69:95)';
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
subplot(4,4,9);  plot(Re*xk(:,22:23));      ylabel('dLtiD dLgiD(m)');
subplot(4,4,10); plot(xk(:,24));               ylabel('dHD(m)');
subplot(4,4,11); plot(xk(:,25));               ylabel('dKD');            
subplot(4,4,13); plot(xk(:,26:28));            ylabel('dVnxS dVnyS dVnzS(m/s)');
subplot(4,4,14); plot(Re*xk(:,29:30));      ylabel('dLtiS dLgiS(m)');
subplot(4,4,15); plot(xk(:,31));               ylabel('dHS(m)');           
figure;%观测值 
subplot(2,3,1);  plot(zk(:,1:3));              ylabel('z11 z12 z13');       
subplot(2,3,2);  plot(1/min*zk(:,4:5));        ylabel('z14 z15');
subplot(2,3,3);  plot(zk(:,6));                ylabel('z16');
subplot(2,3,4);  plot(zk(:,7:9));              ylabel('z21 z22 z23');       
subplot(2,3,5);  plot(1/min*zk(:,10:11));      ylabel('z24 z25');
subplot(2,3,6);  plot(zk(:,12));               ylabel('z16');
figure;%滤波值1
xk1(:,1:21) = xk1(:,1:21)-xk(:,1:21);
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
figure;%滤波值2
xk2(:,1:21) = xk2(:,1:21)-xk(:,1:21);
subplot(4,4,1);  plot(1/min*xk2(:,1:3));        ylabel('fx fy fz(min)');       
subplot(4,4,2);  plot(xk2(:,4:6));              ylabel('dVnx dVny dVnz(m/s)');
subplot(4,4,3);  plot(Re*xk2(:,7:8));        ylabel('dLti dLgi(min)');      
subplot(4,4,4);  plot(xk2(:,9));                ylabel('dH(m)');
subplot(4,4,5);  plot(1/ppm*xk2(:,10:12));      ylabel('dKGx dKGy dKGz(ppm)'); 
subplot(4,4,6);  plot(1/dph*xk2(:,13:15));      ylabel('ebx eby ebz(deg/h)');    
subplot(4,4,7);  plot(1/ppm*xk2(:,16:18));      ylabel('dKAx dKAy dKAz(ppm)');
subplot(4,4,8);  plot(1/ug*xk2(:,19:21));       ylabel('dbx dby dbz(ug)');    
subplot(4,4,13); plot(xk2(:,22:24));            ylabel('dVnxS dVnyS dVnzS(m/s)');
subplot(4,4,14); plot(Re*xk2(:,25:26));      ylabel('dLtiS dLgiS(min)');
subplot(4,4,15); plot(xk2(:,27));               ylabel('dHS(m)');   