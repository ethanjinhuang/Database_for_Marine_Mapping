function fff
global Re e wie g0 ppm ug deg min sec hur dph   %全局变量

fid = fopen('e:/ygm/vehicle/federal_rk3.bin','r');
%fid = fopen('e:/ygm/vehicle/kf3.bin','r');
for k=1:1:500
    [data, n] = fread(fid, 31+31, 'real*8');   %fwrite(fout, [xk;zk;Xk1;Xk2], 'real*8');
    x0(k,:) = data(1:31)';
    x1(k,:) = data(32:62)';
end
xk=x0;
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
xk=x1;
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
xk=xk-x0;
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

