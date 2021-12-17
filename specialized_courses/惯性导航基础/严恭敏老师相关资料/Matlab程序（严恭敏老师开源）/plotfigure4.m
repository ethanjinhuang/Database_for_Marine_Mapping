function plottype
global Re e wie g0 ppm ug deg min sec hur dph   %全局变量

fid = fopen('e:/ygm/vehicle/kf5.bin','r');
for k=1:1:4000
    [data, n] = fread(fid, 25+25, 'real*8');   %fwrite(fout, [xk;zk1;Xk1;Xk2], 'real*8');
    err(k,:) = data(1:25)';
    xk(k,:) = data(26:50)';
end
figure;
%xk1(:,1:25) = xk1(:,1:25) - xk(:,1:25);
subplot(4,4,1);  plot(1/min*err(:,1:3));         xlabel('t(s)');ylabel('fx fy fz(min)');       
subplot(4,4,2);  plot(err(:,4:6));               xlabel('t(s)');ylabel('dVnx dVny dVnz(m/s)');
subplot(4,4,3);  plot(Re*err(:,7:8));         xlabel('t(s)');ylabel('dLti dLgi(m)');      
subplot(4,4,4);  plot(err(:,9));                 xlabel('t(s)');ylabel('dH(m)');
subplot(4,4,5);  plot(1/ppm*err(:,10:12));       xlabel('t(s)');ylabel('dKGx dKGy dKGz(ppm)'); 
subplot(4,4,6);  plot(1/dph*err(:,13:15));       xlabel('t(s)');ylabel('ebx eby ebz(deg/h)');    
subplot(4,4,7);  plot(1/ppm*err(:,16:18));       xlabel('t(s)');ylabel('dKAx dKAy dKAz(ppm)');
subplot(4,4,8);  plot(1/ug*err(:,19:21));        xlabel('t(s)');ylabel('dbx dby dbz(ug)');
subplot(4,4,9);  plot(Re*err(:,22:23));       xlabel('t(s)');ylabel('dLtiD dLgiD(m)');
subplot(4,4,10); plot(err(:,24));                xlabel('t(s)');ylabel('dHD(m)');
subplot(4,4,11); plot(err(:,25));                xlabel('t(s)');ylabel('dKD');      
figure;
err=err-xk;
subplot(4,4,1);  plot(1/min*err(:,1:3));         xlabel('t(s)');ylabel('err-fx fy fz(min)');       
subplot(4,4,2);  plot(err(:,4:6));               xlabel('t(s)');ylabel('err-dVnx dVny dVnz(m/s)');
subplot(4,4,3);  plot(Re*err(:,7:8));         xlabel('t(s)');ylabel('err-dLti dLgi(m)');      
subplot(4,4,4);  plot(err(:,9));                 xlabel('t(s)');ylabel('err-dH(m)');
subplot(4,4,5);  plot(1/ppm*err(:,10:12));       xlabel('t(s)');ylabel('err-dKGx dKGy dKGz(ppm)'); 
subplot(4,4,6);  plot(1/dph*err(:,13:15));       xlabel('t(s)');ylabel('err-ebx eby ebz(deg/h)');    
subplot(4,4,7);  plot(1/ppm*err(:,16:18));      xlabel('t(s)');ylabel('err-dKAx dKAy dKAz(ppm)');
subplot(4,4,8);  plot(1/ug*err(:,19:21));        xlabel('t(s)');ylabel('err-dbx dby dbz(ug)');
subplot(4,4,9);  plot(Re*err(:,22:23));       xlabel('t(s)');ylabel('err-dLtiD dLgiD(m)');
subplot(4,4,10); plot(err(:,24));                xlabel('t(s)');ylabel('err-dHD(m)');
subplot(4,4,11); plot(err(:,25));                xlabel('t(s)');ylabel('err-dKD');      
