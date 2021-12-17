fid = fopen('e:/ygm/vehicle/trace3_.bin','r'); 
glb;
%status = fseek(fid,21*8*400*100,0);
for k=1:1:5400
    [data, n] = fread(fid, 15, 'real*8'); %fwrite(fid, [rk; k4(10:15)], 'real*8');
    ninf(:,k) = data;
    [data, n] = fread(fid, 15*99, 'real*8');
end
figure
title('¹ì¼£²ÎÊý');
subplot(3,3,1); plot(ninf(1,:)/deg);    xlabel('t(s)'); ylabel('pitch(deg)');       
subplot(3,3,2); plot(ninf(2,:)/deg);    xlabel('t(s)'); ylabel('roll(deg)');   
subplot(3,3,3); plot(ninf(3,:)/deg);    xlabel('t(s)'); ylabel('yaw(deg)');   
subplot(3,3,4); plot(ninf(4,:));    xlabel('t(s)'); ylabel('VE(m/s)');   
subplot(3,3,5); plot(ninf(5,:));    xlabel('t(s)'); ylabel('VN(m/s)');   
subplot(3,3,6); plot(ninf(6,:));    xlabel('t(s)'); ylabel('VU(m/s)');   
subplot(3,3,7); plot(ninf(7,:)/deg);    xlabel('t(s)'); ylabel('latitude(deg)');   
subplot(3,3,8); plot(ninf(8,:)/deg);    xlabel('t(s)'); ylabel('longitude(deg)');   
subplot(3,3,9); plot(ninf(9,:));    xlabel('t(s)'); ylabel('height(deg)');   