function trajectory
%生成轨迹
global Re e wie g0 ppm ug deg min sec hur dph PKG PKA  %全局变量

    glb;
    th = 0.002;
    
    att = pi/180*[0; 0; 0];	 %[pitch roll azimuth]  初始姿态角
	vb = [0; 0; 0];          %[0 vby 0]       初始速度
	pos = [34*deg+14.76289014*min; 108*deg+54.57983*min; 380];  %[latitude longitude height]初始位置

    %tr(1:3) att(pitch roll azimuth);   姿态角
    %tr(2:6) vn(vnE vnN vnU);           速度
    %tr(7:9) pos(lti lgi hgt)           位置
    %tr(10:12) wm;                      角增量
    %tr(13:15) vm;                      比力增量    
    tr = [att; Att2Mat(att)*vb; pos; zeros(6,1)];
    
    fid = fopen('e:/ygm/vehicle/trace3__.bin','wb');  %4 good
    wat2 = getwat(0);
    for k=2:2:5400*1000
        wat0 = wat2;        k1 = getdtr(tr,         wat0);
        wat1 = getwat(k-1); k2 = getdtr(tr+th/2*k1, wat1);
                            k3 = getdtr(tr+th/2*k2, wat1);
        wat2 = getwat(k);   k4 = getdtr(tr+  th*k3, wat2);
        tr = tr + th/6*(k1+2*k2+2*k3+k4);
        if mod(k,10)==0                               %10 ms
            fwrite(fid, tr, 'real*8');   %保存数据
            if mod(k,1000)==0
                step=k/1000,                          %进度显示
            end
        end
    end
    fclose(fid);
        
function wat = getwat(k)          %轨迹设置
%姿态角变化率 wat(1:3) [wPitch wRoll wAzimuth]; 轨迹加速度 wat(4:6) [atx aty atz]
    k=mod(k,1830*1000);
%阶段   动作      起始时间(s)     航向azimuth(deg)  俯仰pitch(deg)  阶段初速度(m/s)
%                         持续时间(s)       倾斜roll(deg)   加速度(m/s^2)   
%       停止       0      30     0          0       0       0          0
        if k<=30*1000
            wat=[0; 0; 0;    0; 0; 0];
%       加速       30     10      0         0       0       1          0
        elseif k<=40*1000
            wat=[0; 0; 0;     0; 1; 0];
%       匀速       40     460      0         0       0       0          10
        elseif k<=500*1000
            wat=[0; 0; 0;    0; 0; 0];
%       右转弯     500     10       0/9      0       0        0          10
        elseif k<=510*1000
            v=10; w=9*pi/180; a=v*w;
            wat=[0; 0; w;     a; 0; 0]; 
%       匀速       510     490      90        0       0         0         10
        elseif k<=1000*1000
            wat=[0; 0; 0;    0; 0; 0];
%       左转弯     1000     10    90/-9       0        0        0          10
        elseif k<=1010*1000
            v=10; w=9*pi/180; a=v*w;
            wat=[0; 0 ; -w;     -a; 0; 0]; 
%       匀速       1010     490      0        0       0         0         10
        elseif k<=1500*1000
            wat=[0; 0; 0;    0; 0; 0];
%       进入爬升   1500     10       0        0       0/1.5       0         10
        elseif k<=1510*1000
            v=10; w=1.5*pi/180; a=v*w;
            wat=[w; 0; 0;   0; 0; a];
%       匀速       1510    100      0          0       15        0         10
        elseif k<=1610*1000
            wat=[0; 0; 0;    0; 0; 0];
%       退出爬升   1610    10       0           0       15/-1.5     0         10
        elseif k<=1620*1000
            v=10; w=1.5*pi/180; a=v*w;
            wat=[-w; 0; 0;   0; 0; -a];
%       匀速       1620    200      0          0        0        0          10
        elseif k<=1820*1000
            wat=[0; 0; 0;    0; 0; 0];
%       减速       1820     10      0          0        0        -1         10
        elseif k<=1830*1000
            wat=[0; 0; 0;    0; -1; 0];
%       停止       1830     ...      0          0        0        0         0
        else
            wat=[0; 0; 0;    0; 0; 0];
        end
        
