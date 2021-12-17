function AlignmentVerVt
%验证初始对准中失准角、速度的变化规律
global Re e wie g0 ppm ug deg min sec hur dph   %全局变量

glb;

	pos = [34*deg+14.76289014*min; 108*deg+54.57983*min; 380];  %[latitude longitude height]
    sl = sin(pos(1));  cl = cos(pos(1));         g = getg(pos);

    fi=[1; 1; 0]*min;
    fiE0=fi(1);  fiN0=fi(2);  fiU0=fi(3);
    
    uE=fiN0*wie*sl-fiU0*wie*cl;
    uN=-fiE0*wie*sl;
    uU=fiE0*wie*cl;
    
    for k=1:1:100
        t=100*k;    t2=t^2;     t3=t^3;
        %正向解真值
        fiE=fiE0+uE*t+t2/2*wie*(uN*sl-uU*cl);               %失准角
        fiN=fiN0+uN*t-t2/2*uE*wie*sl;
        fiU=fiU0+uU*t+t2/2*uE*wie*sl;
        
        vE=-g*fiN0*t-t2/2*g*uN+t3/6*g*wie*uE*sl;           %东向、北向速度
        vN=g*fiE0*t+t2/2*g*uE+t3/6*g*wie*(uN*sl-uU*cl);

        xx(k,:)=[[fiE,fiN,fiU]./min,vE,vN];
    end
    
    figure;
    subplot(2,3,1), plot(xx(:,1));    subplot(2,3,2), plot(xx(:,2));        subplot(2,3,3), plot(xx(:,3)); 
    subplot(2,3,4), plot(xx(:,4));    subplot(2,3,5), plot(xx(:,5)); 
