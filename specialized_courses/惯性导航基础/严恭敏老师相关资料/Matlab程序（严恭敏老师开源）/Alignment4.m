function Alignment
%参数辨识精对准（方位失准角主要从北向速度的二次方向中提取!）
global Re e wie g0 ppm ug deg min sec hur dph   %全局变量

    glb;
    
Tm=2;

    att = pi/180*[0; 0; 0];	                   %[pitch roll azimuth]  
	vb = [0; 0; 0];          vn = Att2Mat(att)*vb;  %[0 vby 0]       
	pos = [34*deg+14.76289014*min; 108*deg+54.57983*min; 380];  %[latitude longitude height]
    
    tr = [att; vn; pos; zeros(6,1)];
    wat=zeros(6,1); 
    tr = getdtr(tr, wat);    
    dwvm=tr(10:15)*Tm/2;            % IMU 输出
    
    qnb = Att2Quat(att);  
    fi=[1; 1; 0]*min;  
    qnb = QuatMul(Rv2Quat(-fi), qnb);  %带误差四元数
   % pos = pos + dpos; 
    
    theE=zeros(3,1); theN=theE;
    PE=0.1*eye(3);    PN=PE;
    
    qnb0=qnb;                       %保存四元数初值
    
    for k=1 :1: 5000
        
 %       ns = randn(3,2);
  %      eb = e_tao(1:3).* eb  + sqw(1:3).*ns(:,1);
   %     dwvm1 = [KG*dwvm(1:3); KA*dwvm(4:6)] + [eb; db]*Tm/2;
    %    eb = e_tao(1:3).* eb  + sqw(1:3).*ns(:,2);
     %   dwvm2 = [KG*dwvm(1:3); KA*dwvm(4:6)] + [eb; db]*Tm/2;        %刻度系数、偏置、晃动
     dwvm1=dwvm; dwvm2=dwvm;
        [qnb,vn,posm]=sins1(qnb, vn, pos, dwvm1, dwvm2, Tm);  % more than --  vn=vn+Cnb*vm;
        
        H=[k*Tm, (k*Tm)^2, (k*Tm)^3];
        vE=vn(1); vN=vn(2);
        
        theE=theE+PE*H'*(vE-H*theE);        %东向
        tao=1/(1+H*PE*H');
        PE=PE-tao*PE*H'*H*PE;
        theN=theN+PN*H'*(vN-H*theN);        %北向
        tao=1/(1+H*PN*H');
        PN=PN-tao*PN*H'*H*PN;
        
        if mod(k,50)==0
            step=k/50,
            
            sl = sin(pos(1)); cl = cos(pos(1)); tl = sl/cl;            g = getg(pos);
            a1E=theE(1); a2E=theE(2); a3E=theE(3);           a1N=theN(1); a2N=theN(2); a3N=theN(3);
    
            fiN0=a1E/(-g);
            %uN=a2E/(-1/2*g);
            %uE1=a3E/(1/6*wie*sl);       
            fiE0=a1N/g;
            uE=a2N/(1/2*g);
            %uU=(a3N/(1/6*g*wie)-uN*sl)/(-cl);
            fiU0=fiN0*tl-uE/(wie*cl);
            
            fi1=[fiE0, fiN0, fiU0];
            xx(step,:)=[vn', fi1./min];
        end
    end
    figure;
    subplot(2,3,1), plot(xx(:,1));    subplot(2,3,2), plot(xx(:,2));    subplot(2,3,3), plot(xx(:,3));
    xx(1:10,:)=[];
    subplot(2,3,4), plot(xx(:,4));    subplot(2,3,5), plot(xx(:,5));    subplot(2,3,6), plot(xx(:,6));

    qnb=QuatMul(Rv2Quat(fi1'),qnb0);     %姿态四元数修正
    
    att1=Quat2Att(qnb),
    err=(att1*deg-att)/min,
    
function [qnb, vnm, posm] = sins1(qnbm_1, vnm_1, posm_1, dwvm1, dwvm2, Tm)
%捷联解算
global Re e g0 wie 
        wm1 = dwvm1(1:3);  wm2 = dwvm2(1:3);        vm1 = dwvm1(4:6);  vm2 = dwvm2(4:6);
        wm = wm1 + wm2;  vm = vm1 + vm2;
        %通用变量计算
	    slti = sin(posm_1(1)); clti = cos(posm_1(1)); %tlti = slti/clti; slti2 = slti^2; slti4 = slti2^2;
        %RM = Re*(1-2*e+3*e*slti2); RN = Re*(1+e*slti2); RMh = RM + posm_1(3); RNh = RN + posm_1(3);
        wnie = wie * [0; clti; slti];    
        %wnen = [-vnm_1(2)/RMh; vnm_1(1)/RNh; vnm_1(1)/RNh*tlti];
    wnen = [0; 0; 0];
        Cnb = Quat2Mat(qnbm_1);
        %姿态更新
        fim = wm + 2/3*cross(wm1,wm2);
	    wnin = wnie + wnen;
        qbm_1bm = Rv2Quat(fim - Cnb'*wnin*Tm);
        qnb = QuatMul(qnbm_1, qbm_1bm);
        %速度更新
        %g = g0*(1+5.27094e-3*slti2+2.32718e-5*slti4) - 3.086e-6*posm_1(3);
  	    %gn = [0; 0; -g];
        %dvG_Corm = (gn - cross(2*wnie+wnen,vnm_1))*Tm;
    dvG_Corm = [0; 0; 0];      
        dvm = vm; dvrotm = 1/2*cross(wm, vm); dvsculm = 2/3*(cross(wm1,vm2)+cross(vm1,wm2));
        Cnn = Rv2Mat(wnin*Tm/2); 
        vnm = vnm_1 + Cnn*Cnb*(dvm+dvrotm+dvsculm) + dvG_Corm;  
        %位置更新
        %vnm_12 = 1/2*(vnm_1+vnm);
        %posm = posm_1 + Tm*[vnm_12(2)/RMh; vnm_12(1)/(RNh*clti); vnm_1(3)];
    posm = [0; 0; 0];       

