function [qnb, vnm, posm] = sins(qnbm_1, vnm_1, posm_1, dwvm1, dwvm2, Tm)
%捷联解算
global Re e g0 wie 
        wm1 = dwvm1(1:3);  wm2 = dwvm2(1:3);        vm1 = dwvm1(4:6);  vm2 = dwvm2(4:6);
        wm = wm1 + wm2;  vm = vm1 + vm2;
        %通用变量计算
	    slti = sin(posm_1(1)); clti = cos(posm_1(1)); tlti = slti/clti; slti2 = slti^2; slti4 = slti2^2;
        RM = Re*(1-2*e+3*e*slti2); RN = Re*(1+e*slti2); RMh = RM + posm_1(3); RNh = RN + posm_1(3);
        wnie = wie * [0; clti; slti];    wnen = [-vnm_1(2)/RMh; vnm_1(1)/RNh; vnm_1(1)/RNh*tlti];
        Cnb = Quat2Mat(qnbm_1);
        %姿态更新
        fim = wm + 2/3*cross(wm1,wm2);
	    wnin = wnie + wnen;
        jm = wnin*Tm;  ym = fim - Cnb'*jm; 
        qbm_1bm = Rv2Quat(ym);
        qnb = QuatMul(qnbm_1, qbm_1bm);
        %速度更新
        g = g0*(1+5.27094e-3*slti2+2.32718e-5*slti4) - 3.086e-6*posm_1(3);
  	    gn = [0; 0; -g];
        dvG_Corm = (gn - cross(2*wnie+wnen,vnm_1))*Tm;
        dvm = vm; dvrotm = 1/2*cross(wm, vm); dvsculm = 2/3*(cross(wm1,vm2)+cross(vm1,wm2));
        Cnn = Rv2Mat(1.0/2*jm); %[1, 0, 0; 0, 1, 0; 0, 0, 1] - 1.0/2*Asym(jm);
        %qnn = Rv2Quat(-jm);        Cnn = Quat2Mat(qnn);
        vnm = vnm_1 + Cnn*Cnb*(vm+dvrotm+dvsculm) + dvG_Corm;   
        %位置更新
        vnm_12 = 1/2*(vnm_1+vnm);
        posm = posm_1 + Tm*[vnm_12(2)/RMh; vnm_12(1)/(RNh*clti); vnm_1(3)];
