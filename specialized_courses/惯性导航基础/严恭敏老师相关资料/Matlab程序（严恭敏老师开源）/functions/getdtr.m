function dtr = getdtr(tr, wat)        
%tr(1:3) att(pitch roll azimuth);   tr(2:6) vn(vnE vnN vnU);   tr(7:9) pos(lti lgi hgt)
%tr(10:12) wm;  tr(13:15) vm;
%姿态角变化率 wat(1:3,i) [wPitch wRoll wAzimuth]; 轨迹加速度 wat(4:6,i) [atx aty atz]
global Re e g0 wie 
    %通用变量计算
	si = sin(tr(1)); ci = cos(tr(1)); sj = sin(tr(2)); cj = cos(tr(2)); sk = sin(tr(3)); ck = cos(tr(3));
	slti = sin(tr(7)); clti = cos(tr(7)); tlti = slti/clti; slti2 = slti^2; slti4 = slti2^2;
    RM = Re*(1-2*e+3*e*slti2); RN = Re*(1+e*slti2); RMh = RM + tr(9); RNh = RN + tr(9);
    wnie = wie * [0; clti; slti];	wnen = [-tr(5)/RMh; tr(4)/RNh; tr(4)/RNh*tlti];
    Cnb = [ cj*ck+si*sj*sk, ci*sk,  sj*ck-si*cj*sk;
           -cj*sk+si*sj*ck, ci*ck, -sj*sk-si*cj*ck;
           -ci*sj,          si,     ci*cj ];
    wt = wat(1:3);  at = wat(4:6);
	%姿态增量
	datt = wt;
	%速度增量
	Cnt = [ ck, ci*sk, -si*sk; 
           -sk, ci*ck, -si*ck; 
            0,  si,     ci ];
	dvn = Cnt*at;
	%位置增量
    dpos = [tr(5)/RMh; tr(4)/(RNh*clti); tr(6)];
	%角增量
    wnin = wnie + wnen;
	a2w = [ cj, 0,  sj*ci; 
            0,  1, -si; 
            sj, 0, -cj*ci ];
    wbnb = a2w*wt;
	dwm = Cnb'*wnin + wbnb;
	%比力增量
    g = g0*(1+5.27094e-3*slti2+2.32718e-5*slti4) - 3.086e-6*tr(9);
	gn = [0; 0; -g];
  	dvm = Cnb' * (dvn+cross(2*wnie+wnen,tr(4:6))-gn);
    %rk增量
    dtr = [datt; dvn; dpos; dwm; dvm];