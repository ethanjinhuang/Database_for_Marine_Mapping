clear;close all;clc;
Re = 6378160;   e = 1/298.3;   wie = 7.2921151467e-5;   g0 = 9.7803267714;
Tm  = 0.02; 
alpha = 36.001542377777781;
slti = sind(alpha); clti = cosd(alpha); tlti = slti/clti; slti2 = slti^2; slti4 = slti2^2;
RM = Re*(1-2*e+3*e*slti2)
RN = Re*(1+e*slti2)