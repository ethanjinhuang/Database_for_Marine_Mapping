function g = getg(pos)
%求重力加速度 g
global g0 e

    sl=sin(pos(1));  cl=cos(pos(1));  tl=sl/cl; sl2=sl^2; sl4=sl2^2;
    g = g0*(1+5.27094e-3*sl2+2.32718e-5*sl4) - 3.086e-6*pos(3);

