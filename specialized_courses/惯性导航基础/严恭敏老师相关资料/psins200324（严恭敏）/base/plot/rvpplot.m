function rvpplot(rvp)
% Standard deviation of Kalman vel&pos measurement noise plot.
%
% Prototype: rvpplot(rvp)
% Input: rvp - =[Var(vel], Var(pos), t] 
%
% See also  kfplot, inserrplot.

% Copyright(c) 2009-2020, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 02/03/2020
global glv
    myfigure
    subplot(211), plot(rvp(:,end), sqrt(rvp(:,1:3))); xygo('std(vel) / m/s')
    subplot(212), plot(rvp(:,end), sqrt([rvp(:,4:5)*glv.Re^2,rvp(:,6)])); xygo('std(pos) / m')
    