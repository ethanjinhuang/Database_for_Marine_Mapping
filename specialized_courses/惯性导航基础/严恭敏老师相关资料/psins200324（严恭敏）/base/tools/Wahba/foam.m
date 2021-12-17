function [A, k, lambda, A1] = foam(B, iter)
% Fast Optimal Attitude Matrix
% [A, signB, k, lambda, A1] = foam(B)
%
% See also  tr3, det3, inv3, adj3, svd3, foam, B33M44, svdest, vortech

% Copyright(c) 2009-2020, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/01/2020
    if nargin<2, iter=100; end
    detB = det(B);
    adjBp = adj3(B');  adjBp2 = tr3(adjBp*adjBp');
    BBp = B*B';  B2 = tr3(BBp);
    lambda = sqrt(3*B2);
    for k = 1:iter
        kappa = (lambda^2-B2)/2;
        zeta = kappa*lambda - detB;
		Psi = (lambda^2-B2)^2 - 8*lambda*detB - 4*adjBp2;
		dPsi = 8*zeta;
		dlambda = Psi / dPsi;
        lambda = lambda - dlambda;
		if dlambda<1e-15, break; end
    end
    A = ((kappa+B2)*B+lambda*adjBp-BBp*B)/zeta;
    if nargout==4
        [u, s, v] = svd(B);   A1 = u*v';
    end
