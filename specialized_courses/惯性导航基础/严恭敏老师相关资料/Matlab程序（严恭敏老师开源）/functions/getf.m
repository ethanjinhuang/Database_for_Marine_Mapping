function Ft = getf(Cnb, vn, pos, wb, fb, vnD, posD, tao)
%«ÛœµÕ≥æÿ’ÛFt
global Re e g0 wie
	sl = sin(pos(1)); cl = cos(pos(1)); tl = sl/cl; secl = 1/cl; secl2 = secl^2; sl2 = sl^2;
    RM = Re*(1-2*e+3*e*sl2); RN = Re*(1+e*sl2); 
    f_RMh = 1/(RM + pos(3)); f_RNh = 1/(RN + pos(3)); f_RMh2 = f_RMh^2; f_RNh2 = f_RNh^2;
    wnie = wie * [0; cl; sl];    wnen = [-vn(2)*f_RMh; vn(1)*f_RNh; vn(1)*f_RNh*tl];
    wnin = wnie + wnen;    
    %%%
    M1 = [0, 0, 0; -wie*sl, 0, 0; wie*cl, 0, 0];
    M2 = [0, -f_RMh, 0; f_RNh, 0, 0; f_RNh*tl, 0, 0];
    M3 = [0, 0, vn(2)*f_RMh2; 0, 0, -vn(1)*f_RNh2; vn(1)*secl2*f_RNh, 0, -vn(1)*tl*f_RNh2];
    M13 = M1+M3;
    M4 = Asym(vn)*M2 - Asym(2*wnie+wnen);                                                                         
    M5 = Asym(vn)*(2*M1+M3);
    M6 = [0, f_RMh, 0; secl*f_RNh, 0, 0; 0, 0, 1];
    M7 = [0, 0, -vn(2)*f_RMh2; vn(1)*secl*tl*f_RNh, 0, vn(1)*secl*f_RNh2; 0, 0, 0];
    aG = diag(-1./tao(1:3));
    S1 = Asym(wnin); S2 = Cnb*diag(wb); S3 = Asym(Cnb*fb); S4 = Cnb*diag(fb); 
    %%%
    slD = sin(posD(1)); clD = cos(posD(1)); tlD = slD/clD; seclD = 1/clD; slD2 = slD^2;
    RMD = Re*(1-2*e+3*e*slD2); RND = Re*(1+e*slD2); 
    f_RMhD = 1/(RM+posD(3)); f_RNhD = 1/(RN+posD(3)); f_RMhD2 = f_RMhD^2; f_RNhD2 = f_RNhD^2;
    MD6 = [0, f_RMhD, 0; seclD*f_RNhD, 0, 0; 0, 0, 1];
    MD7 = [0, 0, -vnD(2)*f_RMhD2; vnD(1)*seclD*tlD*f_RNhD, 0, vnD(1)*seclD*f_RNhD2; 0, 0, 0];
    SD1 = MD6*Asym(vnD); SD2 = MD6*vnD; aD = -1/tao(4);
    %%%
    aSv = diag(-1./tao(5:7)); aSp = diag(-1./tao(8:10));
        
    o13 = [0, 0, 0]; o31 = [0; 0; 0]; o3 = [0, 0, 0; 0, 0, 0; 0, 0, 0];   

    %%%%%  fi    dvn   dpos   dKG   eb    dKA   db  22  dposD dKD  26  dvnS  dposS
    Ft = [ -S1   M2    M13    -S2   -Cnb  o3    o3      o3    o31      o3    o3;
           S3    M4    M5     o3    o3    S4    Cnb     o3    o31      o3    o3;
           o3    M6    M7     o3    o3    o3    o3      o3    o31      o3    o3;
           o3    o3    o3     o3    o3    o3    o3      o3    o31      o3    o3;
           o3    o3    o3     o3    aG    o3    o3      o3    o31      o3    o3;
           o3    o3    o3     o3    o3    o3    o3      o3    o31      o3    o3;
           o3    o3    o3     o3    o3    o3    o3      o3    o31      o3    o3; 
       
           SD1   o3    o3     o3    o3    o3    o3      MD7   SD2      o3    o3; 
           o13   o13   o13    o13   o13   o13   o13     o13   aD       o13   o13;
           
           o3    o3    o3     o3    o3    o3    o3      o3    o31      aSv   o3;
           o3    o3    o3     o3    o3    o3    o3      o3    o31      o3    aSp ]; 
