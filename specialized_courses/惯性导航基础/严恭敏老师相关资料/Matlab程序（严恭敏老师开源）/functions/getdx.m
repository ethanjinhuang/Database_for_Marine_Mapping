function dx = getdx(x, Cnb, vn, pos, wb, fb, vnD, posD, tao)
global Re e wie  

%%%%%  fi    dvn   dpos   dKG   eb    dKA   db  22  dposD dKD  26  dvnS  dposS
    fi=x(1:3);  dvn=x(4:6); dpos=x(7:9);  dKG=x(10:12);  eb=x(13:15);  dKA=x(16:18); db=x(19:21);
    dposD=x(22:24); dKD=x(25); dvnS=x(26:28); dposS=x(29:31);
    
    %%
    sl = sin(pos(1)); cl = cos(pos(1)); tl = sl/cl; secl = 1/cl; secl2 = secl^2; sl2 = sl^2;
    RM = Re*(1-2*e+3*e*sl2); RN = Re*(1+e*sl2); 
    f_RMh = 1/(RM + pos(3)); f_RNh = 1/(RN + pos(3)); f_RMh2 = f_RMh^2; f_RNh2 = f_RNh^2;
    wnie = wie * [0; cl; sl];    wnen = [-vn(2)*f_RMh; vn(1)*f_RNh; vn(1)*f_RNh*tl];
    wnin = wnie + wnen;    

    dwnib = Cnb*(dKG.*wb+eb);
    dfn   = Cnb*(dKA.*fb+db);
    dwnie = [0; -wie*sl*dpos(1); wie*cl*dpos(1)];
    dwnen = [-f_RMh*dvn(2)+f_RMh2*vn(2)*dpos(3); f_RNh*dvn(1)-f_RNh2*vn(1)*dpos(3); f_RNh*tl*dvn(1)-f_RNh2*vn(1)*tl*dpos(3)+f_RNh*vn(1)*secl2*dpos(1)];
    dwnin = dwnie+dwnen;
    
    dfi = -cross(wnin,fi)+dwnin-dwnib; 
    ddvn = cross(Cnb*fb,fi)-cross(2*wnie+wnen,dvn)+cross(vn,2*dwnie+dwnen)+dfn;
    ddLti = f_RMh*dvn(2)-f_RMh2*vn(2)*dpos(3);
    ddLgi = f_RNh*secl*dvn(1)+f_RNh*vn(1)*secl*tl*dpos(1)+f_RNh2*vn(1)*secl*dpos(3);
    ddHgt = dvn(3);
    ddKG = [0; 0; 0];
    deb = -1./tao(1:3).*eb;% + w(1:3);
    ddKA = [0; 0; 0];
    ddb = [0; 0; 0];
    
    %%
    slD = sin(posD(1)); clD = cos(posD(1)); tlD = slD/clD; seclD = 1/clD; slD2 = slD^2;
    RMD = Re*(1-2*e+3*e*slD2); RND = Re*(1+e*slD2); 
    f_RMhD = 1/(RM+posD(3)); f_RNhD = 1/(RN+posD(3)); f_RMhD2 = f_RMhD^2; f_RNhD2 = f_RNhD^2;
 
    dvnD = cross(vnD,fi)+vnD*dKD;
    
    ddLtiD = f_RMhD*dvnD(2)-f_RMhD2*vnD(2)*dposD(3);
    ddLgiD = f_RNhD*seclD*dvnD(1)+f_RNhD*vnD(1)*seclD*tlD*dposD(1)+f_RNhD2*vnD(1)*seclD*dposD(3);
    ddHgtD = dvnD(3);
    ddKD = -1/tao(4)*dKD;% + w(4);
    
    %%
    ddvnS = -1./tao(5:7).*dvnS;% + w(5:10);
    ddposS = -1./tao(8:10).*dposS;% + w(5:10);
    
    dx = [dfi; ddvn; ddLti;ddLgi;ddHgt; ddKG; deb; ddKA; ddb;    ddLtiD;ddLgiD;ddHgtD; ddKD;   ddvnS; ddposS] ;
    
    
