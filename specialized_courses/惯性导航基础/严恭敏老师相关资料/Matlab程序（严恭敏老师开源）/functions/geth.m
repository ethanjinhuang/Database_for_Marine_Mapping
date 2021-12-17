function Ht = geth(vnD)
%«Ûπ€≤‚æÿ’ÛHt
    o3 = [0,0,0; 0,0,0; 0,0,0]; o31 = [0; 0; 0]; I3=[1,0,0; 0,1,0; 0,0,1];
    SD1 = Asym(vnD);
    %%%%%  fi    dvn   dpos   dKG   eb    dKA   db %22  dposD dKD %26  dvnS  dposS
    Ht = [ -SD1  I3    o3     o3    o3    o3    o3      o3    -vnD     o3    o3;
           o3    o3    I3     o3    o3    o3    o3      -I3   o31      o3    o3;
           o3    I3    o3     o3    o3    o3    o3      o3    o31      -I3   o3;
           o3    o3    I3     o3    o3    o3    o3      o3    o31      o3    -I3 ]; 
