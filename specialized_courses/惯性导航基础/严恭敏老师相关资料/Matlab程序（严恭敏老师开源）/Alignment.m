function Alignment
%提供外界航向角基准的初始对准
global Re e wie g0 ppm ug deg min sec hur dph   %全局变量

glb;

KA = [100; -99; 98]*ppm;	 db = [-100; -99; 98]*ug;
%KA = [0; 0; 0]*ppm;	
%db = [0; 0; 0]*ug;
th = 0.01;

    att = pi/180*[3.; 10; 0];	 %[pitch roll azimuth]  
	vb = [0; 0; 0];          %[0 vby 0]       
	pos = [34*deg+14.76289014*min; 108*deg+54.57983*min; 380];  %[latitude longitude height]
    %rk(1:3) att(pitch roll azimuth);   rk(2:6) vn(vnE vnN vnU);   rk(7:9) pos(lti lgi hgt)
    %rk(10:12) wm;  rk(13:15) vm;         
    tr = [att; Att2Mat(att)*vb; pos; zeros(6,1)];
    wt=zeros(3,1); at=zeros(3,1);
    tr = getdtr(tr, wt, at);    
    dvm=tr(13:15)*th;
    
    vm=[0;0;0];
    
    Tm=19.3;
    for k=1:1:Tm/th  %20s
        vm = vm+([1;1;1]+KA).*dvm + db*th + th*0.001*g0*sin(6.36*1*k*th);        %刻度系数、偏置、晃动
    end
    
    gvAzi = att(3)+3*min;
%    gvAzi = att(3);
    lti=pos(1)+10/Re;     lgi=pos(2)+10/Re;     hgt=pos(3)+10;
%    lti=pos(1);     lgi=pos(2);     hgt=pos(3);
    slti = sin(lti); slti2 = slti^2; slti4 = slti2^2;
    g = g0*(1+5.27094e-3*slti2+2.32718e-5*slti4) - 3.086e-6*hgt;

    si=vm(2)/(g*Tm);    ci=sqrt(1-si*si);
    tj=-vm(1)/vm(3);     cj=1/sqrt(1+tj*tj);     sj=tj*cj;
    sk=sin(gvAzi);     ck=cos(gvAzi);
    Cnb = [ cj*ck+si*sj*sk, ci*sk,  sj*ck-si*cj*sk;
           -cj*sk+si*sj*ck, ci*ck, -sj*sk-si*cj*ck;
           -ci*sj,          si,     ci*cj ];
       
    att1 = Mat2Att(Cnb),
    err = (att1-att/deg)*60,

        