syms csf csr cz ch c d f ksf ksr ktf ktr kz kh mwf mwr mb mz mh Ib p1 p2 p3 p4
mwf=40;mwr=45.5;mb=690;mz=70;ktf=195;ktr=195;ksf=24;ksr=24;Ib=1222;mh=11.1;kz=8;csf=1.5;csr=1.5;cz=0.15;c=1.15;d=1.35;f=0.55;kh=4;ch=0.15;p1=0.078;p2=0.31*10^-4;p3=p2;p4=1;
ch=0.15;kh=4;
%% calculation of model coefficient
C=[csf,0,-csf,0,0,c*csf;0,csr,-csr,0,0,-d*csr;-csf,-csr,csf+csr+cz,-cz,0,-c*csf+d*csr+f*cz;0,0,-cz,cz+ch,-ch,-f*cz;0,0,0,-ch,ch,0;c*csf,-d*csr,-c*csf+d*csr+f*cz,-f*cz,0,c.^2*csf+d.^2*csr+f.^2*cz]
K=[ktf+ksf,0,-ksf,0,0,c*ksf;0,ktr+ksr,-ksr,0,0,-d*ksr;-ksf,-ksr,ksf+ksr+kz,kz,0,-c*ksf+d*ksr+f*kz;0,0,kz,kz+kh,-kh,-f*kz;0,0,0,-kh,kh,0;c*ksf,-d*ksr,-c*ksf+d*ksr+f*kz,-f*kz,0,c.^2*ksf+d.^2*ksr+f.^2*kz]
M=[mwf,0,0,0,0,0;0,mwr,0,0,0,0;0,0,mb,0,0,0;0,0,0,mz,0,0;0,0,0,0,mh,0;0,0,0,0,0,Ib;]
Kt=[ktf,0;0,ktr;0,0;0,0;0,0;0,0;]
KT=inv(M)*K
CT=inv(M)*C
QT=inv(M)*Kt
%% STATESPCACE model
A=[-CT,-KT;eye(6),zeros(6)];
B=[zeros(6,2);QT];
%% H=eye(12);
H=[0 0 0 0 0 0 0 0 0 -1 1 0;0 0 0 0 0 0 -1 0 0 0 1 0;0 0 0 0 0 0 0 -1 0 0 1 0;0 0 0 ch/mh -ch/mh 0 0 0 0 kh/mh -kh/mh 0];
D=zeros(12,2);
%% LQG control coefficient
Q=[zeros(3,12);0 0 0 p4*ch^2/mh^2 -p4*ch^2/mh^2 0 0 0 0 p4*kh*ch/mh^2 -p4*kh*ch/mh^2 0;0 0 0 -p4*ch^2/mh^2 p4*ch^2/mh^2 0 0 0 0 -p4*kh*ch/mh^2 p4*kh*ch/mh^2 0;zeros(1,12);0 0 0 0 0 0 p2 0 0 0 -p2 0;0 0 0 0 0 0 0 p3 0 0 -p3 0;zeros(1,12);0 0 0 p4*kh*ch/mh^2 -p4*kh*ch/mh^2 0 0 0 0 p4*kh^2/mh^2+p1 -(p4*kh^2/mh^2+p1) 0;0 0 0 -p4*kh*ch/mh^2 p4*kh*ch/mh^2 0 -p2 -p3 0 -(p4*kh^2/mh^2+p1) (p4*kh^2/mh^2+p1+p2+p3) 0;zeros(1,12)]
R=eye(2);
%% Optimization Gain
J=lqr(A,B,Q,R);
