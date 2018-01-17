clear,clc

[T, rho, vis] = atmosphere(5000);

V=125;              % velocity m/s
q=0.5*rho*V^2;
W=5648*9.81;        % Weight
%x1=0.71;            % distance wing to CG
x1=0.3;
x2=6.93;            % distance tail to cg

ct=1.47;            % cord tail
Clt0=0.113;        % cl van de tail
At=6.067;           % Area tail

cw=1.7;             % cord wing
Aw=30.95;           % Area wing
Clw=0.3114;         % cl wing
Cmw=0.082;          % moment coefficient wing
Lw=Clw*q*Aw;        % lift wing

aw = 2*pi/(1+57.3*2*pi/(pi*0.9*10));
at = 2*pi/(1+57.3*2*pi/(pi*0.9*3));

Vt = (x1+x2)*At/(Aw*cw)
xcg=4.20;
xn=0.1*cw; %+xcg from the nose
x1 = Vt*at/aw-xn;

Mwcm=Cmw*q*Aw*cw;      % moment by the cm of the wing  At Cg
Mwx1=-Lw*x1;            % moment by the lift of the wing  at cg

Lt=-(Mwx1+Mwcm)/x2
alpha = Lt/(Clt0*q*At)

dMda = aw*(x1/cw-Vt*at/aw)