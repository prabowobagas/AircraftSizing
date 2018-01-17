% total drag wings
clear,clc

%% General numbers
[T, rho, vis] = atmosphere(5000);
%rho=0.7364;         % density kg/m3
V=139;              % velocity m/s
%vis=1.628*10^-5;    % dynamic viscosity Pas
c=1.7;
L=8;
S = 2*L*c;
Cl=0.3041;
e=0.8;
AR=9.7;
xtr = 0.57; % from xfoil

%% Wing Drag
Re=(rho*V*c)/vis;  % all TURBULENT SAME AS CF_total
Recr = Re*xtr;
Cfxl = 1.328/sqrt(Recr);
Cfct = 0.074/Re^0.2;
Cfxt = 0.074/Recr^0.2;
Cf_plate=xtr*Cfxl+Cfct-xtr*Cfxt; 

Cf_wing = Cf_plate*2; % 2 sides, pressure

q=0.5*rho*V^2;  % Dynamic pressure
Cdi=Cl^2/(pi*e*AR);
CD=Cf_wing+Cdi;

Drag_wing=CD*q*S

powerrequired(AR, S, 6000*9.8, V, e, 5000, 'coord_seligFmt/naca652415.dat', 11.2, 1.4);
%% Fuselage Drag
C_fuse=11.2;

Re=(rho*139*C_fuse)/vis;  % all TURBULENT SAME AS CF_total
Cf_fuse=0.074/(Re^(1/5));

A=pi*1.4*C_fuse;

Drag_fuselage=q*A*Cf_fuse*2  % Total drag of fuselage in Nwe


%% Total Drag

Total_drag=Drag_wing+Drag_fuselage

LD=5600*9.81/Total_drag