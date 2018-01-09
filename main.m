clear ;
close all;
clc;
%% Parameters
h = 5000; % m (cruising altitude)

AR = 10;
L_D = 17.7; %L/D Guess
m_payload=1000; %Set payload weight
Range = 1000;% range in km

V_max_kmh = 500; %  maximum velocity not cruise (km/h)
V_cruise_kmh = 450;
V_stall_kmh = 145;
V_stall_mps = V_stall_kmh/3.6;
C_L_max = 1.8; %With flaps used during landing
Taper_rat = 0.4;  %Main wing taper ratio
e0=0.98;
per_a = 0.5; %Percentage of aileron span respect to total wing span 
g=9.81; 

F_rat = 8; %Fineness Ratio
T_vt = 0.7; %Vertical tail taper ratio
AR_vt = 0.6; %Aspect Ratio vertical tail
T_ht = 0.3; %Horizontal taper ratio
AR_ht = 2.8; %Aspect Ratio horizontal tail


%% Calculation
[T_sl, rho_sl, mu_sl] = atmosphere(0); % Sealevel
[T_cruise, rho_cruise, mu_cruise] = atmosphere(h); % Cruise altitude

[M_to,f_b] = Mass_Iteration1(Range,m_payload,L_D);
P_W = Power_Weight(V_max_kmh,L_D);
[W_S,S] = Wing_Loading(V_stall_kmh,V_cruise_kmh,C_L_max,AR,M_to,h,e0);

V_cruise_mps = V_cruise_kmh/3.6;

C_L_clean = W_S*g * 1/(0.5*V_stall_mps^2 *rho_sl);

[M_to_refined,M_e_refined,M_bat] = Mass_Iteration2(V_max_kmh,P_W,W_S,AR,f_b,m_payload);

[L_fuse, W_fuse,b,C_root,MAC,S_vt,b_vt,Cr_vt,S_ht,b_ht,Cr_ht,C_a,b_a,C_e,b_e,C_r,b_r,L_vt,L_ht]= Sizing(M_to_refined,F_rat,Taper_rat,AR,S,T_vt,AR_vt,T_ht,AR_ht,per_a);
[W_S,S] = Wing_Loading(V_stall_kmh,V_cruise_kmh,C_L_max,AR,M_to_refined,h,e0);
C_L_cruise = M_to_refined*g/(0.5*V_cruise_mps^2*rho_cruise*S);
C_L_landing = M_to_refined*g/(0.5*(V_stall_kmh/3.6)^2*rho_sl*S);

%Semi-Empirical case for now
S_fuse = pi*W_fuse^2*0.25 * 2 + L_fuse*pi*W_fuse; %Surface area of fuselage

%% costs
figure();
Q = linspace(300,1000);
RDTEflywaway = costs(Q, M_to_refined, M_bat, V_max_kmh);
plot(Q, RDTEflywaway./Q/1e6)
title('Plane Cost')
xlabel('Aircraft produced')
ylabel('Cost per plane [Million USD]')


%% airfoil selection
% [nameacc, Sacc, Pracc] = bestAirfoil(AR, S, M_to_refined*g, V_cruise_mps, e0, h, Cdf);
% 
%  figure()
%  scatter(Sacc, Pracc);
% xlabel('Wing area [m^2]')
% ylabel('Power required [Watt]')
% 
% idx = Sacc < 30 & Sacc > 0;
% Ssmall = Sacc(idx);
% Prsmall = Pracc(idx);
% namesmall = nameacc(idx);
% [Prtop,Pridx] = sort(Prsmall,'ascend');
% Stop = Ssmall(Pridx);
% nametop = namesmall(Pridx);
% 
% figure()
% scatter(Stop(1:20), Prtop(1:20))
% text(Stop(1:20), Prtop(1:20), nametop(1:20))
% xlabel('Wing area [m^2]')
% ylabel('Power required [Watt]')

%airfoil = 'coord_seligFmt/naca652415.dat';
airfoil = 'coord_seligFmt/ua79sf18.dat';
%airfoil = 'coord_seligFmt/nlf414f.dat';

figure()
[Re, M, Cl] = nondimensionalize(AR, S, V_cruise_kmh/3.6, M_to_refined, h);
pol = xfoil(airfoil, 'alfa', 0:0.5:15, Re, M, 'ppar n 200', 'oper iter 500');



% Cd3d = pol.CD + (pol.CL.^2)/(pi*e0*AR)+Cdf;
% LD = pol.CL./Cd3d;
% plot(pol.alpha,LD)
% title('L/D vs alpha');
% xlabel('Alpha');
% ylabel('L/D');

%% wing loading
figure();
hold on;
for AR = 8:12
    [wingloading, Pr, min_aoa] = bestWingLoading(AR, 10:2:40, M_to_refined*g, V_cruise_mps, e0, h, airfoil,S_fuse);
end
title('Power required at cruise')
xlabel('Wing area [m^2]')
ylabel('Power required [Watt]')
legend('show')
grid on;
print -depsc powerrequired

%% Output
f = figure('position',[500 250 400 600]);
t = uitable(f);
data = { 'MTOW (kg)',sprintf('%.f',M_to_refined);
        'Wing Area m^2',sprintf('%.2f',S) ;
        'W/S (kg/m^2)',sprintf('%.2f',W_S);
        'P/W (kwh/kg)',sprintf('%.3f',P_W);
        'Cl cruise',sprintf('%.3f',C_L_cruise);
        'Fuselage Length (m)',sprintf('%.3f',L_fuse);
        'Fuselage Width (m)', sprintf('%.3f',W_fuse);
        'Wing Span (m)', sprintf('%.3f',b); 
        'C root (m)', sprintf('%.3f',C_root);
        'VT area (m^2)', sprintf('%.3f',S_vt);
        'VT span (m)',sprintf('%.3f',b_vt) ;
        'VT chord (m)', sprintf('%.3f',Cr_vt);
        'L VT placement',L_vt;
        'HT area (m^2)',sprintf('%.3f',S_ht);
        'HT span (m)',sprintf('%.3f',b_ht);
        'HT chord (m)', sprintf('%.3f',Cr_ht);
        'Aileron span (m)', sprintf('%.3f',b_a);
        'Aileron Chord (m)', sprintf('%.3f',C_a);
        'Elevator span (m)',sprintf('%.3f',b_e);
        'Elevator chord (m)',sprintf('%.3f',C_e);
        'Rudder span (m)', sprintf('%.3f',b_r);
        'Rudder chord (m)', sprintf('%.3f',C_r)
        };
t.Data = data;
t.ColumnWidth = {200,'auto'};
t.Position = [25 25 380 600];
t.FontSize = 11; 
%PARAMETERS
f2 = figure('position',[500 250 400 600]);
t2 = uitable(f2);
data = {'AR',AR;
        'L/D',L_D;
        'Range (km)',1000;
        'Max Velocity km/h',V_max_kmh;
        'Cruise Velocity km/h',V_cruise_kmh;
        'Stall Velocity km/h',V_stall_kmh;
        'CL max (flaps)',C_L_max;
        'Taper Ratio',Taper_rat
        };
t2.Data = data;
t2.ColumnWidth = {200,'auto'};
t2.Position = [25 25 380 600];
t2.FontSize = 11; 


% For Turboprop, Historical values are :
% P/W = 0.33 (kW/kg)
% W/S = 195 (kg/m^2)
% AR = 9.2

%For reference Beech C90GTx
% Aspect Ratio = 9.7
% Wing Area S = 27.41
% MTOW = 4756 kg
% V cruise max = 504 km/h
% V stall = 140 km/h
% TO distance = 605 m
% Landing distance = 640 m 





