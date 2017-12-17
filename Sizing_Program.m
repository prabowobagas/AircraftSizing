clear all;
close all;
clc; 
%% Parameters for planes
AR = 9.7;
L_D = 17; %L/D Guess
m_payload=1000; %Set payload weight
Range = 1000;% range in km

V_max_kmh = 500; %  maximum velocity not cruise (km/h)
V_cruise_kmh = 450;
V_stall_kmh = 145; 
C_L_max = 2; %With flaps used
Taper_rat = 0.4;  %Main wing taper ratio
e0=0.8; %Oswald efficiency
per_a = 0.5; %Percentage of aileron span respect to total wing span 
g=9.81; 

F_rat = 6; %Fineness Ratio
T_vt = 0.7; %Vertical tail taper ratio
AR_vt = 0.6; %Aspect Ratio vertical tail
T_ht = 0.3; %Horizontal taper ratio
AR_ht = 2.8; %Aspect Ratio horizontal tail

%% Calculation
[M_to,f_b] = Mass_Iteration1(Range,m_payload,L_D);

P_W = Power_Weight(V_max_kmh,L_D);

[W_S,S] = Wing_Loading(V_stall_kmh,V_cruise_kmh,C_L_max,AR,M_to);

V_cruise_mps = V_cruise_kmh*0.27778;
C_L_cruise = M_to*g/(0.5*V_cruise_mps^2*0.73643*S);


[M_to_refined,M_e_refined,M_bat] = Mass_Iteration2(V_max_kmh,P_W,W_S,AR,f_b,m_payload);

[Length_fuse, Width_fuse,b,C_root,MAC,S_vt,b_vt,Cr_vt,S_ht,b_ht,Cr_ht,b_a,C_a,C_e,b_e,C_r,b_r] = Sizing(M_to_refined,F_rat,Taper_rat,AR,S,T_vt,AR_vt,T_ht,AR_ht,per_a);

%% Output
f = figure;
t = uitable(f);
data = { 'MTOW (kg)',sprintf('%.f',M_to_refined);
        'Wing Area m^2',sprintf('%.2f',S) ;
        'W/S (kg/m^2)',sprintf('%.2f',W_S);
        'P/W (kwh/kg)',sprintf('%.3f',P_W);
        'Cl cruise',sprintf('%.3f',C_L_cruise);
        'Fuselage Length (m)',sprintf('%.3f',Length_fuse);
        'Fuselage Width (m)', sprintf('%.3f',Width_fuse);
        'Wing Span (m)', sprintf('%.3f',b); 
        'C root (m)', sprintf('%.3f',C_root);
        'VT area (m^2)', sprintf('%.3f',S_vt);
        'VT span (m)',sprintf('%.3f',b_vt) ;
        'VT chord (m)', sprintf('%.3f',Cr_vt);
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
t.ColumnWidth = {250, 'auto', 'auto', 250};
t.Position = [25 25 380 600];
t.FontSize = 11; 

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



