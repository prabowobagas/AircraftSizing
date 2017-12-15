clear all;
close all;
clc; 
%% Parameters for planes
AR = 9.2;
L_D = 16.5; %L/D Guess
m_payload=1000; %Set payload weight
Range = 1000;% range in km

V_max_kmh = 500; %  maximum velocity not cruise (km/h)
V_cruise_kmh = 450;
V_stall_kmh = 140; 
C_L_max = 2; %With flaps used



%% Variables 
g = 9.8;
eta = 0.8; %efficiency
E = 3.6*10^6; %Energy Density in Joules
initial_m = 3000;
R = Range*10^3;
%% Iteration
f_e = 1.4*initial_m^(-0.1); %inital sizing from 3000 kg
f_b = 1.05*(R*g)/(eta.*L_D.*E); %Battery mass fraction

fo_previous = zeros(12,1);
fe_previous = zeros(12,1);
fe_calc = zeros(12,1); %Calculated empty mass
fo_calc = zeros(12,1); %Calculated MTOW 

fo_previous(1,1) = initial_m; %Start with 3000
fe_previous(1,1) = f_e;
fo_calc(1,1) = (m_payload)./(1-fe_previous(1,1)-f_b);
fe_calc(1,1) = fe_previous(1,1).*fo_calc(1,1);

for i=2:12
    fo_previous(i,1) = (fo_previous(i-1,1)+fo_calc(i-1,1))/2;
    fe_previous(i,1) = 1.4*fo_previous(i,1)^(-0.1);
    fo_calc(i,1) = m_payload/(1-f_b-fe_previous(i,1));
    fe_calc(i,1) = fe_previous(i,1)*fo_calc(i,1);
end
M_to = fo_calc(12,1); 

%%
%% Conversions
V_max_fps = 0.91134*V_max_kmh; %Velocity feet per second
V_max_mps = 0.277778*V_max_kmh; %Velocity meter per second
V_cruise_mps = V_cruise_kmh*0.27778;
V_stall_mps = 0.277778*V_stall_kmh;
V_stall_fps = V_stall_mps*3.28084;
%% Thrust to Weight Specification (Power to Weight for Propeller), Pick Highest P/W ratio
P_W_stat = 0.016*V_max_kmh^0.5; % Watt/g for turboprop
T_W_cruise = 1/L_D;
T_W_to = T_W_cruise*1/0.6;
P_W_cruise_fps = T_W_to/(eta*550/V_max_fps); % P_W required for cruise at retard units
P_W_cruise = P_W_cruise_fps*(0.74570/0.45359); %Convert to metric kW/kg
if P_W_stat>P_W_cruise % Maximum P/W case
    P_W = P_W_stat;
else
    P_W = P_W_cruise;
end
P_W_fps=P_W*(1.34102/2.20462); % in hp/lb

%% Wing Loading Specification (Pick Lowest W/S), Different scenario of W/S
W_S_approx = 195*9.8; %W/S ratio kg/m^2 
W_S_stall = 1/2 * 1.23 * V_stall_mps^2 * C_L_max;

W_S_cruise = 1/2*0.69747*V_cruise_mps^2 * (pi*AR*0.8*0.024)^0.5;

W_S_dat = [W_S_approx W_S_stall W_S_cruise];
W_S = min(W_S_dat)/g; %(kg/m^2)
S = (M_to)/W_S; % Area of the wing for needed W/S (m^2)
W_S_fps =  W_S * 0.20482; % lb/ft^2
T_W = (P_W*0.8)/(1.2*V_stall_mps*0.7); %Thrust to weight conversion


%Lift coefficient for C_L_Cruise
C_L_cruise = M_to*g/(0.5*V_cruise_mps^2*0.73643*S);
S_landing = 5*W_S_stall/9.8*1/C_L_max + 183; %landing distance


% For Turboprop, Historical values are :
% P/W = 0.33 (kW/kg)
% W/S = 195 (kg/m^2)
% AR = 9.2

%% For reference Beech C90GTx
% Aspect Ratio = 9.7
% Wing Area S = 27.41
% MTOW = 4756 kg
% V cruise max = 504 km/h
% V stall = 140 km/h
% TO distance = 605 m
% Landing distance = 640 m 


%% Variables
V_max_knot = V_max_kmh*0.53995;
Wo_initial =initial_m*2.20462;
W_payload = m_payload*2.20462;

%% W_e/W_o fraction refined (Using P/W and W/S from previous script) Based on stats

% Only for empirical sizing 
% Statistical fit only accepts retard units and only for Turboprop
We_initial = 0.37 + 0.09*(Wo_initial)^0.06 * AR^0.08 * P_W_fps^0.8 * W_S_fps^-0.05 * V_max_knot.^0.30; 

Wo_previous = zeros(12,1);
We_previous = zeros(12,1);
We_calc = zeros(12,1); %Calculated empty mass
Wo_calc = zeros(12,1); %Calculated MTOW 

Wo_previous(1,1) = Wo_initial; 
We_previous(1,1) = We_initial;
Wo_calc(1,1) = (W_payload)/(1-f_b-We_previous(1,1));
We_calc(1,1) = We_previous(1,1).*Wo_calc(1,1);

for j=2:12
    Wo_previous(j,1) = (Wo_previous(j-1,1) + Wo_calc(j-1,1))/2;
    We_previous(j,1) = 0.37 + 0.09*(Wo_previous(j,1))^0.06 * AR^0.08 * P_W_fps^0.8 * W_S_fps^-0.05 * V_max_knot.^0.30;
    Wo_calc(j,1) = W_payload/(1-f_b-We_previous(j,1));
    We_calc(j,1) = We_previous(j,1)*Wo_calc(j,1);
end
M_to_refined = Wo_calc(12,1)*0.45359; %Convert back to kg
M_e_refined = We_calc(12,1)*0.45359;
M_bat = f_b*M_to_refined;
S_takeoff = (1.44*M_to_refined^2)/(g*1.23*S*C_L_max*T_W*M_to_refined)


%% Variables
%Fineness ratio of 6 to 8 (Fineness ratio => length of fuselage/maximum width
%of fuselage.
Taper_rat = 0.4; 
F_rat = 6; %Fineness Ratio
%% Initial estimation sizing on Fuselage, Wing, Tail Volume, Control Surfaces.

Length_fuse = 0.169*(M_to_refined)^0.51; %Fuselage Length based on statistics (m)
Width_fuse = F_rat/Length_fuse; % Fuselage Width (m)
b = (AR*S)^0.5; %Wing span
C_root = (2*S)/(b*(1+Taper_rat)); %Length of Root
MAC = 2/3 * C_root * (1+ Taper_rat + Taper_rat^2)/(1 + Taper_rat); %Mean Aerodynamic Chord
S_vt = (0.08*b*S)/(0.55*Length_fuse); %Area of vertical tail
S_ht = (0.9*MAC*0.25*S)/(0.55*Length_fuse); %Area of horizontal tail

T = table(M_to,S,W_S,P_W,S_landing,C_L_cruise,M_to_refined);

T.Properties.VariableNames = {'MTOW_Initial' 'Wing_Area' 'Wing_Loading' 'Power_to_Weight' 'Landing_Distance' 'Cruise_Lift_Coef' 'MTOW_Final'}

T2 = table(Length_fuse,Width_fuse,b,C_root,S_vt,S_ht);
T2.Properties.VariableNames = {'Fuse_Length' 'Fuse_Width' 'Wing_Span' 'Root_Length' 'Vertical_Stabilizier_Area' 'Horizontal_Stabilizer_Area'}



Q = linspace(300,1000); % production quantity
FTA = 2; % number of flight test aircraft (typical 2-6)
Neng = Q*2; % total production quantity times number of engines per aircraft
Tmax = 300; % engine maximum thrust (kw)
Ca = 60000*Q; % avionics cost

He = 5.18*M_e_refined^0.777*V_max_kmh^0.894*Q.^0.163; % Engineering hours
Ht = 7.22*M_e_refined^0.777*V_max_kmh^0.969*Q.^0.263; % tooling hours
Hm = 10.5*M_e_refined^0.820*V_max_kmh^0.484*Q.^0.641; % Mfg hours
Hq = 0.076*Hm; % QC hours if cargo
Cd = 67.4*M_e_refined^0.630*V_max_kmh^1.3; % Devel support cost
Cf = 1947*M_e_refined^0.325*V_max_kmh^0.822*FTA^1.21; % Flt test cost
Cm = 31.2*M_e_refined^0.921*V_max_kmh^0.621*Q.^0.79; %  Mfg materials cost
Ce = 50*Tmax; % engine cost
Kwh_cost = M_bat*1000*0.5;

Re = 115;
Rt = 118;
Rm = 98;
Rq = 108;

RDTEflywaway = He*Re+Ht*Rt+Hm*Rm+Hq*Rq+Cd+Cf+Ce*Neng+Ca+Kwh_cost.*Q;
plot(Q, RDTEflywaway./Q/1e6)

%title()
xlabel('Aircraft produced')
ylabel('Cost per plane [Million USD]')
print -depsc planecosts

