function [W_S, S] = Wing_Loading(V_stall_kmh,V_cruise_kmh,C_L_max,AR,M_to)
%WING_LOADING 
%   INPUT:  V stall (km/h)
%           V cruise (km/h)
%           Cl max at takeoff
%           AR aspect ratio
%           MTOW (kg)
%   OUTPUT: W/S ratio (kg/m^2)
%           S Surface area of wing

rho_sealevel=1.225;
rho_5000m = 0.69747; 
V_cruise_mps = V_cruise_kmh*0.27778;
V_stall_mps = 0.277778*V_stall_kmh;
e0=0.8;
g=9.8;

W_S_approx = 195.*g; %W/S ratio from statistics 
W_S_stall = 1/2 .* rho_sealevel .* V_stall_mps.^2 .* C_L_max; %rho
W_S_cruise = 1/2.*rho_5000m.*V_cruise_mps.^2 .* (pi.*AR.*e0.*0.024).^0.5;

W_S_dat = [W_S_approx W_S_stall W_S_cruise];
W_S = min(W_S_dat)./g; %(kg/m^2)
S = (M_to)./W_S; % Area of the wing for needed W/S (m^2)
end

