function P_W = Power_Weight(V_max_kmh,L_D)
%Power to Weight 
%   INPUT:  Maximum Velocity (km/h)
%           L/D
%   OUTPUT: P/W (Kwh/Kg)

eta = 0.8;
V_max_fps = 0.91134.*V_max_kmh; %Velocity feet per second

P_W_stat = 0.016.*V_max_kmh.^0.5; % Watt./g for turboprop
T_W_cruise = 1./L_D;
T_W_to = T_W_cruise.*1./0.6;
P_W_cruise_fps = T_W_to./(eta.*550./V_max_fps); % P_W required for cruise at retard units
P_W_cruise = P_W_cruise_fps.*(0.74570./0.45359); %Convert to metric kW./kg

if P_W_stat > P_W_cruise %% Limitation: Doesn't work if vector
    P_W = P_W_stat;
else
    P_W = P_W_cruise;
end
end

