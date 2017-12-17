function P_W = Power_Weight(V_max_kmh,L_D)
eta = 0.8;
V_max_fps = 0.91134*V_max_kmh; %Velocity feet per second

%V_max_mps = 0.277778*V_max_kmh; %Velocity meter per second
%V_cruise_mps = V_cruise_kmh*0.27778;
%V_stall_mps = 0.277778*V_stall_kmh;
%V_stall_fps = V_stall_mps*3.28084;

P_W_stat = 0.016*V_max_kmh^0.5; % Watt/g for turboprop
T_W_cruise = 1/L_D;
T_W_to = T_W_cruise*1/0.6;
P_W_cruise_fps = T_W_to/(eta*550/V_max_fps); % P_W required for cruise at retard units
P_W_cruise = P_W_cruise_fps*(0.74570/0.45359); %Convert to metric kW/kg

if P_W_stat>P_W_cruise 
    P_W = P_W_stat;
else
    P_W = P_W_cruise;
end
end

