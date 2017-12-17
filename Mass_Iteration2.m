function [M_to_refined,M_e_refined,M_bat] = Mass_Iteration2(V_max_kmh,P_W,W_S,AR,f_b,m_payload)
%MASS_ITERATION2 Summary of this function goes here
%   Detailed explanation goes here
W_S_fps =  W_S * 0.20482; % lb/ft^2
P_W_fps=P_W*(1.34102/2.20462); % in hp/lb
V_max_knot = V_max_kmh*0.53995;
Wo_initial =3000*2.20462;
W_payload = m_payload*2.20462;

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


end

