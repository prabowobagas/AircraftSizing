function [M_to,f_b,M_e] = Mass_Iteration1(Range,m_payload,L_D)
%MASS_ITERATION1 
%   INPUT:  Range (km)
%           Payload Mass (kg)
%           L/D
%   OUTPUT: MTOW (kg)
%           Battery Mass fraction Mb/Mo 
%           Empty Mass (Kg)

R = Range.*10.^3;
E = 3.6.*10.^6; %Energy Density in Joules
g = 9.8;
eta = 0.8; %efficiency
f_e = 1.4.*3000.^(-0.1); %inital sizing from 3000 kg
f_b = 1.05.*(R.*g)./(eta.*L_D.*E); %Battery mass fraction


fo_previous = zeros(12,1);
fe_previous = zeros(12,1);
fe_calc = zeros(12,1); %Calculated empty mass
fo_calc = zeros(12,1); %Calculated MTOW 

fo_previous(1,1) = 3000; %Start with 3000
fe_previous(1,1) = f_e;
fo_calc(1,1) = (m_payload)./(1-fe_previous(1,1)-f_b);
fe_calc(1,1) = fe_previous(1,1).*fo_calc(1,1);

for i=2:12
    fo_previous(i,1) = (fo_previous(i-1,1)+fo_calc(i-1,1))./2;
    fe_previous(i,1) = 1.4.*fo_previous(i,1).^(-0.1);
    fo_calc(i,1) = m_payload./(1-f_b-fe_previous(i,1));
    fe_calc(i,1) = fe_previous(i,1).*fo_calc(i,1);
end
M_to = fo_calc(12,1); 
M_e = fe_calc(12,1);
end

