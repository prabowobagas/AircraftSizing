function [ RDTEflywaway ] = costs( Q, M, M_bat, V)
%costs Estimate plane costs
%   Q  quantity
%   M  mass (kg)
%   M_bat battery mass (kg)
%   V speed (km/h)

FTA = 2; % number of flight test aircraft (typical 2-6)
Neng = Q*2; % total production quantity times number of engines per aircraft
Tmax = 300; % engine maximum thrust (kw)
Ca = 60000*Q; % avionics cost

He = 5.18*M^0.777*V^0.894*Q.^0.163; % Engineering hours
Ht = 7.22*M^0.777*V^0.969*Q.^0.263; % tooling hours
Hm = 10.5*M^0.820*V^0.484*Q.^0.641; % Mfg hours
Hq = 0.076*Hm; % QC hours if cargo
Cd = 67.4*M^0.630*V^1.3; % Devel support cost
Cf = 1947*M^0.325*V^0.822*FTA^1.21; % Flt test cost
Cm = 31.2*M^0.921*V^0.621*Q.^0.79; %  Mfg materials cost
Ce = 50*Tmax; % engine cost
Kwh_cost = M_bat*1000*0.5;

Re = 115;
Rt = 118;
Rm = 98;
Rq = 108;

RDTEflywaway = He*Re+Ht*Rt+Hm*Rm+Hq*Rq+Cd+Cf+Ce*Neng+Ca+Kwh_cost.*Q;
end

