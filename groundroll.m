function [ Sl ] = groundroll( Cl_max, AR, S, W, e0, airfoil, L_fuse, W_fuse)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    [~, rho, ~] = atmosphere(0);
    g = 9.8;
    
    %W=L=0.5*rho*S*V^2
    V_stall = sqrt(W/(0.5*rho*S*Cl_max));
    Vt = 0.7*1.3*V_stall
    
    L = 0.5*rho*Vt^2*S*Cl_max
    
    [Re, M, ~] = nondimensionalize(AR, S, Vt, W, 0);

    %pol = xfoil(airfoil, "cl", Cl_max, Re, M, 'ppar n 200', 'oper iter 1000', ...
    %            'gdes flap 0.8 0.01 35 exec');
    %Cd2d = pol.CD\
    Cd2d =  0.07170;
    Cdf = fuselage_drag(L_fuse, W_fuse, S, 0, Vt);

    %if Cd2d
        Cdi = (Cl_max^2)/(pi*e0*AR);
        Cdw = Cd2d + Cdi;
        Cd = Cdw + Cdf;
        D = 0.5*rho*Vt^2*S*Cd;
        %Pr = sqrt(2*W^3*Cd.^2./(rho*S*Cl_max.^3));
    mu = 0.4;
    Sl = 1.69*W^2/(g*rho*S*Cl_max*(D+mu*(W-L)));
end

