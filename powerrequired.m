function [ Pr ] = powerrequired( AR, S, W, V, e0, h, airfoil)
%powerrequired calculate the powe required for cruise flight
%   AR aspect ratio
%   S  wing area (m^2)
%   W  weight (newton)
%   V  cruise velocity (m/s)
%   e0 wing efficincy
%   h  altitude (m)

    [T, rho, mu] = atmosphere(h);
    [Re, M, Cl] = nondimensionalize(AR, S, V, W, h);

    pol = xfoil(airfoil, Cl, Re, M, 'ppar n 300', 'oper iter 100');
    %pol = xfoil('NACA2412', Cl, Re, M, 'ppar n 300', 'oper iter 1000');
    Cd2d = pol.CD;
    if Cd2d
        % fuselage drag currently missing!!!
        Cd = Cd2d + (Cl^2)/(pi*e0*AR);
        Pr = sqrt(2*W^3*Cd.^2./(rho*S*Cl.^3));
    else
        Pr = [];
    end
end

