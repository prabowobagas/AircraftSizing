function [ Pr, aoa, Cd] = powerrequired( AR, S, W, V, e0, h, airfoil, L_fuse, W_fuse)
%powerrequired calculate the powe required for cruise flight
%   AR aspect ratio
%   S  wing area (m^2)
%   W  weight (newton)
%   V  cruise velocity (m/s)
%   e0 wing efficincy
%   h  altitude (m)

    [T, rho, mu] = atmosphere(h);
    [Re, M, Cl] = nondimensionalize(AR, S, V, W, h);

    pol = xfoil(airfoil, "cl", Cl, Re, M, 'ppar n 200', 'oper iter 100');
    aoa = pol.alpha;
    %pol = xfoil('NACA2412', Cl, Re, M, 'ppar n 300', 'oper iter 1000');
    Cd2d = pol.CD;
    Cdf = fuselage_drag(L_fuse, W_fuse, S, h, V)

    if Cd2d
        Cdi = (Cl^2)/(pi*e0*AR);
        Cdw = Cd2d + Cdi;
        Cd = Cdw + Cdf;
        Dw = 0.5*rho*V^2*S*Cd;
        Pr = sqrt(2*W^3*Cd.^2./(rho*S*Cl.^3));
    else
        Pr = [];
    end
end

