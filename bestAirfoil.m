function [ nameacc, Sacc, Pracc ] = bestAirfoil(AR, S0, W, V, e0, h)
%bestAirfoil Let xfoil calculate minimum drag
%   Loops over all NACA 6 airfoils in thcoord_seligFmt
%   Running xfoil on them for a given Cl
%   Returns the airfoil with the lowest Cd

    [T, rho, mu] = atmosphere(h);
    [Re, M, Cl] = nondimensionalize(AR, S0, V, W, h);

    %airfoils = [dir('coord_seligFmt/naca6*.dat'); dir('coord_seligFmt/n6*.dat')] ;
    airfoils = dir('coord_seligFmt/*.dat');
    nameacc = cell(0, length(airfoils));
    Pracc = zeros(0, length(airfoils));
    Sacc = zeros(0, length(airfoils));
    
    for coord = 1:length(airfoils)
        str = strcat('coord_seligFmt/', airfoils(coord).name)
        try
            pol = xfoil(str, 'alfa', 0:0.5:5, Re, M, 'ppar n 200', 'oper iter 500');
        catch
            continue
        end
        if ~isempty(pol.alpha)
            Cd3d = pol.CD + (pol.CL.^2)/(pi*e0*AR);
            LD = pol.CL./Cd3d;
            [~, idx] = max(LD);

            Cl = pol.CL(idx);
            Cd = Cd3d(idx);
            alpha = pol.alpha(idx);

            S = W/(0.5*rho*V^2*Cl);

            Pracc(coord) = sqrt(2*W^3*Cd.^2./(rho*S*Cl.^3));
            Sacc(coord) = S;
            nameacc(coord) = cellstr(pol.name);
        end
    end
end

