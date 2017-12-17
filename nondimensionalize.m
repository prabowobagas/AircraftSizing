function [ Re, M, Cl ] = nondimensionalize( AR, S, V, W, h)
%nondimensionalize does what it says on the tin

    [T, rho, mu] = atmosphere(h);
    
    g = 9.791;

    b = sqrt(AR*S); % m (span)
    c = sqrt(S/AR); % m (mean cord lenght)        
    Re = rho*V*c/mu;
    M = V/320;

    Cl = W/(0.5*rho*V^2*S);
end

