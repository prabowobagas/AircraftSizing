function [ Re, M, Cl ] = nondimensionalize( AR, S, V, W)
%nondimensionalize does what it says on the tin

    rho = 0.73643; % Pa (air density at 5km)
    mu = 1.628e-5; % N s/m^2
    T = 255.69; % Kelving (temperature at 5km)
    g = 9.791;

    b = sqrt(AR*S); % m (span)
    c = sqrt(S/AR); % m (mean cord lenght)        
    Re = rho*V*c/mu;
    M = V/320;

    Cl = W/(0.5*rho*V^2*S);
end

