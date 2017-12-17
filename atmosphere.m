function [ T, rho, mu ] = atmosphere(h)
%atmosphere

    g0 = 9.81; %m/s
    p0 = 101325; %Pa
    T0 = 288.16; %K
    rho0 = 1.225; %kg/m^3
    R = 287; %J/kgK
    a1 = -6.5e-3; %K/m

    if h <= 11000
        T = T0 + a1.*h;

        %p1 = p0.*(T1./T0).^(-g0/(a1*R));
        rho = rho0*(T/T0).^(-g0/(a1*R)-1);
    else
        T = T0 + a1.*11000;
        rho1 = rho0.*(T./T0).^(-g0/(a1*R)+1);
        
        %p = p1(end)*exp(-g0/(R*T(end))*(h2-h1(end)));
        rho = rho1*exp(-g0/(R*T)*(h-11000));
    end
    
    T0 = 291.15;
    C = 120;
    mu0 = 18.27e-6;
    mu = mu0*(T/T0)^(3/2)*((T0+C)/(T+C));
end

