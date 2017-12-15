clear all;
close all;
hold on;

rho = 0.73643; % Pa (air density at 5km)
mu = 1.628e-5; % N s/m^2
T = 255.69; % Kelving (temperature at 5km)
g = 9.791;
 
e0 = 0.8; % Wing efficiency
W = 6147.8*g; % N (weight)
V = 450/3.6; % m/s (cruise velocity)
S = 20; % m^2 (wing area)
AR = 10; % aspect ratio

for AR = 9:12
    Pr = [];
    for S = 20:40
        b = sqrt(AR*S); % m (span)
        c = sqrt(S/AR); % m (mean cord lenght)        
        Re = rho*V*c/mu;
        M = V/320;

        Cl = W/(0.5*rho*V^2*S);

        pol = xfoil('NACA2412', Cl, Re, M, 'oper iter 1000');
        Cd2d = pol.CD;

        Cd = Cd2d + (Cl^2)/(pi*e0*AR);
        Pr = [Pr sqrt(2*W^3*Cd.^2./(rho*S*Cl.^3))];
    end
    
    idxmin = find(Pr == min(Pr));
    plot(20:40, Pr,'-o',...
        'DisplayName',sprintf("AR=%d", AR),...
        'MarkerIndices',[idxmin])
end

title("Power required at cruise")
xlabel("Wing area [m^2]")
ylabel("Power required [Watt]")
legend('show')
grid on;

print -depsc powerrequired