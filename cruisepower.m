clear all;
close all;

rho = 0.73643; % Pa (air density at 5km)
mu = 1.628e-5; % N s/m^2
T = 255.69; % Kelving (temperature at 5km)
g = 9.791;
 
e0 = 0.8; % Wing efficiency
W = 6147.8*g; % N (weight)
V = 450/3.6; % m/s (cruise velocity)
S = 30; % m^2 (wing area)
AR = 9; % aspect ratio

b = sqrt(AR*S); % m (span)
c = sqrt(S/AR); % m (mean cord lenght)        
Re = rho*V*c/mu;
M = V/320;
Cl = W/(0.5*rho*V^2*S);

%% find best profile

airfoils = dir('coord_seligFmt/*.dat');
alpha = zeros(0, length(airfoils));
Cd = zeros(0, length(airfoils));
for coord = 1:length(airfoils)
    coord
    str = strcat('coord_seligFmt/', airfoils(coord).name);
    try
        pol = xfoil(str, Cl, Re, M, 'ppar n 200', 'oper iter 200');
    catch
        continue
    end
    if pol.alpha
        alpha(coord) = pol.alpha;
        Cd(coord) = pol.CD;
    end
end
 
%% plot results
myalpha = alpha(Cd~=0);
myairfoils = airfoils(Cd~=0);
myCd = Cd(Cd~=0);

outliers = isoutlier(myCd);
outlierCd = myCd(outliers);
outlieralpha = myalpha(outliers);
myalpha(outliers) = [];
myCd(outliers) = [];
myairfoils(outliers) = [];

scatter(myalpha, myCd)
title("Drag of various airfoils at cruise condition")
xlabel("Alpha")
ylabel("Cd")
print -depsc airfoilsummary

%% find porfile with lowest drag

[~,Cdidx] = sort(myCd,'ascend');
topNames = {myairfoils(Cdidx(1:10)).name}';
topCd = myCd(Cdidx(1:10))';
topAlpha = myalpha(Cdidx(1:10))';
table(topNames, topCd, topAlpha)

topAirfoil = 'coord_seligFmt/naca663418.dat';

%% find optimal wing loading
hold on;
for AR = 8:12
    Pr = [];
    usedS = [];
    for S = 20:40
        b = sqrt(AR*S); % m (span)
        c = sqrt(S/AR); % m (mean cord lenght)        
        Re = rho*V*c/mu;
        M = V/320;

        Cl = W/(0.5*rho*V^2*S);

        %pol = xfoil(topAirfoil, Cl, Re, M, 'ppar n 300', 'oper iter 1000');
        pol = xfoil('NACA2412', Cl, Re, M, 'ppar n 300', 'oper iter 1000');
        Cd2d = pol.CD;
        if Cd2d
            Cd = Cd2d + (Cl^2)/(pi*e0*AR);
            Pr = [Pr sqrt(2*W^3*Cd.^2./(rho*S*Cl.^3))];
            usedS = [usedS S];
        end
    end
    AR
    Pr
    usedS
    idxmin = find(Pr == min(Pr));
    plot(usedS, Pr,'-o',...
        'DisplayName',sprintf("AR=%d", AR),...
        'MarkerIndices',[idxmin])
end

title("Power required at cruise")
xlabel("Wing area [m^2]")
ylabel("Power required [Watt]")
legend('show')
grid on;

print -depsc powerrequired