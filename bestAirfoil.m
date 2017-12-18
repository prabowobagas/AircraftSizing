function [ airfoil ] = bestAirfoil( Cl, Re, M )
%bestAirfoil Let xfoil calculate minimum drag
%   Loops over all NACA 6 airfoils in thcoord_seligFmt
%   Running xfoil on them for a given Cl
%   Returns the airfoil with the lowest Cd

airfoils = dir('coord_seligFmt/naca6*.dat');
alpha = zeros(0, length(airfoils));
Cd = zeros(0, length(airfoils));
for coord = 1:length(airfoils)
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

%outliers = isoutlier(myCd);
%outlierCd = myCd(outliers);
%outlieralpha = myalpha(outliers);
%myalpha(outliers) = [];
%myCd(outliers) = [];
%myairfoils(outliers) = [];

figure();
scatter(myalpha, myCd)
title('Drag of various airfoils at cruise condition')
xlabel('Alpha')
ylabel('Cd')
print -depsc airfoilsummary

%% find porfile with lowest drag

[~,Cdidx] = sort(myCd,'ascend');
topNames = {myairfoils(Cdidx(1:10)).name}';
topCd = myCd(Cdidx(1:10))';
topAlpha = myalpha(Cdidx(1:10))';
table(topNames, topCd, topAlpha)

%topAirfoil = 'coord_seligFmt/naca663418.dat';
airfoil = strcat('coord_seligFmt/', myairfoils(Cdidx(1)).name);
end

