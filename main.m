clear all;
close all;

g = 9.791;
e0 = 0.8; % Wing efficiency
W = 6147.8*g; % N (weight)
W_bat = 1000; % kg
V = 450/3.6; % m/s (cruise velocity)

%% costs
Q = linspace(300,1000);
RDTEflywaway = costs(Q, W/g, W_bat, V*3.6);

figure();
plot(Q, RDTEflywaway./Q/1e6)

xlabel('Aircraft produced')
ylabel('Cost per plane [Million USD]')
print -depsc planecosts

%% airfoil selection
AR = 9;
S = 30;
[Re, M, Cl] = nondimensionalize(AR, S, V, W);
airfoil = bestAirfoil(Cl, Re, M);
%airfoil = 'coord_seligFmt/naca663418.dat';

%% wing loading
figure();
hold on;
for AR = 8:12
    [wingloading, Pr] = bestWingLoading(AR, 20:40, W, V, e0, airfoil)
end

title("Power required at cruise")
xlabel("Wing area [m^2]")
ylabel("Power required [Watt]")
legend('show')
grid on;

print -depsc powerrequired