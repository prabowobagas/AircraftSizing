function [ wing_loading, Pr ] = bestWingLoading(AR, Sp, W, V, e0, h, airfoil)
%bestWingLoading 
%   AR aspect ratio
%   Sp vector of proposed areas to check ([m^2])
%   W  weight (N)
%   V  Velocity (m/s)
%   e0 span efficiency
%
%   wing_loading (N/m^2)
%   Pr power required (N)

    allPr = [];
    allS = [];
    for S = Sp
        Pr = powerrequired(AR, S, W, V, e0, h, airfoil);
        if Pr
            allPr = [allPr Pr];
            allS = [allS S];
        end
    end
    idxmin = find(allPr == min(allPr));
    wing_loading = W/allS(idxmin);
    Pr = allPr(idxmin);
    plot(allS, allPr,'-o',...
        'DisplayName',sprintf('AR=%d', AR),...
        'MarkerIndices',[idxmin])
end

