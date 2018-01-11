function [ wing_loading, Pr, min_aoa ] = bestWingLoading(AR, Sp, W, V, e0, h, airfoil)
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
    all_aoa = [];
    for S = Sp
        [Pr, aoa]= powerrequired(AR, S, W, V, e0, h, airfoil);
        if Pr
            allPr = [allPr Pr];
            allS = [allS S];
            all_aoa = [all_aoa aoa];
        end
    end
    idxmin = find(allPr == min(allPr));
    min_aoa = all_aoa(idxmin); 
    wing_loading = W/allS(idxmin);
    Pr = allPr(idxmin);
    plot(allS, allPr,'-o',...
        'DisplayName',sprintf('AR=%d', AR),...
        'MarkerIndices',[idxmin])
   aoa_text = num2str(min_aoa);
   text(allS(idxmin)+0.1, allPr(idxmin)+0.1, aoa_text);
end

