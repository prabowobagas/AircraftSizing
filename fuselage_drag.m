function [ Cd ] = fuselage_drag( C_fuse, D_fuse, S_wing, h, V_cruise)
    [T, rho, mu] = atmosphere(h);
    Re=(rho*V_cruise*C_fuse)/mu;  % all TURBULENT SAME AS CF_total
    Cf_fuse=0.074/(Re^(1/5));

    A=pi*D_fuse*C_fuse;
    
    Df = 0.5*rho*V_cruise^2*A*Cf_fuse*2;
    
    Cd = A/S_wing*Cf_fuse*2; % Total drag
end

