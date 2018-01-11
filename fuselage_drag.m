function [ Cd ] = fuselage_drag( C_fuse, D_fuse, S_wing, h )
    [T, rho, mu] = atmosphere(h);
    Re=(rho*139*C_fuse)/mu;  % all TURBULENT SAME AS CF_total
    Cf_fuse=0.074/(Re^(1/5));

    A=pi*D_fuse*C_fuse;
    
    Cd = A/S_wing*Cf_fuse*2;  % Total drag of fuselage in Nwe
end

