function [Length_fuse, Width_fuse,b,C_root,MAC,S_vt,b_vt,Cr_vt,S_ht,b_ht,Cr_ht,b_a,C_a,C_e,b_e,C_r,b_r] = Sizing(M_to_refined,F_rat,Taper_rat,AR,S,T_vt,AR_vt,T_ht,AR_ht,per_a)
%SIZING Summary of this function goes here
%   Detailed explanation goes here
%
%
%
%
%
%
%
%
%

Length_fuse = 0.169*(M_to_refined)^0.51; %Fuselage Length based on statistics (m)
Width_fuse = 1/(F_rat/Length_fuse); % Fuselage Width (m)
b = (AR*S)^0.5; %Wing span
C_root = (2*S)/(b*(1+Taper_rat)); %Length of Root
MAC = 2/3 * C_root * (1+ Taper_rat + Taper_rat^2)/(1 + Taper_rat); %Mean Aerodynamic Chord

S_vt = (0.08*0.95*b*S)/(0.55*Length_fuse); %Area of vertical tail
b_vt = (AR_vt*S_vt)^0.5; %Span of vertical and horizontal tail
Cr_vt = (2*S_vt)/(b_vt*(1+T_vt)); %chord length of vt and ht

S_ht = (0.9*0.95*MAC*0.25*S)/(0.55*Length_fuse); %Area of horizontal tail
b_ht = (AR_ht*S_ht)^0.5;
Cr_ht = (2*S_ht)/(b_ht*(1+T_ht));

b_a = per_a * b; %Aileron span assumed to be 50% of whole span 
C_a = 0.18 * C_root; %Aileron chord

C_e = 0.36 * Cr_ht; 
b_e = 0.9 * b_ht; 

C_r = 0.46 * Cr_vt; 
b_r = 0.9 * b_vt; 


end

