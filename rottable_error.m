function [ d_phi,d_theta,d_psi ] = rottable_error( f_input,omega_input )
%ROTTABLE_ERROR Function
%   Die Parameter muessen bereits gemittelt sein, damit eine Berechnung
%   durchgefuehrt werden kann.
%
% [d_phi,d_theta]=rottable_error(f_input,omega_input)
%
%   f_input=[fx,fy,fz]                      =>G_vektoren
%   omega_input=[omega_x,omega_y,omega_z]   =>Earth_vektoren
%   Rueckgabewert ist in [°]
if nargin==1
    d_psi=0;
f_abs=sqrt(f_input(1)^2+f_input(2)^2+f_input(3)^2);
d_phi=atand(-f_input(2)/-f_input(3));
d_theta=asind(f_input(1)/f_abs);
elseif nargin==2
f_abs=sqrt(f_input(1)^2+f_input(2)^2+f_input(3)^2);
d_phi=atand(-f_input(2)/-f_input(3));
d_theta=asind(f_input(1)/f_abs);
d_psi=atand((-omega_input(2)*cosd(d_phi)+omega_input(3)*sind(d_phi))/...
    (omega_input(1)*cosd(d_theta)+omega_input(2)*sind(d_phi)*sind(d_theta)+...
    omega_input(3)*cosd(d_phi)*sind(d_theta)));
end
end

