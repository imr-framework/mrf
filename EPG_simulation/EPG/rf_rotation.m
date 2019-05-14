function T_phi_alpha = rf_rotation (phi,alpha)
%RF_ROTATION(phi,alpha)
% 3 x 3 rotation matrix in complex representation for RF pulse
% with phase = phi (from the x axis) and filp angle = alpha
% See Weigel 2015, Eq. 15

phi = deg2rad(phi);
alpha = deg2rad(alpha);

T_phi_alpha = [cos(alpha/2).^2                     exp(2*1i*phi)*sin(alpha/2).^2       -1i*exp(1i*phi)*sin(alpha);
           exp(-2*1i*phi)*sin(alpha/2).^2       cos(alpha/2).^2                     1i*exp(-1i*phi)*sin(alpha);
           -1i*0.5.*exp(-1i*phi)* sin(alpha)    1i*0.5.*exp(1i*phi)* sin(alpha)     cos(alpha)];