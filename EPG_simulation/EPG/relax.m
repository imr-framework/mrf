function omega_new = relax(tau,T1,T2,omega)
%relax(tau,T1,T2,omega)
%relax.m: updates omega with relaxation effects
% INPUTS
%   tau: duration of relaxation in ms
%   T1,T2: time constants in ms
%   omega: the input (3 x n) matrix of k-states

% Gehua Tong, Oct 7 2018

if size(omega,1) ~= 3
    error('Size of k-state matrix incorrect. Input needs to be (3 x n)')
end
if T1 ~= 0 && T2 ~=0
    E1 = exp(-tau/T1);
    E2 = exp(-tau/T2);
    Emat = [E2  0  0;
             0 E2  0;
             0  0 E1];
    omega_new = Emat*omega;
    omega_new(:,1) = omega_new(:,1) + [0,0,1-E1]'; % M0 = 1 default
else
    omega_new = omega;
end

end
 
