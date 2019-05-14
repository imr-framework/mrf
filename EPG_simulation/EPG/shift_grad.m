function omega_new = shift_grad(delk,omega)
%SHIFT_GRAD(delk,omega)
% Modified by: Gehua Tong
% Date: 15 Oct 2018
% Modified by: Sachin A B Anchan
% Date: 30 June 2014
% Shift applies to only F+ and F-* as it does not dephase in z
% check size of previous omega to determine the effect - test multiple
% times

% delk: integer indicating discrete change in k
% omega: inputted omega matrix (with columns of [F+,F-,Z]')

[m,n] = size(omega); %previous time point
% if(m~=3)
%     error('Still implementing equation 26, please use 3xk');
% end
if delk == 0
    omega_new = omega;
else
if(n>1) % typical case: an RF pulse has happened and we have transverse components
    F = [fliplr(omega(1,:)) squeeze((omega(2,2:end)))]; %arrange to make it like eq 27
        % Negative shift
        if(delk < 0)
              F = [zeros(1,abs(delk)) F]; %negative shift moves the population downwards - 2n-1 + delk
              Z = [squeeze(omega(3,:)) zeros(1,abs(delk))]; %No change in z due to grads
              Fp = [fliplr(F(1:n)) zeros(1,abs(delk))]; 
              Fm = F(n:end);
              % Here, V(k=1) moves into V'(k=+0), 
              %       so V'(k=-0) is the conjugate of V'(k=+0) 
              Fm(1) = conj(Fm(1));
        % Positive shift
        else
              F = [F zeros(1,delk)]; %positive shift pushes the population upwards
              Z = [squeeze(omega(3,:)) zeros(1,delk)];
              Fp = fliplr(F(1:n+delk)); 
              Fm = [F(n+delk:end) zeros(1,delk)];
              % Here, V(k=-1) moves into V'(k=-0),
              %       so V'(k=+0) is the conjugate of V'(k=-0)
              Fp(1) = conj(Fp(1));
        end


else % n = 1; this happens if pulse sequence starts with nonzero transverse
     %        components and no RF pulse at t = 0 - gradient happens first
     % omega(1) (=F+(0)) and omega(2)(= F-(0)) must be complex conjugates!!!
     % omega(3) = Z(0)
     if(delk > 0)
       Fp = [zeros(1,abs(delk)) omega(1)];
       Fm = [0  zeros(1,abs(delk))];
        Z = [squeeze(omega(3)) zeros(1,abs(delk))];
     else
        Fp = [0  zeros(1,abs(delk))];
        Fm =  [zeros(1,abs(delk)) omega(2)];
        Z = [squeeze(omega(3)) zeros(1,abs(delk))];
        
    end
end

omega_new = [Fp;Fm;Z];
end
