function [B1,A,w1] = make_HSI(A0,beta,mu,t, disp)
% This script implements a hyperbolic secant (HSI) pulse introduced in [1].
%         A0  Maximum pulse amplitude (rad/s)
%       beta  Constants
%         mu  Constants
%          t  Time
%
% OUTPUT
%         w1  Frequency as a function of time
%          A  Amplitude as a function of time
%         B1  B1 field as a function of time

% Created by Sairam
% 9/2019 Modified by Enlin Qian
% # Copyright of the Board of Trustees of Columbia University in the City of New York
% [1] Garwood, M., & DelaBarre, L. (2001). The Return of the Frequency Sweep: 
% Designing Adiabatic Pulses for Contemporary NMR. Journal Of Magnetic Resonance, 
% 153(2), 155-177. doi: 10.1006/jmre.2001.2340

A = A0.*sech(beta.*t); % Amplitude modulation function
w1 = -mu.*beta.*tanh(beta.*t);
B1 = A.*exp(-1i*w1.*t);

%%
if(disp)
figure(211); subplot(511); stem(t,A,'k'); xlabel('time(ms)'); ylabel('AM envelope'); hold all;
subplot(512); plot(t,w1,'k');xlabel('time(ms)'); ylabel('FM envelope');hold all;
subplot(513);plot(t,real(B1)); hold on; plot(t, imag(B1)); xlabel('time(ms)'); ylabel ('B1'); legend('B1x','B1y');
subplot(514);plot(real(B1), imag(B1)); xlabel('B1x'); ylabel('B1y');
subplot(515); plot(A,w1); xlabel(''); ylabel('myuB');hold on;
% figure(2); plot(a);
end