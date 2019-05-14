function [Gx, Gy, Gz] = k2g(k,dt)
% k = k.*1e3; %conversion to m-1
% gammabar = 42.57e6;%Hz/T
% fact = 1/gammabar;

fact =1;% for units in Hz/m; gammabar gets cancelled
Gx = (fact./dt).*(diff(squeeze(k(1,:))));
Gy = (fact./dt).*(diff(squeeze(k(2,:))));
if(size(k,1) >2)
    Gz =  (fact./dt).*(diff(squeeze(k(3,:))));
else
    Gz =0;
end



Gx = [0 Gx];
Gy = [0 Gy];
Gz = [0 Gz];
%Convert to Hz/m