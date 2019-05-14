
%% GenerateSpiralTraj.m
%
%[k]=GenerateSpiralTraj_v1()
% Generates a spiral k-space trajectory with a method adapted from [1]
% (the paper is a bit buggy).
%
% INPUTS:   * field of view in meters
%           * resolution (Nyquist distance) in meters
%           * undersampling factor along frequency encoding direction (normalized units, positive, may be lower than one)
%           * undersampling factor (normalized units, positive, may be lower than one)
%           * number of interleaves
%           * variable density factor (alpha in [1])
%           * maximum amplitude (in T/m)
%           * maximum slew rate (in T/m/s)
%           * resampling trajectory (true or false)
%           * analysis option (true or false)
%
% [1]   "Simple Analytic Variable Density Spiral Design"
%       Dong-hyun Kim, Elfar Adalsteinsson, and Daniel M. Spielman
%       Magnetic Resonance in Medicine 50:214-219 (2003)
%
% Copyright, Matthieu Guerquin-Kern, 2012
% Modified by : Pavan Poojar
% Modified by: Sairam Geethanath, MIRC & CMRRP - integrated with Pulseq -
% still flaky - heavily dependent on choice of alpha
 
function [kshot,Gn,lambda]= vds2D_pulseq_v1(FOV,N,Nshots,alpha,lims)
gamma =  42576000; % in Hz  %Determined from Pulseq - do not change
% Gmax=33e-3;%T/m
% SRmax=170 ; %T/m/s - 100+
Gmax=lims.maxGrad/gamma;             %T/m
SRsafety = 0.6;
SRmax=SRsafety*lims.maxSlew/gamma;             %T/m/s - safety measures
res=FOV/N;
%% Generating first interleave

lambda = .5/res; % in m^(-1)
n = (1/(1-(1-Nshots/FOV/lambda)^(1/alpha)));
w = 2*pi*n;
Tea = lambda*w/gamma/Gmax/(alpha+1); % in s
%Tea1 = (lambda*w)/((gamma.*gm)*(alpha+1));
Tes = sqrt(lambda*w^2/SRmax/gamma)/(alpha/2+1); % in s
Ts2a = (Tes^((alpha+1)/(alpha/2+1))*(alpha/2+1)/Tea/(alpha+1))^(1+2/alpha); % in s

% Ts2a = Ts2a - mod(Ts2a,1e-4);
% Tes = Tes - mod(Tes,1e-4);

if Ts2a<Tes
    tautrans = (Ts2a/Tes).^(1/(alpha/2+1));
    tau = @(t) (t/Tes).^(1/(alpha/2+1)).*(0<=t).*(t<=Ts2a)+((t-Ts2a)/Tea + tautrans^(alpha+1)).^(1/(alpha+1)).*(t>Ts2a).*(t<=Tea).*(Tes>=Ts2a);
    Tend = Tea;
else
  tau = @(t) (t/Tes).^(1/(alpha/2+1)).*(0<=t).*(t<=Tes);
  Tend = Tes;
end
 
k = @(t) lambda*tau(t).^alpha.*exp(1i*w*tau(t));
dt = Tea*1e-4; % in s
Dt = dt/FOV/abs(k(Tea)-k(Tea-dt)); % in s
t = 0:Dt:Tend; % in s - Tend
kt = k(t); % in m^-1

%% Resample at gradRaster
DT=lims.gradRasterTime;
tn=0:DT:Tend;
 kn= k(tn);
 
 %%Obtain gradient waveforms
 ktn(1,:) = real(kn);
 ktn(2,:) = imag(kn);
 
[Gx, Gy, ~] = k2g(ktn,DT);

% figure(103); plot(Gx/gamma); hold on; plot(Gy/gamma);

%%  Should start from here again
figure(1030);
  if(mod(length(kn),2)~=0)
                kn = [0  kn];
  end
      kshot = zeros(length(kn), Nshots);
      Gn = zeros(length(kn), Nshots);



            for s=1:Nshots
                kshot(:,s) =kn.*exp(2*pi*1i*s/Nshots);%m-1
%                 Gn(:,s) =[0  (diff(squeeze(kn(:,s))).*1e3)/dt]; %mT/m
                km(1,:) = real(squeeze(kshot(:,s)));
                km(2,:) = imag(squeeze(kshot(:,s)));
                [Gx, Gy] = k2g(km,DT);
                 if(mod(length(Gx),2)~=0)
                    Gx = [0 Gx];
                    Gy = [0; Gy];
                end
%                 Gx =[0 Gx];Gy =[0 Gy];
                Gn(:,s) = complex(Gx,Gy);%Hz/m
                %% Check for SR violations here
                gx = mr.makeArbitraryGrad('x', squeeze(Gx),lims);
                gy = mr.makeArbitraryGrad('x', squeeze(Gy),lims);
                plot(squeeze(kshot(:,s)));hold on;grid on;xlabel('kx (m-1)'); ylabel('ky(m-1)');
                disp(s);
            end


%         if(( max(abs(Gn(:))) - mod(max(abs(Gn(:))),0.1)) > (Gmax*1e3))
%         %mT/ 
%             error('Design error')
%         end

end
                
 
 
 
 
 


