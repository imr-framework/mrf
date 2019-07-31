function [gx,gy,x,y] = make_epi_jfn(fov, npix, etl, prefix, issoft, zres)
% function [gx,gy] = makeEPI(fov, npix, etl, prefix, [issoft, zres])
%
% Design EPI readout gradients for ssfpepi psd.
%
% INPUTS:
%  fov   - cm
%  npix  - etl*nl, where etl = echo-train length and nl = # leaves (shots)
%  etl   - echo-train length (must be an odd number)
%  prefix  - filename prefix (any string)
%  issoft  - (optional) back off on gradient amp/slew 
%  zres    - spatial resolution along z (mm)
%
% example:
%  etl = 17;
%  numshots = 3;
%  fov = 24;
%  makeEPI(fov,etl*numshots,etl,'myepireadout');
%
% $Id: makeEPI.m,v 1.13 2010-07-30 02:25:46 jfnielse Exp $

if nargin < 5
	issoft = false;
end

dt = 4e-3;                           % gradient/daq sample duration (msec)
if issoft
	mxg = 0.9*4.0;
	mxs = 0.9*15.0;                 % make gradients quieter
else
	mxg = 4.0;
	mxs = 15.0;
end
s = mxs * dt;                        % max change in g per sample (G/cm)

% gradient limits for phase/partition-encode blips. 
if nargin < 6
	zres = 10*fov/npix;
end


% ===================== make the various gradient waveform elements ===================== %

% calculate readout amplitude
gamma = 4.2575;                     % kHz/Gauss
g = (1/dt)/(gamma*fov);             % Gauss/cm
if g > 4.0
	g = 4.0;
	fprintf(1, 'max g reduced to %.2f \n', g);
end

% readout trapezoid
gxro = g*ones(1,npix);                         % plateau of readout trapezoid
areapd = sum(gxro)*4e-6;                       % area of partial dephaser gradient (to be used below)
ramp = s/sqrt(2):s/sqrt(2):(g-s/sqrt(2)) ;     % both x and y gradients are playing out during EPI direction switching
gxro = [ramp g*ones(1,npix) fliplr(ramp)];

% x prewinder. make sure res_kpre is even. 
area = sum(gxro)*4e-6;
gxprew = -trapwave(area/2, 4e-6, mxg/sqrt(2), mxs/sqrt(2)*1000);    % 1 mm resolution phase-encode blip
gxprew = [zeros(1,mod(length(gxprew)+length(ramp),2)) gxprew];

% partial dephaser (one cycle of phase across each voxel)
gxpd = -trapwave(areapd/2, 4e-6, mxg/sqrt(3), mxs/sqrt(3)*1000);       % y and z phase-encode blips playing at the same time 
gxpd = [zeros(1,mod(length(gxpd),2)) gxpd];

% phase-encode trapezoids before/after gx
gyprew = -trapwave(areapd/2*(10*fov/npix/zres), 4e-6, mxg/sqrt(3), mxs/sqrt(3)*1000);
gyprew = [zeros(1,mod(length(gyprew),2)) gyprew];

% gy readout gradient elements
areagyblip = areapd/etl;
gyblip = trapwave(areagyblip, 4e-6, mxg/sqrt(2), mxs/sqrt(2)*1000);    % tiny blip to move to next phase-encode line
gyro = [zeros(1,length(gxro)-length(gyblip)) gyblip];



 
% ======================= put together gx and gy ======================= %

gx = gxprew;
gy = zeros(1,length(gxprew)+length(gyblip)/2);
for e = 1:(etl-1)
	gx = [gx (-1)^(e+1)*gxro];
	gy = [gy gyro];
end
if etl == 1
	e = 1;
else
	e = e+1;
end
gx = [gx (-1)^(e+1)*gxro ]; 
gy = [gy zeros(1,length(gx)-length(gy))];  % make sure gy remains zero past end of DAQ

% add rephasers at end of gx and gy readout
% gx = [gx gxprew]; 
areagx = sum(gx)*4e-6;
gxrep = trapwave(-areagx, 4e-6, mxg/sqrt(2), mxs/sqrt(2)*1000);
gx = [gx gxrep];

areagy = sum(gy)*4e-6;   % units = G/cm*sec
gyrep = trapwave(-areagy, 4e-6, mxg/sqrt(2), mxs/sqrt(2)*1000);
gy = [gy gyrep];

% make sure length of gx and gy are same, and even
lx = numel(gx);
ly = numel(gy);
if lx > ly
	gy = [gy zeros(1,lx-ly)];
else
	gx = [gx zeros(1,ly-lx)];
end

gx = [gx zeros(1,mod(length(gx),2))];
gy = [gy zeros(1,mod(length(gy),2))];


figure; plot(gx, 'b');
hold on; plot(gy, 'r'); 

x = cumsum(gx);
y = cumsum(gy);

% figure, plot(x,y);

% ====================== write gradients to .jfnwav files ==================== %

prefix = [prefix '_fov' num2str(fov) '_res' num2str(npix) '_etl' num2str(etl) '_'];
if issoft
	prefix = [prefix 'soft_'];
end

% gx and gy
im = find(gx==g | gx==-g);
hdr.res_g = length(gx);
hdr.res_kpre = im(1)-1;
hdr.res_k = im(end) - im(1) + 1;
hdr.res_k = hdr.res_k + mod(hdr.res_k,2);
hdr.npix = npix;  
hdr.nechoes = etl;  
hdr.nechosep = round(length(gxro)/2);  
% hdr
writejfnwav(gx, [prefix 'gx.jfnwav'], hdr, 1.0, 4.0)
writejfnwav(gy, [prefix 'gy.jfnwav'], hdr, 1.0, 4.0)

% gxpd
hdr.res_g = length(gxpd);
hdr.res_kpre = 0;
hdr.res_k = 0;
hdr.npix = 0;  
hdr.nechoes = 0;  
hdr.nechosep = 0;
writejfnwav(gxpd, [prefix 'gxpd.jfnwav'], hdr, 1.0, 4.0)

% gyprew
hdr.res_g = length(gyprew);
hdr.res_kpre = 0;
hdr.res_k = 0;
hdr.npix = 0;  
hdr.nechoes = 0;  
hdr.nechosep = 0;
writejfnwav(gyprew, [prefix 'gyblip.jfnwav'], hdr, 1.0, 4.0)

% done
 
return;

% EOF
