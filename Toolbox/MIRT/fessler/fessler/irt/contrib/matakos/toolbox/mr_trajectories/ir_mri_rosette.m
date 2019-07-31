 function [kspace, omega, grad, slew, ts, idx, wi, rlen] = ...
	ir_mri_rosette(varargin)
%function [kspace, omega, grad, slew, ts, idx, wi, rlen] = ...
%|	ir_mri_rosette(varargin)
%|
%| Design rosette readout gradients based on GE 3T scanner constraints.
%|
%| in
%|
%| options
%|	fov		Field of view in cm (square)
%|	N		size of reconstructed image (square)
%|	dt		sampling interval in sec
%|	Nt		# of time points for each shot
%|	Tr		total readout time for each shot
%|	shot		# of shots
%|	gamp		gradient amplitude G/cm
%|	gslew		gradient slew rate G/cm/ms
%|	out			file output (flag)
%|	prefix		filename prefix (any string)
%|	rot			trajectory rotation
%|	chat		show info during run
%|
%| out
%|	kspace		kspace trajectory [kx ky] in cycles/cm
%|	omega		kspace trajectory in radians
%|	grad		gradient waveforms [gx gy]
%|	slew		gradient slew rate waveforms [sx sy]
%|	ts		vector of sampling times
%|	idx		index vector of sampled points
%|
%| example:
%|	ir_mri_rosette('fov', 24, 'N', 64, 'dt', 4e-6);
%|
%| Based on rosette3 from mri_trajectory
%| Copyright 2010-10-28, Antonis Matakos, University of Michigan
%| 2014-12-12, tweaks by JF

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(varargin{1}, 'test'), ir_mri_rosette_test; return; end

% defaults for rosette trajectory
arg.fov = 24;		% cm
arg.N = 64;
arg.dt = 4e-6;		% sec
arg.Nt = [];		% # of time sampling points per shot
arg.Tr = [];		% readout time
arg.shot = 1;		% # of shots
arg.gamp = 4;		% gradient amplitude (2014-12-10, doug noll say 5 g/cm)
arg.gslew = 15;		% gradient slew (2014-12-10, doug noll says 20 g/cm/ms)
arg.out = 0;
arg.prefix = 'test_rosette';
arg.rot = 0;		% trajectory rotation
arg.chat = 0;

arg = vararg_pair(arg, varargin);

fov = arg.fov;
if length(fov) > 1
	warn('FOV square, keeping only first value')
	fov = fov(1);
end

N = arg.N;
if length(N) > 1
	warn('square matrix, keeping only first value')
	N = N(1);
end

if isempty(arg.Nt)
	arg.Nt = 2*N^2;		% # of time sampling points per shot
end

dt = arg.dt;

%% Specify rosette parameters (freq, samples per sweep, etc...)
gambar = 4257.5; % gyromagnetic ratio in Hz/G

fmax1 = gambar*arg.gamp*fov/N/pi;
fmax2 = sqrt(gambar*1000*arg.gslew*fov/N/pi^2/2.2);

fmax = min([fmax1 fmax2]);

% Find # samples per sweep and make it a power of 2 if possible
ns = ceil(1/fmax/dt);
% keep 2 largest powers of 2
ns = de2bi(ns,'left-msb');
ns(3:end) = 0;
ns = bi2de(ns,'left-msb') + bi2de(eye(1,length(ns)-1),'left-msb');

% fast frequency of rosette
f1 = 1/ns/dt;

if ~isempty(arg.Tr)
	arg.Nt = floor(arg.Tr/dt);
end

% Find # of sweeps
nw = ceil(arg.Nt/ns);
nw = nw + mod(nw,2);		% make nw even

% update # of samples
arg.Nt = ns*nw;

% find m relatively prime to nw s.t. m/nw approx 0.25
m = floor(nw/4);
for ii=0:m-1
	if isrelprime(m+ii,nw)
		m = m + ii;
		break;
	elseif isrelprime(m-ii,nw)
		m = m - ii;
		break;
	end
end

% find slow frequency of rosette at approx 25% of f1
f2 = f1*m/nw;


%% Create rosette trajectory
kmax = N/fov/2;
ti = (0:arg.Nt-1)' * dt;


tmp = 2 * pi * ti;
p1 = f1 * tmp;
p2 = f2 * tmp;

kx = kmax * sin(p1) .* cos(p2);
ky = kmax * sin(p1) .* sin(p2);

Gx = diff(kx)/dt/gambar;
Gy = diff(ky)/dt/gambar;

wk = kmax^2*abs(sin(p1).*cos(p1)); % Analytical weights as in bucholz:08:miw
wi = repmat(wk,[1 arg.shot]);

% Fix gradients to start from 0 from both ends
areax = sum(Gx)*dt;
gxrew = trapwave(-areax, dt, arg.gamp, arg.gslew*1000,Gx(end));
Gx = [Gx;gxrew.'];

areay = sum(Gy)*dt;
gyrew = trapwave(-areay, dt, arg.gamp, arg.gslew*1000,Gy(end));
Gy = [Gy;gyrew.'];

lx = numel(Gx);
ly = numel(Gy);
if lx > ly
	Gy = [Gy;zeros(lx-ly,1)];
	len1= length(gxrew);
else
	Gx = [Gx;zeros(ly-lx,1)];
	len1 = length(gyrew);
end

areax = sum(Gx)*dt;
gxprew = fliplr(trapwave(-areax, dt, arg.gamp, arg.gslew*1000,Gx(1)));
Gx = [gxprew.';Gx];

areay = sum(Gy)*dt;
gyprew = fliplr(trapwave(-areay, dt, arg.gamp, arg.gslew*1000,Gy(1)));
Gy = [gyprew.';Gy];

lx = numel(Gx);
ly = numel(Gy);
if lx > ly
	Gy = [zeros(lx-ly,1);Gy];
	len2 = length(gxprew);
else
	Gx = [zeros(ly-lx,1);Gx];
	len2 = length(gyprew);
end

kxo = cumsum(Gx)*gambar*dt;
kyo = cumsum(Gy)*gambar*dt;

sxo = [0 ;diff(Gx)/dt/1000];
syo = [0 ;diff(Gy)/dt/1000];

% Apply rotations
kx = zeros(length(kxo),arg.shot);
ky = zeros(length(kyo),arg.shot);

gx = zeros(length(Gx),arg.shot);
gy = zeros(length(Gy),arg.shot);

sx = zeros(length(sxo),arg.shot);
sy = zeros(length(syo),arg.shot);

idx = [false(1,len2-1) true(1,arg.Nt) false(1,len1)];
idx = repmat(idx,[1 arg.shot]);


kx(:,1) = kxo;
ky(:,1) = kyo;

gx(:,1) = Gx;
gy(:,1) = Gy;

sx(:,1) = sxo;
sy(:,1) = syo;

for is=1:(arg.shot-1) % n-shot, rotate kx,ky by 2 pi / N
	ang = is * pi / arg.shot / nw;
	c = cos(ang);
	s = sin(ang);
	kx(:,is+1) = c * kxo + s * kyo;
	ky(:,is+1) = -s * kxo + c * kyo;

	gx(:,is+1) = c * Gx + s * Gy;
	gy(:,is+1) = -s * Gx + c * Gy;

	sx(:,is+1) = c * sxo + s * syo;
	sy(:,is+1) = -s * sxo + c * syo;

end

kspace = [kx(:) ky(:)];
omega = 2*pi*fov*kspace / N;
grad = [gx(:) gy(:)];
slew = [sx(:) sy(:)];
ts = (0:length(kspace)-1)'*dt;
idx = find(idx(:)>0);
wi = wi(:)*4*pi/arg.Nt/arg.shot;

rlen = idx(end) - idx(1) + 1;
rlen = rlen + mod(rlen,2);

if arg.chat
	figure, plot(kspace(:,1), kspace(:,2), '.-'), axis square
	figure, plot(omega(:,1), omega(:,2), '.-'), axis square
	figure, plot(ts,grad(:,1), 'b.-',ts,grad(:,2), 'r.-')
	figure, plot(ts,slew(:,1), 'b.-',ts,slew(:,2), 'r.-')
end


% ====================== write gradients to .jfnwav files ==================== %
if arg.out

	arg.prefix = [arg.prefix '_sh' num2str(arg.shot) ...
		'_fov' num2str(fov) '_res' num2str(N) '_'];

	% gx and gy
	hdr.res_g = length(gx);
	hdr.npix = N;
	hdr.res_kpre = idx(1) - 1;
	hdr.res_k = rlen;
	hdr.nechoes = N;
	hdr.nechosep = 0;

	% hdr
	writejfnwav(gx, [arg.prefix 'gx.jfnwav'], hdr, 1.0, 4.0)
	writejfnwav(gy, [arg.prefix 'gy.jfnwav'], hdr, 1.0, 4.0)

	% gxpd
	gxpd = zeros(1,100);
	hdr.res_g = length(gxpd);
	hdr.res_kpre = 0;
	hdr.res_k = 0;
	hdr.N = 0;
	hdr.nechoes = 0;
	hdr.nechosep = 0;
	writejfnwav(gxpd, [arg.prefix 'gxpd.jfnwav'], hdr, 1.0, 4.0)

	% gyprew
	gyprew = zeros(1,100);
	hdr.res_g = length(gyprew);
	hdr.res_kpre = 0;
	hdr.res_k = 0;
	hdr.N = 0;
	hdr.nechoes = 0;
	hdr.nechosep = 0;
	writejfnwav(gyprew, [arg.prefix 'gyblip.jfnwav'], hdr, 1.0, 4.0)

end

end


function res = isrelprime(m,n)

if m > n
	tmp = m;
	m = n;
	n = tmp;
end


while m > 0
	mold = m;
	m = rem(n,m);
	n = mold;
end

if mold == 1
	res = true;
else
	res = false;
end

end


% ir_mri_rosette_test
function ir_mri_rosette_test

[ks om gr sl ts idx] = ir_mri_rosette(...
	'fov', 24, 'N', 64, 'dt', 4e-6, 'Nt', 8192, 'chat', 0, 'shot', 1);

ts = 1000*ts;

if im
	pl = @(i) subplot(230+i); clf
	pl(1); plot(ts, gr(:,1), 'b.-', ts, gr(:,2), 'r.-')
	pl(2); plot(ts, sl(:,1), 'b.-', ts, sl(:,2), 'r.-')
	pl(3); plot(ks(:,1), ks(:,2), '.-'), axis square
	pl(4); plot(om(:,1), om(:,2), '.-'), axis square
end

if ~isempty(idx)
	ks = ks(idx,:);
	gr = gr(idx,:);
	ts = ts(idx);
	if im
		pl(5), plot(ks(:,1), ks(:,2), '.-'), axis square
		pl(6), plot(ts, gr(:,1), 'b.-', ts, gr(:,2), 'r.-')
	end
end

end % ir_mri_rosette_test
