function [kspace,omega,grad,slew,ts,idx,rlen] = make_epi(fov, N, etl, dt, varargin)
%| function [gx,gy] = make_epi(fov, npix, etl, dt, options)
%|
%| Design EPI readout gradients for ssfpepi psd.
%|
%| in:
%|	fov		Field of view in cm
%|	N		# of pixels (square). N = etl*nl, where 
%|			etl = echo-train length and nl = # leaves (shots)
%|	etl		echo-train length
%|	dt		sampling rate in sec
%|
%| options:
%|	type	type of EPI trajectory
%|		'std'	standard EPI
%|		'usmp'	standard undersampled EPI (single-shot)
%|		'int'	interleaved EPI w/ 2 interleaves going bottom up twice
%|		'intud'	interleaved EPI going bottom up and then top bottom 
%|		'mod'	modified interleaved EPI to avoid N/2 ghost
%|	offset	used for multi-shot EPI goes from 0 to #shots-1
%|	dirx	x direction of EPI -1 left to right and 1 right to left
%|	diry	y direction of EPI -1 bottom-up 1 top-bottom
%|	smpr	sampling on ramps (flag)
%|	gamp	gradient amplitude
%|	glsew	gradient slew rate
%|	out		file output (flag)
%|	prefix	filename prefix (any string)
%|	issoft	(optional) back off on gradient amp/slew 
%|	zres	spatial resolution along z (mm)
%|	chat	show info during run
%|
%| out:
%|	kspace	kspace trajectory [kx ky] in cycles/cm
%|	omega	kspace trajectory in radians
%|	grad	gradient waveforms [gx gy]
%|	slew	gradient slew rate waveforms [sx sy]
%|	ts		vector of sampling times
%|	idx		index vector of sampled points
%|
%| example:
%|	etl = 32;
%|	fov = 24;
%|	N = 64;
%|	dt = 4e-6;
%|	make_epi(fov,N,etl,dt);
%|
%| Based on files from J. F. Nielsen
%| Copyright 2010-08-07, Antonis Matakos, The University of Michigan

if nargin < 1
	make_epi_test;
	return;
end

if ~isvar('fov')
	fov = 24;
end

if ~isvar('N')
	N = 64;
end

if ~isvar('etl')
	etl = N;
end

if ~isvar('dt')
	dt = 4e-6;
end

arg.type = 'int';
arg.offset = 0;
arg.dirx = -1;
arg.diry = -1;
arg.smpr = 0;
arg.gamp = 4;
arg.gslew = 15;
arg.out = 0;
arg.prefix = 'test_epi';
arg.issoft = 0;
arg.zres = 10*fov/N;
arg.chat = 0;

arg = vararg_pair(arg, varargin);

nshot = N/etl;
if arg.offset < 0 || arg.offset >= nshot
	warning('make_epi:OffsetError', 'Incorrect value for offset, reset to 0');
	arg.offset = 0;
end


if arg.issoft
	mxg = 0.9*arg.gamp;
	mxs = 0.9*arg.gslew;	% make gradients quieter
else
	mxg = arg.gamp;
	mxs = arg.gslew;
end
s = mxs * dt * 1000;

% Why is scaling necessary?
% scale = sqrt(2);
scale = 1;

% scale2 = sqrt(3);
scale2 = 1;

scaley = 20;

% =============== make the various gradient waveform elements =============== %

% calculate readout amplitude
gamma = 4.2575;					% kHz/Gauss
g = (1/(1000*dt))/(gamma*fov);	% Gauss/cm
if g > mxg
	g = mxg;
	if arg.chat
		fprintf(1, 'max g reduced to %.2f \n', g);
	end
end

% readout trapezoid
gxro = g*ones(1,N);		% plateau of readout trapezoid
areapd = sum(gxro)*dt;	% area of partial dephaser gradient (to be used below)

% Added sampling on ramps option - recon using NUFFT
if arg.smpr
	gxro = trapwave(areapd, dt, mxg/scale, mxs/scale*1000);
else
	ramp = s/scale:s/scale:g; % change here to restrict slew rate
	gxro = [0 ramp g*ones(1,N) fliplr(ramp)]; % added 0 to restrict slew rate
end

% x prewinder. make sure res_kpre is even. 
% Handle even N by changing prewinder
if mod(N,2) == 0
	area = (sum(gxro)-arg.dirx*g)*dt;
else
	area = sum(gxro)*dt;
end
gxprew2 = arg.dirx*trapwave(area/2 + arg.dirx*g*dt, dt, mxg/scale, mxs/scale*1000);
gxprew = arg.dirx*trapwave(area/2, dt, mxg/scale, mxs/scale*1000);
if arg.smpr
	gxprew = [zeros(1,mod(length(gxprew),2)) gxprew];
else
	gxprew = [zeros(1,mod(length(gxprew)+length(ramp),2)) gxprew];
end

% partial dephaser (one cycle of phase across each voxel)
gxpd = -trapwave(areapd/2, dt, mxg/scale2, mxs/scale2*1000);
gxpd = [zeros(1,mod(length(gxpd),2)) gxpd];

% phase-encode trapezoids before/after gx
% gyprew = -trapwave(areapd/2*(10*fov/N/arg.zres), dt, mxg/scale2, mxs/scale2*1000);
% Handle even N by changing prewinder
if mod(N,2) == 0
	areayprew = areapd/2 - arg.offset*g*dt;
else
	areayprew = (areapd - g*dt)/2 - arg.offset*g*dt;
end
if strcmpi(arg.type, 'int')
	areayprew = areayprew - 3*g*dt;
	gyprew = -arg.diry*trapwave(areayprew, dt, mxg/scale2, mxs/scaley*1000);
else
	gyprew = arg.diry*trapwave(areayprew, dt, mxg/scale2, mxs/scaley*1000);
end
gyprew = [zeros(1,mod(length(gyprew),2)) gyprew];

lx = numel(gxpd);
ly = numel(gyprew);
if lx > ly
	gyprew = [gyprew zeros(1,lx-ly)];
else
	gxpd = [gxpd zeros(1,ly-lx)];
end

% gy readout gradient elements
% Changed readout patterns to create interleaved EPIs
areagyblip = areapd/etl;
if strcmpi(arg.type, 'std') || strcmpi(arg.type, 'usmp')
	gyblip = trapwave(areagyblip, dt, mxg/scale, mxs/scaley*1000);
	gyro = [zeros(1,length(gxro)-length(gyblip)) gyblip];
	gyro2 = 0;
elseif strcmpi(arg.type, 'int')
	gyblip = trapwave(2*areagyblip, dt, mxg/scale, mxs/scaley*1000);
	gyro = [zeros(1,length(gxro)-length(gyblip)) gyblip];

	armult = 2*ceil((etl-2)/2) - 1;
	gyblip2 = trapwave(armult*areagyblip, dt, mxg/scale, mxs/scaley*1000);
	gyro2 = [zeros(1,length(gxro)-length(gyblip2)) gyblip2];
elseif strcmpi(arg.type, 'mod')
	gyblip = trapwave(areagyblip, dt, mxg/scale, mxs/scaley*1000);
	gyro = [zeros(1,length(gxro)-length(gyblip)) gyblip];
	
	gyblip1 = trapwave(3*areagyblip, dt, mxg/scale, mxs/scaley*1000);
	gyro1 = [zeros(1,length(gxro)-length(gyblip1)) gyblip1];

	armult = 2*ceil((etl-2)/2) - 3;
	gyblip2 = trapwave(armult*areagyblip, dt, mxg/scale, mxs/scaley*1000);
	gyro2 = [zeros(1,length(gxro)-length(gyblip2)) gyblip2];
elseif strcmpi(arg.type, 'intud')
	gyblip = trapwave(2*areagyblip, dt, mxg/scale, mxs/scaley*1000);
	gyro = [zeros(1,length(gxro)-length(gyblip)) gyblip];

	gyblip2 = trapwave(areagyblip, dt, mxg/scale, mxs/scaley*1000);
	gyro2 = [zeros(1,length(gxro)-length(gyblip2)) gyblip2];
else
	error 'Invalid trajectory type'
end


 
% ======================= put together gx and gy ======================= %

gxro = -arg.dirx*gxro;
gx = gxprew;

gyro = -arg.diry*gyro;
gyro2 = -arg.diry*gyro2;
% gy = gyprew;
gy = 0;

lx = numel(gx);
ly = numel(gy);
if lx > ly
	gy = [gy zeros(1,lx-ly)];
else
	gx = [gx zeros(1,ly-lx)];
end

im(1) = length(gx)+1;

gy = [gy zeros(1,length(gyblip)/2)];

if strcmpi(arg.type, 'std') || strcmpi(arg.type, 'usmp')
	for e = 1:(etl-1)
		gx = [gx (-1)^(e+1)*gxro];
		gy = [gy gyro];
	end
elseif strcmpi(arg.type, 'int')
	gx = [-gxprew2 -gxro 0*gyro2];
	gy = [gy 0*gyro -gyro2];
	for e = 1:(etl-1)
		if e == floor((etl+1)/2)
			gx = [gx (-1)^(e+1)*gxro 0*gyro2];
			gy = [gy 0*gyro -gyro2];
		else
			gx = [gx (-1)^(e+1)*gxro];
			gy = [gy gyro];
		end
	end
elseif strcmpi(arg.type, 'mod')
	for e = 1:(etl-1)
		if e == floor((etl+1)/2)
			gx = [gx (-1)^(e+1)*gxro 0*gyro2];
			gy = [gy 0*gyro -gyro2];
		elseif mod(e,2)
			gx = [gx (-1)^(e+1)*gxro];
			gy = [gy gyro];
		else
			gx = [gx (-1)^(e+1)*gxro];
			gy = [gy gyro1];
		end
	end
elseif strcmpi(arg.type, 'intud')
	for e = 1:floor((etl-1)/2)
		gx = [gx (-1)^(e+1)*gxro];
		gy = [gy gyro];
	end

	e = floor((etl+1)/2);
	gx = [gx (-1)^(e+1)*gxro];

	if mod(etl,2) == 0
		gy = [gy gyro2];
	else
		gy = [gy -gyro2];
	end

	for e = floor((etl+3)/2):(etl-1)
		gx = [gx (-1)^(e+1)*gxro];
		gy = [gy -gyro];
	end
else
	error 'Invalid trajectory type'
end

if etl == 1
	e = 1;
else
	e = e+1;
end
gx = [gx (-1)^(e+1)*gxro 0]; % added 0 to restrict slew rate
gy = [gy zeros(1,length(gx)-length(gy))];  % make sure gy remains zero past end of DAQ

im(2) = length(gx);

% add rephasers at end of gx and gy readout
areagx = sum(gx)*dt;
gxrep = trapwave(-areagx, dt, mxg/scale, mxs/scale*1000);
gx = [gx gxrep];

areagy = sum(gy)*dt;   % units = G/cm*sec
gyrep = trapwave(-areagy, dt, mxg/scale, mxs/scaley*1000);
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

sx = diff(gx)./(dt*1000);
sy = diff(gy)./(dt*1000);

% Use interpolation to create accurate sampling grid avoiding offsets
% caused by the way cumsum works
Gx = [zeros(size(gyprew)) ...
	interp1(1:length(gx),gx,[1 (2:length(gx))-.5]) zeros(size(gyprew))];
% Gy = interp1(1:length(gy),gy,[1 (2:length(gy))-.5]);
Gy = [gyprew gy -gyprew];

% kx = cumsum(gx)*gamma*dt;
% ky = cumsum(gy)*gamma*dt;

% change here!
kx = cumsum(Gx)*gamma*dt*1000;
ky = cumsum(Gy)*gamma*dt*1000;


kspace = [kx.' ky.'];
omega = 2*pi*fov*kspace / N;
grad = [Gx.' Gy.'];
slew = [sx.' sy.'];

if arg.chat
	figure, plot(kx, ky, '.');
	figure, plot(1:length(gx), gx, 'b', 1:length(gy), gy, 'r');
	figure, plot(1:length(sx), sx, 'b', 1:length(sy), sy, 'r');
end

ts = (0:length(kspace)-1)*dt;
ts = ts';

if arg.smpr
	idx = im(1):im(2);
else
	idx = find(gx==g | gx==-g); % bug here!!! FIXED!
end
idx = idx';

if strcmpi(arg.type, 'int')
	idx = idx(N+1:end);
end

rlen = idx(end) - idx(1) + 1;
rlen = rlen + mod(rlen,2);

% ====================== write gradients to .jfnwav files ==================== %
if arg.out
	
	arg.prefix = [arg.prefix '_' arg.type '_fov' num2str(fov) ...
		'_res' num2str(N) '_etl' num2str(etl) '_'];
	if arg.issoft
		arg.prefix = [arg.prefix 'soft_'];
	end

	% gx and gy
	hdr.res_g = length(gx);
	hdr.npix = N;
	hdr.res_kpre = idx(1)-1;
	hdr.res_k = rlen;
	if strcmpi(arg.type, 'usmp')
		hdr.nechoes = N;
	else
		hdr.nechoes = etl;
	end
	hdr.nechosep = round(length(gxro)/2);

	% hdr
	writejfnwav(gx, [arg.prefix 'gx.jfnwav'], hdr, 1.0, 4.0)
	writejfnwav(gy, [arg.prefix 'gy.jfnwav'], hdr, 1.0, 4.0)

	% gxpd
	hdr.res_g = length(gxpd);
	hdr.res_kpre = 0;
	hdr.res_k = 0;
	hdr.npix = 0;
	hdr.nechoes = 0;
	hdr.nechosep = 0;
	writejfnwav(gxpd, [arg.prefix 'gxpd.jfnwav'], hdr, 1.0, 4.0)

	% gyprew
% 	gyprew = [0 0];
	hdr.res_g = length(gyprew);
	hdr.res_kpre = 0;
	hdr.res_k = 0;
	hdr.npix = 0;
	hdr.nechoes = 0;
	hdr.nechosep = 0;
	writejfnwav(gyprew, [arg.prefix 'gyblip.jfnwav'], hdr, 1.0, 4.0)

	% done
end

idx = idx + length(gyprew);

end

% EOF

function make_epi_test

fov = 24;
N = 64;
etl = 64;
dt = 4e-6;

[ks, om, gr, sl, ts, idx] = make_epi(fov,N,etl,dt,'type','int','smpr',0,...
	'dirx',-1,'diry',1,'chat',0);
% [ks2, om2, gr2, sl2, ts2, idx2] = make_epi(fov,N,etl,dt,...
% 	'offset',1,'type','int','smpr',0,'dirx',1);

% tst = ts(1:end-1);
% tst2 = ts2(1:end-1);

% ks = [ks;ks2];
% om = [om;om2];
% gr = [gr;gr2];
% sl = [sl;sl2];
% ts = [ts;ts2];
% tst = [tst;tst2];
% idx = [idx;idx2+length(ks2)];

ts = 1000*ts;

figure, plot(1:length(gr), gr(:,1), '.b', 1:length(gr), gr(:,2), '.r');
figure, plot(1:length(sl), sl(:,1), 'b', 1:length(sl), sl(:,2), 'r');
figure, plot(ks(:,1), ks(:,2), '.'), axis square;
figure, plot(om(:,1), om(:,2), '.'), axis square;

if ~isempty(idx)
	ks = ks(idx,:);
	gr = gr(idx-72,:);
	ts = ts(idx);
	figure, plot(ks(:,1), ks(:,2), 'b.-'), axis square;
	figure, plot(ts, gr(:,1), '.b', ts, gr(:,2), '.r');
end

end

