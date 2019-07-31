function [kspace, omega, grad, slew, ts, idx, rlen] = ...
	make_spiral(fov,N,dt,varargin)
%| function [kspace, omega] = mri_kspace_spiral(fov,N,dt,options)
%| k-space spiral trajectory based on GE 3T scanner constraints
%|
%| in
%|	fov		field of view in cm
%|	N		size of reconstructed image (square)
%|	dt		time sampling interval
%|
%| options (name / value pairs)
%|	type	type of spiral trajectory
%|		'std'	standard spiral
%|		'int'	interleaved spiral w/ 2 interleaves going in twice
%|		'intio'	interleaved spiral w/ 2 interleaves going in and out
%|	Nt		# of time points
%|	shot	# of shots
%|	gamp	gradient amplitude
%|	gslew	gradient slew rate
%|	out		file output (flag)
%|	prefix	filename prefix (any string)
%|	rot		trajectory rotation
%|	rev		trajectory reversal (flag) - used for spiral in
%|	chat	show info during run
%|
%| out
%|	kspace [Nt,2]	kspace trajectory [kx ky] in cycles/FOV, NO:cycles/cm 
%|	omega [Nt,2]	kspace trajectory [kx ky] in radians
%|	grad [Nt ,2]	gradient waveforms
%|	slew [Nt ,2]	gradient slew rate waveforms
%|	ts		vector of sampling times
%|	idx		index vector of sampled points
%|
%| example:
%|	fov = 24;
%|	N = 64;
%|	dt = 4e-6;
%|	make_rosette(fov,N,dt);
%|
%| based on m-files from valur that he got from brad who got them from doug...
%| Edited 2010-08-10, Antonis Matakos, University of Michigan

if nargin < 1
	make_spiral_test;
	return;
end

if ~isvar('fov')
	fov = 24;
end

if ~isvar('N')
	N = 64;
end

if ~isvar('dt')
	dt = 4e-6;
end

% defaults for spiral trajectory
arg.type = 'int';
arg.Nt = [];		% # of time sampling points
arg.shot = 1;
arg.gamp = 4;		% gradient amplitude
arg.gslew = 140;	% gradient slew
arg.out = 0;
arg.prefix = 'test_spiral';
arg.rot = 0;		% trajectory rotation
arg.rev = 1;		% trajectory reversal
arg.chat = 0;

arg = vararg_pair(arg, varargin);

if strcmpi(arg.type, 'std')
	nl = arg.shot;
else
	nl = 2*arg.shot;
end

% provide sensible default Nt
if isempty(arg.Nt)
	if fov == 20
		arg.Nt = 4026;
	elseif fov == 22
		arg.Nt = 3770;
	else
		arg.Nt = 0; % let algorithm choose
% 		warning('spiral:UnknownFOV','unknown FOV: specify Nt?');
	end
end

if fov > 100
	warning('spiral:FOVunits','fov > 100; use cm not mm!');
end

% generate spiral k-space trajectory
[kx,ky,gx,gy,sx,sy,idx] = genkspace(fov, N, nl*arg.Nt, nl,arg.shot,...
	arg.gamp, arg.gslew, dt, arg.rot, arg.rev, arg.type, arg.chat);

gx = [gx;zeros(1,mod(length(gx),2))];
gy = [gy;zeros(1,mod(length(gy),2))];

kspace = [kx ky]/fov;
omega = 2*pi*[kx ky] / N;
grad = [gx gy];
slew = [sx sy];

if arg.chat
	figure, plot(kx, ky, '.');
	figure, plot(gx, 'b', gy, 'r');
	figure, plot(sx, 'b', sy, 'r');
end

ts = (0:length(kspace)-1)*dt;
ts = ts';
idx = find(idx>0);

rlen = idx(end) - idx(1) + 1;
rlen = rlen + mod(rlen,2);

% if max(omega(:)) > pi
% 	error 'bad spiral'
% end

if arg.out
	
	arg.prefix = [arg.prefix '_sh' num2str(arg.shot) '_' arg.type ...
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
 

%
% genkspace
%
function [kx, ky, gx, gy, sx, sy, idx] = genkspace(FOV, N, ld, nint,nshot,...
		gamp, gslew, tsamp, rotamount, rev_flag, type, chat)
% Generate the proper length of k-space trajectory.
% It linearly interpolates the output of genspiral to the correct length
% and takes care of the rotations for the interleaves.
% ld is the length of the data
% nint is the number of interleaves
%
% Brad Sutton, University of Michigan
% Edited 2010/08/10, Antonis Matakos, University of Michigan

flag = 0;	% auto determine number of k-space points

if ~isvar('rotamount')
  rotamount = 0;
end
if ~isvar('rev_flag')
  rev_flag = 0;
end

nk = ld/nint;
if round(nk) ~= nk
  sprintf('Input should have num data pts/number of interleaves must be int')
end

ns = nint/nshot;
if round(ns) ~= ns
  sprintf('Input should have # interleaves/# of shots must be int')
end

if (nk == 0)
     flag = 1;
end

dt = 4e-6;

[Gx,Gy,kxi,kyi,sxi,syi, idx] = genspi(FOV,N,nint,gamp,gslew);

kxt=interp1(0:dt:dt*length(kxi)-dt,kxi,0:tsamp:dt*length(kxi)-tsamp)';
kyt=interp1(0:dt:dt*length(kyi)-dt,kyi,0:tsamp:dt*length(kyi)-tsamp)';

gxt=interp1(0:dt:dt*length(Gx)-dt,Gx,0:tsamp:dt*length(Gx)-tsamp)';
gyt=interp1(0:dt:dt*length(Gx)-dt,Gy,0:tsamp:dt*length(Gx)-tsamp)';

sxt=interp1(0:dt:dt*length(sxi)-dt,sxi,0:tsamp:dt*length(sxi)-tsamp)';
syt=interp1(0:dt:dt*length(syi)-dt,sxi,0:tsamp:dt*length(syi)-tsamp)';

% idx=interp1(0:dt:dt*length(idx)-dt,idx,0:tsamp:dt*length(idx)-tsamp)';

if flag
	nk = length(kxt)-2;
end

kx = zeros(nk,ns);
ky = zeros(nk,ns);

kxo = kxt(1:nk);
kyo = kyt(1:nk);

gx = zeros(nk,ns);
gy = zeros(nk,ns);

gxo = gxt(1:nk);
gyo = gyt(1:nk);

sx = zeros(nk,ns);
sy = zeros(nk,ns);

sxo = sxt(1:nk);
syo = syt(1:nk);

idx = idx(1:nk);
idx = repmat(idx,[1 ns]);

%rotate matrix for proper orientation
phir = -rotamount*pi/2;
kxop = kxo*cos(phir) - kyo*sin(phir);
kyop = kyo*cos(phir) + kxo*sin(phir);

if rev_flag
	kxop = -flipud(kxop);
	kyop = -flipud(kyop);
	idx = flipud(idx);
end

if ns > 1 && chat
	sprintf('Performing %d rotations', nint)
end
kx(:,1) = kxop;
ky(:,1) = kyop;
phi = 2*pi/ns;
for ii = 1:(ns-1)
	kx(:,ii+1) = kxop*cos(ii*phi) - kyop*sin(ii*phi);
	ky(:,ii+1) = kyop*cos(ii*phi) + kxop*sin(ii*phi);
end

gxop = gxo*cos(phir) - gyo*sin(phir);
gyop = gyo*cos(phir) + gxo*sin(phir);

if rev_flag
	gxop = -flipud(gxop);
	gyop = -flipud(gyop);
end

gx(:,1) = gxop;
gy(:,1) = gyop;
for ii = 1:(ns-1)
	gx(:,ii+1) = gxop*cos(ii*phi) - gyop*sin(ii*phi);
	gy(:,ii+1) = gyop*cos(ii*phi) + gxop*sin(ii*phi);
end

sxop = sxo*cos(phir) - syo*sin(phir);
syop = syo*cos(phir) + sxo*sin(phir);

if rev_flag
	sxop = -flipud(sxop);
	syop = -flipud(syop);
end

sx(:,1) = sxop;
sy(:,1) = syop;
for ii = 1:(ns-1)
	sx(:,ii+1) = sxop*cos(ii*phi) - syop*sin(ii*phi);
	sy(:,ii+1) = syop*cos(ii*phi) + sxop*sin(ii*phi);
end

if strcmpi(type, 'intio')
	kx = [kx(:,1);flipud(kx(:,2))];
	ky = [ky(:,1);flipud(ky(:,2))];

	gx = [gx(:,1);flipud(gx(:,2))];
	gy = [gy(:,1);flipud(gy(:,2))];

	sx = [sx(:,1);flipud(sx(:,2))];
	sy = [sy(:,1);flipud(sy(:,2))];
	
	idx = [idx(:,1);flipud(idx(:,2))];
else
	kx = kx(:);
	ky = ky(:);

	gx = gx(:);
	gy = gy(:);

	sx = sx(:);
	sy = sy(:);
	
	idx = idx(:);
end

end



%
% genspi()
% this is translation of C code from scanner, exactly what is played
% out to gradients at 4us. 
%
function [Gx, Gy, kx, ky, sx, sy, idx] = genspi(D, N, nl, gamp, gslew)
%function [Gx, Gy, kx, ky, sx, sy] = genspi(D, N, nl, gamp, gslew)
%   multi- shot spiral design
%    uses Duyn's approximate slewrate limited design
%    augmented with archimedian gmax limit
% inputs (args)
%        D = FOV, cm
%        N = matrix size
%	 Tmax = longest acquisition allowed, s
%	 dts = output sample spacing, s
%        gtype = trajectory type
% inputs (CVs)
%	nl = number of interleaves
%	gamp = design grad max, G/cm
%	gslew = design slew rate, mT/m/ms
%	nramp = number of rampdown points
% out
%	Gx, Gy
%	grev
% time is in sec
%
%	rev 0 12/26/98	original
%	rev 1 4/15/99	little better calc of ts
%
% borrowed from Doug Noll, Univ. of Michigan
% modified to take more input cv's

%%%%%%%%%% Predefined variables

GRESMAX= 21000;
if ~isvar('nl')
    nl=2;   % Number of interleaves
end
if ~isvar('gamp')
      gamp=4; %3.50; % 2.2 for both 1.5 T and 3 T data
end
if ~isvar('gslew')
      gslew=150; % 200 % 180 for 3T data and 120 (150) for 1.5 T data
end

gts = 4e-6;

Tmax = GRESMAX*gts;

dts = gts;
opfov = D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gambar = 4.257e3;
gamma = 2*pi*gambar;

% gx=zeros(1,2*GRESMAX);
% gy=zeros(1,2*GRESMAX);


q = 5;
S0 = gslew*100;
dt = dts*.5;

% slew-rate limited approximation

Ts = .666667 / nl*sqrt(((pi*N)^3)/(gamma*D*S0));
if (Ts > Tmax), display('slew limited readout too long'); return; end

a2 = N*pi/(nl*(Ts^(.666667)));
a1 = 1.5*S0/a2;
beta = S0*gamma*D/nl;
Gmax = a1*(Ts^.333333);
% gmax = 0;

printf('Ts=%g', Ts)

t = 0:dt:Ts;
x = t.^1.333333;
theta = (t.^2).*(.5*beta./(q + .5*beta./a2.*x));
y = q+.5.*beta./a2.*x;
dthdt = t.*(beta.*(q+.166667*beta./a2.*x)./(y.*y));
c = cos(theta);
s = sin(theta);
gx = (nl/(D*gamma)).*dthdt.*(c - theta.*s);
gy = (nl/(D*gamma)).*dthdt.*(s + theta.*c);
gabs = abs(gx+i.*gy);
% cut short if over peak
gmax = abs(gamp./(theta+eps) + i.*gamp);
l1 = length(t) - sum(gabs>gmax);
ts = t(l1);
thetas = theta(l1);


% gmax limited approximation

l3 = 0;
% T=ts;
if Gmax > gamp
  T=((pi*N/nl)*(pi*N/nl) - thetas*thetas)/(2*gamma*gamp*D/nl)+ts;
  if T > Tmax
      sprintf('gmax limited readout too long')
      return;
  end
  t = ts+dt:dt:T;
  theta = sqrt(thetas*thetas + (2*gamma*gamp*D).*(t-ts)./nl);
  c = cos(theta);
  s = sin(theta);
  ind2 = l1+(1:length(t));
  gx(ind2) = gamp.*(c./theta - s);
  gy(ind2) = gamp.*(s./theta + c);
  l3 = length(t);
end

l2 = l1 + l3;
Gx = gx(1:2:l2);
Gy = gy(1:2:l2);

len = length(Gx);

% Bring gradients back
areax = sum(Gx)*dts;
gxrew = trapwave(-areax, dts, gamp/sqrt(3), S0/sqrt(3),Gx(end));
Gx = [Gx gxrew];

areay = sum(Gy)*dts;
gyrew = trapwave(-areay, dts, gamp/sqrt(3), S0/sqrt(3),Gy(end));
Gy = [Gy gyrew];

lx = numel(Gx);
ly = numel(Gy);
if lx > ly
	Gy = [Gy zeros(1,lx-ly)];
else
	Gx = [Gx zeros(1,ly-lx)];
end

Gx = [Gx zeros(1,mod(length(Gx),2))];
Gy = [Gy zeros(1,mod(length(Gy),2))];

g = Gx + i.*Gy;				% grad vector
s = diff(g)./(gts*1000);	% slew rate vector
Kx = cumsum([0 Gx])*gts*opfov*gambar;
Ky = cumsum([0 Gy])*gts*opfov*gambar;
k = Kx + 1i*Ky;				% kspace vector

kx = real(k);
ky = imag(k);
sx = real(s);
sy = imag(s);

idx = [true(1,len+1) false(1,length(Gx)-len)];

end

function make_spiral_test
fov = 24;
N = 64;
dt = 4e-6;

[ks, om, gr, sl, ts, idx] = make_spiral(fov,N,dt,'type','std','shot',2,'rot',0);
[ks2, om2, gr2, sl2, ts2, idx2] = make_spiral(fov,N,dt,'type','std','shot',2,'rot',2);

% tst = ts(1:end-1);
% tst2 = ts2(1:end-1);

ks = [ks;ks2];
om = [om;om2];
gr = [gr;gr2];
sl = [sl;sl2];
ts = [ts;ts2];
% tst = [tst;tst2];
idx = [idx;idx2+length(ks2)];

ts = 1000*ts;


% figure, plot(ts, gr(:,1), '.b', ts, gr(:,2), '.r');
% figure, plot(ts, sl(:,1), 'b', ts, sl(:,2), 'r');
figure, plot(ks(:,1), ks(:,2), '.-'), axis square;
figure, plot(om(:,1), om(:,2), '.-'), axis square;

if ~isempty(idx)
	ks = ks(idx,:);
	gr = gr(idx,:);
	ts = ts(idx);
	figure, plot(ks(:,1), ks(:,2), '.-'), axis square;
% 	figure, plot(ts, gr(:,1), '.b', ts, gr(:,2), '.r');
end

end

