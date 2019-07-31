 function fit = se_ftab_fit(sl, fm, varargin)
%function fit = se_ftab_fit(sl, fm, [options])
%
% Fit to SE table fm(), suitable for subsequent interpolation / extrapolation.
% We can get expressions for derivatives / curvatures
% Uses an experimental exponential model -log(sum_k p_k exp(-m_k . s))
% eq below eq(18)
%
% in
%	sl	{L}		sample locations for each of L materials (cell)
%	fm	[s1,...,sL,M]	DE tables f_m(s_1,...,s_L)
% option
%	'type'	'exp'	type of fit (default: 'exp')
%	'wt'	{M}		fit weighting for each of M energies.
%	'show'	1|0		plot?
% options for exp (to get mac)
%	'kev'			kev (default: [10:5:200])          
%	'mtype'			mtype (required!) eq below eq(18), only alpha is estimated
% out
%	fit	strum		strum object for fitted f_m(s_1,...,s_L)
%				fit.coef is [nbasis, M]
%	methods:
%	fit.fmfun(sll)		fm function evaluation, for stacked array sll
%	fit.fgrad(sll)		fm gradient evaluation, for stacked array sll 
%	fit.show_sp(en, sp)	plot true spectrum vs fitted spectrum
%	fit.show_fm(sl, fm)	mesh plot of fm and its fit
%	fit.show_err(sl, fm)	mesh plot of fit error
%	fit.mac_eff		effective mass atten coef based on fit [MM, LL]
%				(valid for 'exp' only, at key energies)                       
%
% Copyright 2006-3-3, Jeff Fessler, The University of Michigan
% Modified by Joonki Noh 2007-11-3,

if nargin == 1 && streq(sl, 'test'), se_ftab_fit_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.show = false;
arg.type = '';
arg.kev = [10:5:200]';
arg.mtype = {};
arg.wt = {};
arg = vararg_pair(arg, varargin);
if isempty(arg.type), arg.type = 'exp'; end

if ndims(fm) == length(sl)
	MM = 1; %se-ct
elseif ndims(fm) == length(sl) + 1
	MM = size(fm,ndims(fm)); % last dimension is number of spectra
else
	error 'invalid fm size'
end %valid for any Kvp ct with LL = MM

for ll=1:length(sl) % make sure arguments match table dimensions
	if length(sl{ll}) ~= size(fm, ll)
		error('dim mismatch %d', ll)
	end
end

%if length(sl) ~= 2, warning 'todo: only L=2 materials tested', end

%%%%%%%%%%%%%%%%%%%%
% DE Table Fitting %
%%%%%%%%%%%%%%%%%%%%
switch arg.type
case 'exp'
	[fit fmfun fgrad] = se_ftab_fit_exp(sl, fm, MM, arg.wt, arg.mtype, arg.kev);
	fit.mac_eff = zeros(MM, size(fit.mac{1}, 2)); % [M,L] only for 'exp' approx 
	for mm=1:MM
		fit.mac_eff(mm,:) = fit.mac{mm}' * fit.coef{mm}; 
        %meaning of mac_eff : eq(16) on p6 [MM,LL] / indep of ith ray
	end
otherwise
	error('unknown fit type %s', arg.type)
end

fit.MM = MM;
fit.type = arg.type;

meth = {'fmfun', fmfun, '(sll)'; ... %cell 2nd arg is function handle
        'fgrad', fgrad, '(sll)'; ...
        'show_err', @se_ftab_fit_show_err, '(sl, fm)'; ...
        'show_fm', @se_ftab_fit_show_fm, '(sl, fm)'; ...
        'show_sp', @se_ftab_fit_show_sp, '(en, sp)'; ...
       }; %format for strum
fit = strum(fit, meth);  %strum object

if arg.show
	fit.show_fm(sl, fm);
end

end % se_ftab_fit()

%
% se_ftab_fit_exp()
%
function [fit, ffun, fgrad] = se_ftab_fit_exp(sl, fm, MM, wt, mtype, kev)

%sll = ndgrid_jf('mat', sl{:}); %[(Nd) LL] 3D
sll = sl';

if isempty(mtype), error 'mtype required', end
if isempty(wt), wt = {1,1}; end

LL = 1;Ab = 1;
mac = xray_read_atten(mtype, kev); % [ne, L]
for ll=1:LL
	sl = col(stackpick(sll, ll)); %[*Nd, 1]
	Ab = Ab .* exp(-sl * mac(:,ll)'); % [*Nd ne]
end

fit.kev = cell(1,MM);
fit.mac = cell(1,MM);
fit.coef = cell(1,MM);
for mm=1:MM
	dat = fm(:,mm); %for a given i-th index
	y = exp(-dat); %[(Nd)] 2D
	Wh = spdiag(sqrt(wt{mm}(:)), 'nowarn');
	x = wls_simplex(Ab, y(:), Wh); % coefficients x : [ne, 1]
    % min_x ||Wh * (Abx - y(:))||
    % subject to 0 <= x <= 1 and sum(x) = 1

	ie = x > 1e-6; % find key energies
	fit.kev{mm} = kev(ie);
	fit.mac{mm} = mac(ie,:); % [nke,L], nke = # of key energies
    %fit.mac{mm} is sampled from mac but depends on mth ray
    %If it combined with fit.coef, then the approximation of integral is fine
	A = 1;
	for ll=1:LL
		sl = col(stackpick(sll,ll));
		A = A .* exp(-sl * fit.mac{mm}(:,ll)'); % [*Nd, nke]
	end
	fit.coef{mm} = wls_simplex(A, y(:), Wh); %[nke,1] 
    % coeffs depend on mth incident spectrum
end
ffun = @se_ftab_fit_exp_eval;
fgrad = @se_ftab_fit_exp_grad;
end % se_ftab_fit_exp()


%
% se_ftab_fit_exp_eval()
% evaluate 
% in
%	sll	[(Nd),L]	stackup of s1,s2,...,s_L
% out
%	f	[(Nd),M]	stackup of f1,f2,...,f_M
%
function f = se_ftab_fit_exp_eval(fit, sll)
Nd = size(sll); LL = Nd(end); Nd = Nd(1:end-1);
sll = reshape(sll, [], LL); % [*Nd,L]
MM = fit.MM;
f = zeros(prod(Nd),MM);
for mm=1:MM
	A = 1;
	mac = fit.mac{mm};
	for ll=1:LL
		sl = sll(:,ll);
		A = A .* exp(-sl * mac(:,ll)'); % [*Nd,nke] only for key energies
	end
	tmp = -log(A * fit.coef{mm}); % [*Nd,1], summed w.r.t key energy
	f(:,mm) = tmp;
end
f = reshape(f, [Nd MM]);
end % se_ftab_fit_exp_eval()


%
% se_ftab_fit_exp_grad()
% evaluate gradient of f for each of the given s vectors.
% in
%	sll	[(Nd),L]	stackup of s1,s2,...,s_L
% out
%	g	[(Nd),L,M]	stackup of gradients of f(s)
%
function g = se_ftab_fit_exp_grad(fit, sll)
Nd = size(sll); LL = Nd(end); Nd = Nd(1:end-1);
sll = reshape(sll, [], LL); % [*Nd,L]
MM = fit.MM;
g = zeros(prod(Nd), LL, MM);
for mm=1:MM
	A = 1;
	mac = fit.mac{mm}; % [kne,L]
	alf = fit.coef{mm}; % [kne,1]
	for ll=1:LL
		sl = sll(:,ll);
		A = A .* exp(-sl * mac(:,ll)'); % [*Nd,kne]
	end
	vm = A * alf; % [*Nd,1] for den of grad
	tmp = A * (mac .* repmat(alf, [1 LL])); %[*Nd,kne]*[kne,L] = [*Nd,L]
	g(:,:,mm) = tmp ./ repmat(vm, [1 LL]);
end
g = reshape(g, [Nd LL MM]);
end % se_ftab_fit_exp_grad()


%
% se_ftab_fit_show_sp()
% compare true spectra to fitted spectra(fitted coeffs)
% fit.coeff is a rough approx to normalized spectra
function se_ftab_fit_show_sp(fit, en, sp)
if nargin < 3, error 'se_ftab_fit_show_sp(fit, en, sp)', end

if ~streq(fit.type, 'exp'), printm 'show_sp only done for exp', return, end
if im
	clf, pl = (fit.MM+1)*100 + 10 + 1;
	subplot(pl)
	plot(en, sp * diag(1 ./ max(sp)))
	for mm=1:fit.MM
		subplot(pl+mm)
		bar(fit.kev{mm}, fit.coef{mm}) %just bar graph plot
		axis tight, axisx(minmax(en))
	end
end
end % se_ftab_fit_show_sp()


%
% se_ftab_fit_show_fm()
% compare fit to sampled fm
%
function se_ftab_fit_show_fm(fit, sl, fm)
if nargin < 3, error 'se_ftab_fit_show_fm(fit, sl, fm)', end

sll = ndgrid_jf('mat', sl{:});

if fit.MM ~= 2, error 'only M=2 done', end

fh = fit.fmfun(sll); %approxed fm

s1 = sl{1};
s2 = sl{2};
smax(1) = max(s1);
smax(2) = max(s2);
fmax = max(fm(:)); %exact fm

ax = [0 smax(1) 0 smax(2) 0 fmax];
cax = [0 fmax];

im clf
im pl 2 2

show(1, s1, s2, fm(:,:,1), ax, cax, 'f_1(s)')
% text(-60, 10, '[cm^2/g]')
show(2, s1, s2, fm(:,:,2), ax, cax, 'f_2(s)')

show(3, s1, s2, fh(:,:,1), ax, cax, 'f_1 approx')
show(4, s1, s2, fh(:,:,2), ax, cax, 'f_2 approx')
end % se_ftab_fit_show_fm()


%
% se_ftab_fit_show_err()
% show fit errors, and relate to HU
% f = mac s, mac(70kev,H2O) = 0.2 cm^2 / g = 1000 HU
% mh = mac * 1000 HU / (0.2 g/cm^2)
% so Df/Dmh = Df/Dmac * Dmac/Dmh = (50 g/cm^2) (0.2 cm^2/g) / 1000 HU = 1/100 HU
% so Dmh = 100 HU * Df for 50 cm of H2O.
%
function se_ftab_fit_show_err(fit, sl, fm)
if nargin < 3, error 'se_ftab_fit_show_err(fit, sl, fm)', end

sll = ndgrid_jf('mat', sl{:});

fh = fit.fmfun(sll); %approxed fm
err = fh - fm; %[(Nd) MM]
printm('worst model error = %g of %g', max(abs(err(:))), max(fm(:)))
printm('worst model error = %g HU over 50cm H2O', 100*max(abs(err(:))))
for mm=1:fit.MM
	ee = stackpick(err,mm);
	printm('worst error (mm=%d): %g', mm, max(abs(ee(:))))
end

if max(abs(err(:))) == 0, return, end
if fit.MM > 2, error 'only M=2 done', end
if 0
	e1 = err(:,:,1);
	e2 = err(:,:,2);
	disp([minmax(e1); minmax(e2)]')
	printm('worst error1 %g', max(col(abs(err(:,1,:)))))
	printm('worst error2 %g', max(col(abs(err(:,2,:)))))
end
%err = abs(err);

s1 = sl{1};
s2 = sl{2};
smax(1) = max(s1);
smax(2) = max(s2);

elim = minmax(err(:))';
%elim = [-1 1] * 0.01; % +/- 1 HU
ax = [0 smax(1) 0 smax(2) elim];
im plc 1 2
for mm=1:fit.MM
	show(mm, s1, s2, err(:,:,mm), ax, elim, sprintf('f_%d error', mm))
end
end % se_ftab_fit_show_err()


%
% show() : Fessler's version mesh
%
function show(pl, x, y, f, ax, cax, ti)
if ~im, return, end
im('subplot', pl)
if 1
	mesh(x,y,f')
	colormap hsv, caxis(cax), cbar
	axis(ax)
	xtick, ytick%, ztick, zwhite
else
	im(x,y,f), cbar
	xtick, ytick
end
xlabel 's_1', ylabel 's_2', title(ti)
end % show()


%
% se_ftab_fit_test()
%
function se_ftab_fit_test
sl{1} = linspace(0, 50, 26); %For a given lth material,
sl{2} = linspace(0, 30, 31); %makes sinogram of (Nd) size
%stype = 'mono,70';
stype = 'ps1';
xray = xray_read_spectra(stype);
mtype = {'water', 'bone'};
mac = xray_read_atten(mtype, xray.en);
if im
	clf, semilogy(xray.en, mac), legend(mtype{:})
end
if length(sl) == 1
    sll = sl; %cell
elseif length(sl) ~= 1
    sll = ndgrid_jf('mat', sl{:}); %mat
end
fm = de_ftab_fm(sll, mac(:,(1:length(sl))), xray.Ide(:,(1:length(sl))));
fit = se_ftab_fit(sl, fm, 'type', 'exp', 'mtype', mtype); %sl should be cell
g = fit.fgrad(sll);
%fit = se_ftab_fit(sl, fm, 'type', 'poly')

if 1
	fit.show_sp(xray.en, xray.sp)
	prompt
	fit.show_fm(sl, fm)
	prompt
	fit.show_err(sl, fm)
end

end % se_ftab_fit_test()
