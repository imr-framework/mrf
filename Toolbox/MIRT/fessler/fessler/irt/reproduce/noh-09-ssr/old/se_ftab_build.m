 function ftab = se_ftab_build(s_arg, varargin)
%function ftab = se_ftab_build(s_arg, [options])
%
% For multiple-kVp X-ray imaging, we must evaluate functions f_m(s_1, ..., s_L)
% for m=1,..,M, where M is the number of kVp settings,
% and L is the number of material components.
% f_m(s1, s2) = -log( \int exp(- (m1(E) s1 + m2(E) s2)) dI_m(E) / I_m(E) )
% (For "dual-energy" imaging, M=L=2.)
%
% This routine builds tables/models of f_m, its inverse, and its derivatives.
% These are needed for (dual-energy) X-ray CT reconstruction.
%
% in:	(these all have sensible defaults, so try calling with no args)
%	s_arg	cell	arguments for xray_read_spectra ex) 'ps1'
% option
%	'mtype	cell	material types, e.g., {'soft', 'bone'}
%	'ftype	char	fit type for de_ftab_fit()
%	'sl	cell{LL} sample thicknesses for each of LL materials	
%	's_n	[L,1]	number of samples of l'th material integrals (default exist)
%	's_max	[L,1]	maximum material density line integrals, units: g/cm^2
%	'wt_fm	{M}	weighting for fitting fm
%
% out:
%	ftab	struct
%	methods:
%		ftab.fit.fmfun(s1, s2, ...)
%
% Copyright 2001-04-27, Jeff Fessler, The University of Michigan
if nargin < 1, help(mfilename), error(mfilename), end

ftab.show = false;
ftab.sl = {};
ftab.s_n = [];
ftab.s_max = [];
ftab.wt_fm = {};
ftab.stop_after_fit = false;
ftab.mtype = {'soft', 'bone'};	% default materials
ftab.ftype = ''; % defer to de_ftab_fit() default ('exp')
ftab = vararg_pair(ftab, varargin);

if isempty(s_arg)
%	s_arg = {'mono,60,100'}; % monoenergetic for testing
	s_arg = {'ps1'};
% fdir = sprintf('../fig,%s/', s_arg{1});
end

% material integral sampling:
% s_l is samples of the integral of the density of the lth material type.
[ftab.sl ftab.s_n ftab.s_max] = de_ftab_build_sl(ftab.sl, ftab.s_n, ftab.s_max); 
LL = length(ftab.mtype); % # of material components

% Load X-ray spectra, compute effective energy etc.
ftab.xray = xray_read_spectra(s_arg{:}, 'show', ftab.show);
MM = size(ftab.xray.sp, 2);	% # of spectra
if ftab.show, prompt, end

% mac: mass attenuation coefficients
ftab.mac = de_ftab_build_mac(ftab.mtype, ...
	ftab.xray.en, ftab.xray.Ide, ftab.xray.eff, ftab.show); 

% Build tables of f_m(s1, s2)
% f_m looks very linear but it is not quite linear!
ftab.fm_fun = @de_ftab_build_fm;
ftab.fm = ftab.fm_fun(ftab, ftab.sl{:});
if ftab.show % figure showing f1, f2
	de_ftab_build_show_fm(ftab.fm, ftab.sl, 4)
prompt
end

% parametric fit to each fm to form a continuous function!
ftab.fit = de_ftab_fit(ftab.sl, ftab.fm, 'type', ftab.ftype, ...
	'mtype', ftab.mtype{length(ftab.sl)}, 'wt', ftab.wt_fm);
%ftab.fit = de_ftab_fit(ftab.sl, ftab.fm, 'type', 'poly', 'order', 3, 'dc', 1);
if ftab.show % compare fit to sampled fm
	ftab.fit.show_fm(ftab.sl, ftab.fm)
	ftab.fit.show_err(ftab.sl, ftab.fm)
end

% for safefty, remove table, forcing user to use ftab.fit.fmfun!
% ftab = rmfield(ftab, 'fm');

if ftab.stop_after_fit, return, end

warning 'todo: not done after here'
return %Upto here till now

 end % de_ftab_build() - necessary??


%
% de_ftab_build_sl()
%
function [sl, s_n, s_max] = de_ftab_build_sl(sl, s_n, s_max)

% number of samples of the material integrals "s"
if isempty(s_n)
	if isempty(sl)
		s_n = [45 43]; %for spie02.m # of samples defined
	else
		for ll=1:length(sl)
			s_n(1,ll) = length(sl{ll});
		end
	end
end

% maximum material "integrals"
if isempty(s_max)
	if isempty(sl)
		% soft max: 50cm * 1g/cc
		% bone max: 15cm * 2g/cc (for now)
		s_max = [50 30]; %max of s_{il} defined
	else
		for ll=1:length(sl)
			s_max(1,ll) = max(sl{ll});
		end
	end
end

if isempty(sl)
	for ll=1:length(s_n)
		sl{ll} = linspace(0, s_max(ll), s_n(ll))';
	end
end

end % de_ftab_build_sl()


%
% de_ftab_build_mac()
%
function mac = de_ftab_build_mac(type, en, Ide, xeff, show)
mac.type = type;

% interpolate mass atten coef to source energy sampling
mac.mac = xray_read_atten(type, en); % [ne,L]
mac.bar = diag(1 ./ sum(Ide)) * (Ide' * mac.mac);

LL = length(type);
MM = size(Ide,2);

% examine condition number
if 1
	printm 'bmassml = '
	disp(mac.bar)
	% inv(mass.bar)
%	printm('condition = %g', cond(mass.bar))
	pr cond(mac.bar)
	disp( inv(mac.bar' * mac.bar) )
end

if show
	clf
	plot(	en, mac.mac(:,1), 'c:', ...
		en, mac.mac(:,2), 'y:', ...
		xeff(1), mac.bar(1,:), 'g+', ...
		xeff(2), mac.bar(2,:), 'r^')
	axisy(0.1, max(0.7, 1.1*max(mac.bar(:))))
	legend(mac.type{:})
prompt
end

end % de_ftab_build_mac()


%
% de_ftab_build_fm()
% build tables of f_m(s1, s2, ...)
%
function fm = de_ftab_build_fm(ftab, varargin)
if length(varargin) == 1 && iscell(varargin{1})
	sll = varargin{1};
else
	sll = ndgrid_jf('cell', varargin{:});
end
fm = de_ftab_fm(sll, ftab.mac.mac(:,(1:length(varargin))), ftab.xray.Ide(:,(1:length(varargin))));

end % de_ftab_build_fm()


%
% de_ftab_build_show_fm()
%
function de_ftab_build_show_fm(fm, sl, down)
fmax = 1.01 * max(fm(:));
fmax = ceil(fmax);
s_max(1) = max(sl{1}(:));
s_max(2) = max(sl{2}(:));

i1 = 1:down:length(sl{1});
i2 = 1:down:length(sl{2});
s1 = sl{1}(i1);
s2 = sl{2}(i2);
f1 = fm(i1,i2,1);
f2 = fm(i1,i2,2);

clf
subplot(221), mesh(s1, s2, f1')
colormap hsv, caxis([0 fmax]), cbar
axis([0 s_max(1) 0 s_max(2) 0 fmax])
xtick, ytick, ztick, zwhite, xlabel s_1, ylabel s_2, title 'f_1(s)'
text(-60, 10, '[cm^2/g]')

subplot(222), mesh(s1, s2, f2')
colormap hsv, caxis([0 fmax]), cbar
axis([0 s_max(1) 0 s_max(2) 0 fmax])
xtick, ytick, ztick, zwhite, xlabel s_1, ylabel s_2, title 'f_2(s)'

% ir_savefig('c', fdir, 'fig_f1_f2')
end % de_ftab_build_show_fm()
