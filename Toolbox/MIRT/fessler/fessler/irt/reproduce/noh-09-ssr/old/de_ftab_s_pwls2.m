 function [sh cost_eval] = de_ftab_s_pwls2(fit, fh, yim, varargin)
%function sh = de_ftab_s_pwls2(fit, fh, yim, varargin)
% estimate s from fh by PWLS with a modified curvature
% in
%	fit	from de_ftab_fit()
%	fh	[(Nd),M] estimates of f
%   yim [(Nd),M] noisy measurements
% option
%	'niter	# of iterations
%	'init	[(Nd),L] initial estimates of s (default: linear inverse)
%   'beta   parameter for global regulization (default : 1)
% out
%	sh          [(Nd),L] estimates of s
%   cost_eval   evaluation of cost function
% Copyright 2007-07-13, Jeff Fessler, The University of Michigan
if nargin == 1 && streq(fit, 'test'), de_ftab_s_pwls_test, return, end
if nargin < 2, error(mfilename), help(mfilename), end

arg.niter = 1; 
arg.init = [];
arg.beta = [];
arg = vararg_pair(arg, varargin);

LL = size(fit.mac{1}, 2);
Nd = size(fh); MM = Nd(end); Nd = Nd(1:end-1); 
fh = reshape(fh, [], MM); % [*Nd,M]
yim = reshape(yim,[],MM); % [*Nd,M]
if MM ~= fit.MM, error 'M', end

if isempty(arg.init) %mac_eff is indep of ith ray
	sh = (pinv(fit.mac_eff) * fh')'; %inv[M,L]*[M,*Nd] = [L,*Nd]' 
    %sh = max(sh,0); %by me
    % linear inverse (perfect if mono) (Eq (17)) on p6
else    
	sh = max(arg.init, 0);
end

if arg.niter < 1, return, end 
if isempty(arg.beta),   arg.beta = 1;   end

fm = fit.fmfun(sh); % fm : [*Nd,M], sh: [*Nd,L]
% Penalty object
R = Reg1(logical(ones(Nd)),'beta',arg.beta,'offsets',[1],'type_denom','matlab');    

%%%%%%%%%%%%%
% Iteration %
%%%%%%%%%%%%%
ticker reset
for ii=1:arg.niter
    %%%%%%%%%%%%%
    % Curvature %
    %%%%%%%%%%%%%
    [curv1 curv2] = de_ftab_curv2(fit,sh); %[M,1]
    fstep = @(fm, fh) (yim*curv1 + (max(fh - fm, 0).*yim)*curv2); % [*Nd,1]
    %fstep = @(fm, fh) (yim*curv1 + (max(fh, 0).*yim)*curv2); % [*Nd,1]
    
    %%%%%%%%%%%%
    % Gradient %
    %%%%%%%%%%%%
    %LS
	ticker(mfilename, ii, arg.niter)
	fgrad = fit.fgrad(sh); % [*Nd,L,M]
	tmp = repmat((fh - fm).*yim, [1 1 LL]); % [*Nd,M,L]
	tmp = permute(tmp, [1 3 2]); % [*Nd,L,M]
	fgrad = fgrad .* tmp;
    %Penalty 
    penal_cgrad = [R.cgrad(R,sh(:,1)) R.cgrad(R,sh(:,2))]; %[*Nd,L]
    penal_denom = [R.denom(R,sh(:,1)) R.denom(R,sh(:,2))]; %[*Nd,L]  
    
    %%%%%%%%%%%%%
    % Step Size %
    %%%%%%%%%%%%%
	step = fstep(fm, fh); %[*Nd,1]
  	step = repmat(step, [1 LL]); %[*Nd,L]
    step = 1./(step+penal_denom);
%	minmax(step)

    %%%%%%%%%%%%%%%%%%%%%%
    % Cost Function Eval %
    %%%%%%%%%%%%%%%%%%%%%%
    penal_eval = sum(R.penal(R,sh)); %?? material by material 
    ls_eval = sum(sum((yim.*(fh - fm).^2)/2));
    cost_eval(ii) = ls_eval+penal_eval;

    %%%%%%%%%%
    % Update %
    %%%%%%%%%%
	sh = sh + step .* (sum(fgrad,3)-penal_cgrad); 
    %sh = max(sh,0); %non-negativity constraint by me
	fm = fit.fmfun(sh); % [*Nd,M]
%	cost = mean(col(fh - fm).^2)
end

% de_ftab_curv2()
% a modified curvature computation for pwls sinogram restoration
%
function [curv1, curv2] = de_ftab_curv2(fit,sh)
MM = fit.MM;[Nd LL] = size(sh); %[*Nd,LL]
curv1 = zeros(MM,1);curv2 = zeros(MM,1);
for mm=1:MM
    %%%%%%%%%%%%
    % 1st part % ||b_{im}||^2
    %%%%%%%%%%%%
	alf = fit.coef{mm}; % [ne,1] / rough approximation of p_m at key energies
	mac = fit.mac{mm}; % [ne,LL]
	g0 = mac' * alf; % [LL,1] / gradient of fit at s=0 (largest point)
    %eq(16) on p6 w/o the ith index
    curv1(mm) = norm(g0)^2; %vector norm

    %%%%%%%%%%%%
    % 2nd part % max_{s_i} ||-hessian{f_{im}}||
    %%%%%%%%%%%%
    tmp = exp(-mac*sh'); %[ne,LL]*[*Nd,LL]' = [ne,*Nd] 
    vim = tmp'*alf; %[*Nd,1] / def of vim on eq (14) 
    qim = repmat(alf,[1,Nd]).*tmp; %[ne,*Nd]
    qim = qim./ repmat(vim',[size(alf,1),1]); %[ne,*Nd]
    
    %gradient
    g = mac'*qim; %[LL,*Nd] grad of fim on eq(14)
    for ii=1:Nd %hessian
        h = mac'*diag(qim(:,ii))*mac - g(:,ii)*g(:,ii)'; %[LL,LL] / 2nd eq on p7
        curv2_tmp(ii) = norm(h); %matrix norm
    end
    curv2(mm) = max(curv2_tmp); %??
end

% de_ftab_curv()
% curv* are [M,1] upper bounds : Todo : curvature modification
%
function [curv1, curv2] = de_ftab_curv(fit)
MM = fit.MM;
curv1 = zeros(MM,1);
curv2 = zeros(MM,1);
for mm=1:MM
	alf = fit.coef{mm}; % [ne,1] / if 'exp' fit, [nke,1]
	mac = fit.mac{mm}; % [ne,L]
	g0 = mac' * alf; % [L,1] / gradient of fit at s=0 (largest point)
    %eq(16) on p6 w/o the ith index
	h0 = mac' * diag(alf) * mac - g0 * g0'; %[L,L] / hessian at s=0 (largest??)
    %2nd eq on p7
	curv1(mm) = norm(g0)^2; %vector norm
	curv2(mm) = norm(h0); %matrix norm 
    %curv1 and curv2 are independent of i due to eq(16)
end

%
% de_ftab_s_pwls_test
%
function de_ftab_s_pwls_test
sl{1} = linspace(0, 50, 26);
sl{2} = linspace(0, 30, 31);
stype = 'ps1';
xray = xray_read_spectra(stype);
mtype = {'boron', 'iron'};
mac = xray_read_atten(mtype, xray.en);
sll = ndgrid_jf('mat', sl{:});
fm = de_ftab_fm(sll, mac, xray.Ide);
fit = de_ftab_fit(sl, fm, 'type', 'exp', 'mtype', mtype);

sh = de_ftab_s_iter(fit, fm, 'niter', 0); % initialize with linear inv.
%In fact, fm is the place of fhat
st = reshape(sll, [], 2);
pr rms(col(sh - st))
if 1 % picture of error of linear inv. approach vs (s_1,s_2)
	fit.mac_eff
	pr cond(fit.mac_eff)
	tmp = reshape(sqrt(mean((sh - st).^2, 2)), size(sll(:,:,1)));
	if im
		im(sl{:}, tmp), cbar % error is smallest at (0,0)
		xlabel(mtype{1}), ylabel(mtype{2})
	end
prompt
end

% todo: run profiler
sh = de_ftab_s_iter(fit, fm, 'init', sh, 'niter', 3000);
pr rms(col(sh - st))

shp = reshape(sh, size(sll));
if im
	plot(sl{1}, shp(:,:,1), 'c', sl{2}, shp(:,:,2), 'y')
	grid, axis equal, axis square
	xlabel 'true', ylabel 'noiseless estimate'
end
