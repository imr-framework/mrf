 function [sh cost_eval_ls cost_eval_penal NRMS] = de_ftab_s_pwls(fit,fh,yim,varargin)
%function [sh cost_eval NRMS] = de_ftab_s_pwls(fit,fh,yim,varargin)
% estimate s from fh by PWLS
% in
%	fit	from de_ftab_fit()
%	fh	[(Nd),M] estimates of f
%   yim [(Nd),M] noisy measurements
% option
%   'strue   [(Nd),L] true sinogram (only in simulation for DECT NRMS)
%	'niter	 # of iterations
%	'init	 [(Nd),L] initial estimates of s (default: linear inverse)
%   'beta    parameters for regulization (soft / bone) (default:1,1)
%   'regsel  selector for regularization
%            0 for conventional & 1 for quad modified (default: 0)
% out
%	sh                  [(Nd),L] estimates of s
%   cost_eval_ls        evaluation of cost function (LS part)
%   cost_eval_penal     evaluation of cost function (penalty part)
%   NRMS                component-wise nrms error for each iteration [# of iter,L]
% Copyright 2007-07-13, Jeff Fessler, The University of Michigan
% Modified by Joonki Noh 2007-11-04
%
if nargin == 1 && streq(fit, 'test'), de_ftab_s_pwls_test, return, end
if nargin < 2, error(mfilename), help(mfilename), end

arg.niter = 1; 
arg.strue = [];
arg.init = [];
arg.beta = [];  
arg.regsel = 0;
arg = vararg_pair(arg, varargin);

LL = size(fit.mac{1}, 2);Nd = size(fh);
if ndims(fh) == 3
    MM = Nd(end); Nd = Nd(1:end-1);
elseif ndims(fh) == 2
    MM = 1;
else
    error 'something wrong in matrix size'
end

fh = reshape(fh, [], MM); % [*Nd,M]
sm_ker = ones(3,3);sm_ker = sm_ker/sum(sum(sm_ker));
for rr=1:MM %simple smoothing
    ybi_hat(:,:,rr) = conv2(yim(:,:,rr),sm_ker,'same');
end
yim = reshape(yim,[],MM); % [*Nd,M]
ybi_hat = reshape(ybi_hat,[],MM); %[*Nd,M]
if MM ~= fit.MM, error 'M', end

if isempty(arg.init) %mac_eff is indep of ith ray
	sh = (pinv(fit.mac_eff) * fh')'; %inv[M,L]*[M,*Nd] = [L,*Nd]' 
    % linear inverse (perfect if mono) (Eq (17)) on p6
    sh = max(sh,0);
else    
	sh = max(arg.init, 0);
end

if arg.niter < 1, return, end 
if isempty(arg.beta),   arg.beta = ones(1,LL); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighted Least Square Part & Penalty Part %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LS (precomputed curvature)
[curv1 curv2] = de_ftab_curv(fit); %curv1 & curv2 : [M,1] 
fstep = @(fm, fh) (yim*curv1 + (max(fh-fm, 0).*yim)*curv2); %[*Nd,1]
fm = fit.fmfun(sh); %[*Nd,M]

%%%%%%%%%%%%%%%%%%
% Regularization % conventional / quad modified (based on Monochromatic spectrum)
%%%%%%%%%%%%%%%%%%
% Modified penalty ==> D1i, D12i, D2i: [*Nd] 
if MM ~= 1
    %mac_soft_vec = [0.2100 0.1896]; 
    %mac_bone_vec = [0.3393 0.2495];
    mac_soft_vec = fit.mac_eff(:,1)'; %From fitted table
    mac_bone_vec = fit.mac_eff(:,2)';
elseif MM == 1
    %mac_soft_vec = 0.1896;mac_bone_vec = 0.2495;
    mac_soft_vec = fit.mac_eff;mac_bone_vec = 0;
end

for jj=1:LL
    if arg.regsel == 0 %conventional regularizer
        R{jj} = Reg1(logical(ones(Nd)),'beta',arg.beta(jj),'offsets',[1],'type_denom','matlab');
    elseif arg.regsel == 1 %quad modified regularizer
        D1i = 0;D2i = 0;D12i = 0;
        for ii=1:MM %with smoothing of yim
             D1i = D1i + ybi_hat(:,ii)*mac_soft_vec(ii)^2; %soft tissue
             D2i = D2i + ybi_hat(:,ii)*mac_bone_vec(ii)^2; %bone
             D12i = D12i + ybi_hat(:,ii)*mac_soft_vec(ii)*mac_bone_vec(ii); %interaction
        end
        Di(:,:,1) = reshape(sqrt(D1i),Nd(1),Nd(2)); %kappa1
        Di(:,:,2) = reshape(sqrt(D2i),Nd(1),Nd(2)); %kappa2
        R{jj} = Reg1(Di(:,:,jj),'beta',arg.beta(jj),'offsets',[1],'type_denom','matlab');
    end
end %offsets,[1] option gives smoothing in the radial direction only

%%%%%%%%%%%%%%%%%%%%%%
% Main loop for PWLS % 
%%%%%%%%%%%%%%%%%%%%%%
ticker reset
for ii=1:arg.niter
    %%%%%%%%%%%%
    % Gradient %
    %%%%%%%%%%%%
    %WLS
	ticker(mfilename, ii, arg.niter)
	fgrad = fit.fgrad(sh); % [*Nd,L,M]
	tmp = repmat((fh - fm).*yim, [1 1 LL]); % [*Nd,M,L]
	tmp = permute(tmp, [1 3 2]); % [*Nd,L,M]
	fgrad = fgrad.*tmp;
    %Penalty 
    penal_cgrad = [];penal_denom = [];
    for jj=1:LL 
        penal_cgrad = double([penal_cgrad R{jj}.cgrad(R{jj},sh(:,jj))]); %[*Nd,L]
        penal_denom = double([penal_denom R{jj}.denom(R{jj},sh(:,jj))]); %[*Nd,L]
        %This works as "penal_cgrad = [R.soft.cgrad(R.soft,sh(:,1)) R.bone.cgrad(R.bone,sh(:,2))];"
    end
    
    %%%%%%%%%%%%%
    % Step Size %
    %%%%%%%%%%%%%
	step = fstep(fm, fh); %[*Nd,1]
  	step = repmat(step, [1 LL]); %[*Nd,L]
    step = 1./(step+penal_denom);

    %%%%%%%%%%
    % Update % Non-negativity constraint on sinogram !!
    %%%%%%%%%%
	sh = sh + step.*(sum(fgrad,3)-penal_cgrad); %[*Nd,L]
    sh = max(sh,0);
    fm = fit.fmfun(sh); % [*Nd,M]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NRMS error per iteration % from 1st iteration 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(arg.strue)
        sh_tmp = reshape(sh, [size(arg.strue,1) size(arg.strue,2) LL]);NRMS_tmp = [];
        for kk=1:LL
            NRMS_tmp = [NRMS_tmp nrms(sh_tmp(:,:,kk), arg.strue(:,:,kk))];
        end
        NRMS(ii,:) = NRMS_tmp;
        clear sh_tmp;
    else
        NRMS = zeros(arg.niter,LL);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Cost Function Eval % from 1st iteration
    %%%%%%%%%%%%%%%%%%%%%%
    Rtemp = 0;
    for jj=1:LL
        Rtemp = Rtemp + R{jj}.penal(R{jj},sh(:,jj));    
    end
    penal_eval = sum(Rtemp);  
    pwls_eval = sum(sum((yim.*(fh - fm).^2)/2));
    cost_eval_ls(ii) = pwls_eval;
    cost_eval_penal(ii) = penal_eval;
end

function [curv1, curv2] = de_ftab_curv(fit) %[M,1]
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
	curv1(mm) = norm(g0)^2; 
	curv2(mm) = norm(h0); 
    %curv1 and curv2 are independent of i (see eq(16))
end
