 function [sh cost_eval_li cost_eval_penal NRMS] = de_ftab_s_pl(fit,xray,mac,fh,yim,varargin)
%function sh = de_ftab_s_pl(fit, xray, fh, yim, varargin)
% estimate s from fh by PL based on possion model
% in
%	fit	from de_ftab_fit()
%   xray from xray_read_spectra
%   mac full mass atten coeff 
%	fh	[(Nd),M] estimates of f, necessary for initialization
%   yim [(Nd),M] noisy measurements
% option
%   'strue [(Nd),L] true sinogram (only in simulation for DECT NRMS)
%	'niter	   # of iterations
%	'init	   [(Nd),L] initial estimates of s (default: linear inverse)
%   'beta      parameters for regulization (soft / bone) (default : 1, 1)
%   'curvtype  'pc' or 'npc'
%   'regsel     selector for regularization
%               0 for conventional & 1 for quad modified (default: 0)
% out
%	sh              [(Nd),L] estimates of s
%   cost_eval_li    evaluation of cost function (likelihood part)
%   cost_eval_penal evaluation of cost function (penalty part)
%   NRMS            component-wise nrms error for each iteration [# of iter,L]
% Copyright 2007-07-18, Jeff Fessler, The University of Michigan
% Modified by Joonki Noh 2007-11-04
if nargin == 1 && streq(fit, 'test'), de_ftab_s_pwls_test, return, end
if nargin < 2, error(mfilename), help(mfilename), end

arg.niter = 1; 
arg.init = [];
arg.beta = [];
arg.strue = [];
arg.curvtype = '';
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
if isempty(arg.beta),   arg.beta = ones(1,LL);   end

%%%%%%%%%%%%%%
% Likelihood %
%%%%%%%%%%%%%%
% Precomputed curvature 
curvml = de_ftab_curv(fit,xray,mac); %[M,L] using bar{E}_m
%curvml = de_ftab_curv_jeff(fit) %[M,L] using effective mac
step_pc = yim*curvml; % [*Nd,L]
fm = fit.fmfun(sh); %[*Nd,M]

%%%%%%%%%%%%%%%%%%
% Regularization % conventional / quad modified (based on Monochromatic spectrum)
%%%%%%%%%%%%%%%%%% In real data, we may need to apply a spatial smoother ??
if MM ~= 1
    mac_soft_vec = fit.mac_eff(:,1)'; %From fitted table
    mac_bone_vec = fit.mac_eff(:,2)';
elseif MM == 1
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
        Di(:,:,1) = reshape(sqrt(D1i),Nd(1),Nd(2));
        Di(:,:,2) = reshape(sqrt(D2i),Nd(1),Nd(2));
        R{jj} = Reg1(Di(:,:,jj),'beta',arg.beta(jj),'offsets',[1],'type_denom','matlab');
    end
end %offsets,[1] option gives smoothing in the radial direction only
Im = xray.I; %[1,M] total intensity

%%%%%%%%%%%%%%%%%%%%
% Main loop for PL % 
%%%%%%%%%%%%%%%%%%%%
ticker reset
for ii=1:arg.niter
    %%%%%%%%%%%%
    % Gradient %
    %%%%%%%%%%%%
    ticker(mfilename, ii, arg.niter)
  	fgrad = fit.fgrad(sh); % [*Nd,L,M]
    %ybar & its gradient
    ybim = exp(-fm)*diag(Im); %[*Nd,M]*[M,M] = [*Nd,M] 
    ybgrad = repmat(ybim,[1 1 LL]); %[*Nd,M,L]
    ybgrad = permute(ybgrad,[1 3 2]); %[*Nd,L,M]
    ybgrad = ybgrad.*fgrad; %[*Nd,L,M]
    %Likelihood
    Lgrad = (1 - yim./ybim); %[*Nd,M]
    Lgrad = repmat(Lgrad,[1 1 LL]); %[*Nd,M,L]
    Lgrad = permute(Lgrad,[1 3 2]); %[*Nd,L,M]
    Lgrad = sum(Lgrad.*ybgrad,3); %[*Nd,L]
    %Penalty 
    penal_cgrad = [];penal_denom = [];
    for jj=1:LL
        penal_cgrad = double([penal_cgrad R{jj}.cgrad(R{jj},sh(:,jj))]); %[*Nd,L]
        penal_denom = double([penal_denom R{jj}.denom(R{jj},sh(:,jj))]); %[*Nd,L]
        %This works as "penal_cgrad = [R.soft.cgrad(R.soft,sh(:,1)) R.bone.cgrad(R.bone,sh(:,2))];"
    end
    
    %%%%%%%%%%%%%
    % Step size %
    %%%%%%%%%%%%%
    if streq(arg.curvtype,'npc')
        %non-precomputed step size
        curvilm = de_ftab_curv_npc(fit,xray,sh); %[*Nd,L,M]
        step_npc = yim./ybim; %[*Nd,M]
        step_npc = repmat(step_npc,[1 1 LL]); %[*Nd,M,L]
        step_npc = permute(step_npc,[1 3 2]); %[*Nd,L,M]
        step_npc = step_npc.*curvilm; %[*Nd,L,M]
        step = 1./(sum(step_npc,3) + penal_denom); %[*Nd,L]
    elseif streq(arg.curvtype,'pc')
        %precomputed step size
        step = 1./(step_pc + penal_denom); %[*Nd,L]
    end    
    
    %%%%%%%%%%
    % Update % Non-negativity constraint on sinogram !!
    %%%%%%%%%%
	sh = sh + step.*(Lgrad-penal_cgrad); %[*Nd,L]
    sh = max(sh,0);
    fm = fit.fmfun(sh); %[*Nd,M]
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NRMS error per iteration % from 1st iteration 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(arg.strue)
        sh_tmp = reshape(sh, [size(arg.strue,1) size(arg.strue,2) LL]);NRMS_tmp = [];
        for kk=1:LL
            NRMS_tmp = [NRMS_tmp Nrms(sh_tmp(:,:,kk), arg.strue(:,:,kk))];
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
    ybim = exp(-fm)*diag(Im); %[*Nd,M]*[M,M] = [*Nd,M] 
    L_eval = ybim - yim.*log(ybim); %[*Nd,M]
    cost_eval_li(ii) = sum(sum(L_eval));
    cost_eval_penal(ii) = penal_eval;
end

function curvml = de_ftab_curv(fit,xray,mac) %[M,L]
% mac [ne,L] for full energies, not key energies 
% based on p7 of my note before summing over m and yim multiplication
MM = fit.MM;LL = size(mac,2);curvml = zeros(MM,LL);
curv = sum(mac,2); %[ne,1]
curv = repmat(curv,[1,LL]); %[ne,L]
curv = mac.*curv; %[ne,L]
for mm=1:MM
    %For effective energy 
    eff_en = xray.eff(mm);
    tmp = abs(xray.en - eff_en);
    ind = find(min(tmp) == tmp); %index for effective energy
    if ~isscalar(ind),   error 'error to pick effective energy';    end
    curvml(mm,:) = curv(ind,:);
end

% returning integral w.r.t E
function curvilm = de_ftab_curv_npc(fit,xray,sh) %[*Nd,L,M]
MM = fit.MM;LL = size(fit.mac{1},2);
for mm=1:MM
    %pim(E) at key energies
    alf = fit.coef{mm}; %[ne,1] 
    alf = repmat(alf,[1 LL]); %[ne,L]
    %tilde mac
    mac = fit.mac{mm}; %[ne,L]
    tmac = sum(mac,2); %[ne,1]
    tmac = repmat(tmac,[1 LL]); %[ne,L]
    tmac = tmac.*mac; %[ne,L], tilde{\beta}_{l}(E)
    %integral wrt E
    tmp = tmac.*alf; %[ne,L]
    intE = exp(-sh*mac'); %[*Nd,L]*[ne,L]'=[*Nd,ne]
    intE = intE*tmp; %[*Nd,L]
    curvilm(:,:,mm) = intE*xray.I(mm); %[*Nd,L]
end

% de_ftab_curv_jeff() using effective mac instead of effec energy
% based on Jeff's de_pl_denom / no significant difference from de_ftab_curv
function curvml = de_ftab_curv_jeff(fit) %[M,L]
curv = zeros(size(fit.mac_eff));LL = size(fit.mac_eff,2);
curv = sum(fit.mac_eff,2); %[M,1]
curv = repmat(curv,[1 LL]); %[M,L]
curvml = curv.*fit.mac_eff;
