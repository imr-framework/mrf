%% Initialization parameters

% Values for parameter mu
f.mu = 2.^(-2:-1:-8);
% f.mu = 2^-6;
% f.mu = 2.^(-6:-1:-9);
f.nmu = length(f.mu);

% Values for parameter nu
f.nu = 2.^(-1:-1:-9);
% f.nu = 2.^(-7:-4);
f.nnu = length(f.nu);

f.lam = 2^-15; % regularization parameter
f.reg_type = 'tv,aniso';
f.niter = 4000; % # of iterations
f.plot = 0; % plot results
f.chat = 1; % show printouts


%% AL-P2 params
%
% TV reg BSNR 30, lam=2^-12
% mu = 2^-5, nu = 2^-3 (mu = 2^-4, nu = 2^-4 or mu = 2^-3, nu = 2^-5)
% In all cases mu*nu = 2^-8 and lam/(mu*nu) = 2^-4
% 
% TV reg BSNR 40, lam=2^-15
% mu = 2^-4, nu = 2^-6 (mu = 2^-5, nu = 2^-5 or mu = 2^-3, nu = 2^-7)
% In all cases mu*nu = 2^-10 and lam/(mu*nu) = 2^-5
% 
% TV reg BSNR 50, lam=2^-17
% mu = 2^-5, nu = 2^-6 (mu = 2^-4, nu = 2^-7)
% In all cases mu*nu = 2^-11 and lam/(mu*nu) = 2^-6
% 
% 
% Wavelet reg BSNR 30, lam=2^-12 (retune) (not known)
% mu = 2^-4, nu = 2^0 (mu = 2^-5, nu = 2^1 or mu = 2^-3, nu = 2^-1 or mu = 2^-2, nu = 2^-2)
% In all cases mu*nu = 2^-4 and lam/(mu*nu) = 2^-9
% 
% Wavelet reg BSNR 40, lam=2^-15
% mu = 2^-4, nu = 2^-3 (mu = 2^-5, nu = 2^-2 or mu = 2^-3, nu = 2^-4)
% In all cases mu*nu = 2^-7 and lam/(mu*nu) = 2^-8
% 
% Wavelet reg BSNR 50, lam=2^-17 (retune) (not known)
% mu = 2^-5, nu = 2^-2 (mu = 2^-4, nu = 2^-3 or mu = 2^-3, nu = 2^-4)
% In all cases mu*nu = 2^-7 and lam/(mu*nu) = 2^-8

%% AL-P1 params
%
% For all setting use circulant preconditioner for CG
% 2 CG iterations are enough (not much gain from more iter) bc of good precon
% Even only 1 CG iter could work
%
% TV reg BSNR 30, lam = 2^-12
% mu = 2^-8 (lam/mu = 2^-4)
%
% TV reg BSNR 40, lam = 2^-15
% mu = 2^-10 (lam/mu = 2^-5)
%
% TV reg BSNR 50, lam = 2^-17
% mu = 2^-11 (lam/mu = 2^-6)
%
%
% Wavelet reg BSNR 30, lam = 2^-12 (not known)
% mu = 2^-4 (lam/mu = 2^-9)
%
% Wavelet reg BSNR 40, lam = 2^-15
% mu = 2^-7 (lam/mu = 2^-8)
%
% Wavelet reg BSNR 50, lam = 2^-17 (not known)
% mu = 2^??
%

%% SALSA params
%
% For all setting use circulant preconditioner for CG
% 2 CG iterations are enough (not much gain from more iter) bc of good precon
% Could use 4 CG iterations for lam=2^-15
%
% TV reg BSNR 30, lam = 2^-12
% mu = 2^-7 (lam/mu = 2^-5)
% Needs 8-10 Chambolle iterations (more is not beneficial)
%
% TV reg BSNR 40, lam = 2^-15
% mu = 2^-8 (lam/mu = 2^-7)
% Needs 8-10 Chambolle iterations (more is not beneficial)
%
% TV reg BSNR 50, lam = 2^-17
% mu = 2^-10 (lam/mu = 2^-7)
% Needs 8-10 Chambolle iterations (more is not beneficial)
%
%
% Wavelet reg BSNR 30, lam = 2^-12 (not known)
% mu = 2^-4 (lam/mu = 2^-9)
% Needs 8-10 Chambolle iterations (more is not beneficial)
%
% Wavelet reg BSNR 40, lam = 2^-15
% mu = 2^-8 (lam/mu = 2^-7)
% Needs 8-10 Chambolle iterations (more is not beneficial)
%
% Wavelet reg BSNR 50, lam = 2^-17 (not known)
% mu = 2^??
% Needs 8-10 Chambolle iterations (more is not beneficial)
%

%% (M)FISTA params
%
% (M)FISTA needs a lot of chambolle iterations to converge to actual minimizer
% Using 1-4 inner iterations gets stuck to wrong minimizer
%
% For 10 and more inner iter becomes too costly
% 

%% NCG params
% 
% Use 4 line search iterations (not much benefit from using more)
% Use circulant preconditioner. Convergence is slow anyway and there's no need
% for better diag-circ preconditioner
% Use epsilon = 10^-8 for the rounding parameter
% 

%% Create blur and masks, and load true image and data
load ../data/mat/xinf_anal_exact_tv_bsnr40_mot9-30.mat

% xtrue = xres(:,:,1);
xtrue = mean(xres,3);
f.climx = [0 max(xtrue(:))+0.1];

f.nb = 9;
% f.psf = ones(f.nb)/f.nb^2;
% f.psf = fspecial('gaussian', f.nb, 1.5); % Gaussian blur with sigma=2
f.psf = padarray(fspecial('motion', f.nb, 30), [2 0]); % Gaussian blur with sigma=1.5

[nxr, nyr] = size(xtrue);
f.mask_res = true(nxr,nyr);

nx = nxr + f.nb - 1;
ny = nyr + f.nb - 1;
f.mask_obj = true(nx,ny);
f.mask = false(nx,ny);
f.mask((f.nb+1)/2:end-(f.nb-1)/2,(f.nb+1)/2:end-(f.nb-1)/2) = true;


%% Setup and run experiments

nrms_lim = nrms(xres(:,:,1),xtrue);
nrms_dblim = 20*log10(nrms_lim);

printm('NRMSE bound =%6.3fx10^-4, in db =%6.1f\n', nrms_lim*1e4, nrms_dblim);

xres = zeros(nxr,nyr,f.nmu,f.nnu);
err = zeros(f.niter+1,f.nmu,f.nnu);
time = zeros(f.niter+1,2,f.nmu,f.nnu);

for ii=1:f.nmu
	for jj=1:f.nnu
% 		jj=ii;
		[xr out] = image_restore(y, f.psf, 'est_type', 'alp2', ...
			'model_type', 'non-circ', 'lambda', f.lam, 'reg_type', f.reg_type, ...
			'mu', f.mu(ii), 'nu', f.nu(jj), 'xtrue', xtrue, 'niter', f.niter, ...
			'plot', f.plot, 'chat', f.chat);
		
		xr = embed(xr(f.mask),f.mask_res);
		xres(:,:,ii,jj) = xr;
		err(:,ii,jj) = out.err;
		time(:,:,ii,jj) = out.time;
		
		figure(1), im(xtrue,'Original image','cbar',f.climx), drawnow;
		figure(2), im(y,'Blured noisy image','cbar'), drawnow;
		figure(3), im(xr,'Restored image','cbar',f.climx), drawnow;
		figure(4), im(xr-xtrue,'Difference image','cbar'), drawnow;
		
		printm('Image NRMS in db: %6.4f\n',20*log10(nrms(xr(:), xtrue(:))));
	end
end

% Find parameters mu, nu that give the best possible outcome
errm = squeeze(err(end,:,:));
[errmin ind] = min(errm(:));

[mu,nu] = ndgrid(f.mu,f.nu);
mu_min = mu(ind);
nu_min = nu(ind);

printm('Parameters for minimum NRMSE: mu=2^%d, nu=2^%d\n', ...
	log2(mu_min), log2(nu_min));

% plot for different mu values
for ii=1:f.nmu
	figure(ii), plot(0:f.niter,squeeze(err(:,ii,:)),...
		0:f.niter,nrms_dblim*ones (1,f.niter+1),'b--')
	title(sprintf('Error per iteration for mu=2^%d',log2(f.mu(ii))));
end

% plot for different nu values
% for jj=1:f.nnu
% 	figure(f.nmu+jj), plot(0:f.niter,squeeze(err(:,:,jj)),...
% 		0:f.niter,nrms_dblim*ones(1,f.niter+1),'b--')
% 	title(sprintf('Error per iteration for nu=2^%d',log2(f.nu(jj))));
% end


printm('Actual NRMSE for all reconstructions\n');
for ii=1:f.nmu*f.nnu
	printm('mu=2^%d, nu=2^%d => NRMS = %6.3fx10^-4',log2(mu(ii)),log2(nu(ii)),...
		nrms(xres(:,:,ii),xtrue)*1e4);
end
disp(' ');


% plot best cases
% ind_m = [1:3];
% ind_n = [1:3];
% ind = sub2ind([f.nmu f.nnu],ind_m,ind_n);
% figure(1), plot(0:f.niter,err(:,ind),...
% 	0:f.niter,nrms_dblim*ones(1,f.niter+1),'b--')

