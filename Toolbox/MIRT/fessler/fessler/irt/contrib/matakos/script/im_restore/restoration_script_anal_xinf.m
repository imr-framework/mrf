%% Initialization parameters

% # of experiments
f.nexp = 6;

% File names for experiment results
f.fname{1} = '../data/mat/xinf_anal_exact_tv_bsnr30.mat';
f.fname{2} = '../data/mat/xinf_anal_exact_tv_bsnr40_circ.mat';
f.fname{3} = '../data/mat/xinf_anal_exact_tv_bsnr50.mat';
f.fname{4} = '../data/mat/xinf_anal_exact_wave_bsnr30.mat';
f.fname{5} = '../data/mat/xinf_anal_exact_wave_bsnr40_circ.mat';
f.fname{6} = '../data/mat/xinf_anal_exact_wave_bsnr50.mat';

% BSNR levels for all experiments
f.snr = [30 40 50 30 40 50];

% Regularization type for all experiments
f.reg_type{1} = 'tv,aniso';
f.reg_type{2} = 'tv,aniso';
f.reg_type{3} = 'tv,aniso';
f.reg_type{4} = 'l1,wave';
f.reg_type{5} = 'l1,wave';
f.reg_type{6} = 'l1,wave';

% Lagrangian parameters for all experiments
f.mu = [2^-4 2^-4 2^-4 2^-4 2^-4 2^-4];
f.nu = [2^-4 2^-6 2^-7 2^1 2^-1 2^-3];

% Regularization parameter for all experiments
f.lam = [2^-12 2^-15 2^-17 2^-12 2^-15 2^-17];

% Inner iterations for MFISTA
f.chiter = 5*ones(1,6);
f.cgiter = 4*ones(1,6);

% # of part of experiment
f.exper = 2;
f.etype{1} = 'alp1';
f.etype{2} = 'alp2';


f.niter = 50000; % # of iterations
f.plot = 0; % plot results
f.chat = 1; % show printouts

%% Create blur, image, and masks
f.nb = 9; % Size of blur
f.psf = ones(f.nb)/f.nb^2; % Uniform blur
% f.psf = fspecial('gaussian', f.nb, 1.5); % Gaussian blur with sigma=1.5
% f.psf = padarray(fspecial('motion', f.nb, 30), [2 0]); % Gaussian blur with sigma=1.5

% True image, Cameraman scaled at [0,1]
xtrue = double(imread('cameraman.tif')).';
xtrue = xtrue/max(xtrue(:));

[nx ny] = size(xtrue);

f.mask = true(nx,ny);
% f.mask = false(nx,ny);
% f.mask((f.nb+1)/2:end-(f.nb-1)/2,(f.nb+1)/2:end-(f.nb-1)/2) = true;

nxr = nx-f.nb+1;
nyr = ny-f.nb+1;
% f.mask_res = true(nxr,nyr);
f.mask_res = true(nx,ny);


%% Create System matrix and data
Ac = Gblur(true(nx,ny), 'psf', f.psf, 'type', 'fft,same');
% T = Gdiag2(ones(nxr,nyr),'mask_in',f.mask);

reset(RandStream.getDefaultStream);
% ns = randn(nxr,nyr);
ns = randn(nx,ny);

% yb = T*(Ac*xtrue); % Blurred noiseless data
yb = Ac*xtrue;

f.sig_pow = sum(abs(yb(:)-mean(yb(:))).^2)/numel(yb);

% xtrue = T*xtrue; % True image cropped to valid size of convolution
f.climx = [0 max(xtrue(:))+0.1]; % Magnitude display limits


%% Setup and run experiments

for ii=[2 5]
	
	% Add noise to data
	f.sigma = sqrt(f.sig_pow/(10^(f.snr(ii)/10)));
	y = yb + f.sigma*ns;

% 	xres = zeros(nxr,nyr,f.exper);
	xres = zeros(nx,ny,f.exper);
	isnr = zeros(f.niter+1,f.exper);
	err = zeros(f.niter+1,f.exper);
	time = zeros(f.niter+1,2,f.exper);
	
	for jj=1:f.exper
		
		[xr out] = image_restore(y, f.psf, 'est_type', f.etype{jj}, ...
			'model_type', 'circ', 'lambda', f.lam(ii), ...
			'reg_type', f.reg_type{ii}, 'mu', f.mu(ii), 'nu', f.nu(ii), ...
			'xtrue', xtrue, 'niter', f.niter, 'chiter', f.chiter(ii), ...
			'cgiter', f.cgiter(ii), 'plot', f.plot, 'chat', f.chat);
		
		xr = embed(xr(f.mask),f.mask_res);
		
		xres(:,:,jj) = xr;
		isnr(:,jj) = out.isnr;
		err(:,jj) = out.err;
		time(:,:,jj) = out.time;
		
		figure(1), im(xtrue,'Original image','cbar',f.climx), drawnow;
		figure(2), im(y,'Blured noisy image','cbar'), drawnow;
		figure(3), im(xr,'Restored image','cbar',f.climx), drawnow;
		figure(4), im(xr-xtrue,'Difference image','cbar'), drawnow;
		
% 		save('../data/mat/xinf_anal_exact_tv_bsnr40_gauss9_al.mat', 'xres', 'err', 'time', 'y');
		
% 		printm('Image NRMS in db: %6.4f\n',err(end,jj));
	end
	
	figure(1), im(xres,'Restored images','cbar',f.climx);
% 	figure(2), im(xres-repmat(xtrue,[1 1 f.exper]),'Difference images','cbar');
	figure(2), plot(1:length(err),err),title 'Error per iteration'
	figure(3), plot(squeeze(time(:,1,:)),err),title 'Error with CPU time'
	figure(4), plot(squeeze(time(:,2,:)),err),title 'Error with Actual time'
	
	disp([nrms(xres, xtrue) err(end,:).' isnr(end,:).']);
	
	save(f.fname{ii}, 'xres', 'err', 'time', 'y');
	
end



