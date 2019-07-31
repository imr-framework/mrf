%% Initialization parameters

% Experiments for paper timing comparison
% 1) NCG-2
% 2) NCG-5
% 3) ISTA
% 4) FISTA-1
% 5) FISTA-4
% 6) FISTA-10
% 7) FISTA-20
% 8) MFISTA-1
% 9) MFISTA-4
% 10) MFISTA-10
% 11) MFISTA-20
% 12) SALSA-1-1
% 13) SALSA-1-4
% 14) SALSA-1-10
% 15) SALSA-1-15
% 16) SALSA-4-1
% 17) SALSA-4-4
% 18) SALSA-4-10
% 19) SALSA-4-15
% 20) SALSA-10-1
% 21) SALSA-10-4
% 22) SALSA-10-10
% 23) SALSA-10-15
% 24) ALP1-1
% 25) ALP1-4
% 26) ALP1-10
% 27) ALP2
% 28) BRR-1
% 29) BRR-4
% 30) BRR-10

% In total: 2xNCG, 1xISTA, 3xFISTA, 3xMFISTA, 9xSALSA, 3xALP1, 1xALP2
% Can skip all of FISTA, ALP1-1, SALSA-1-X and SALSA-X-1 for a total of 13 exper
% Use 2500 iterations for all "slow" methods and 6000 for "faster" methods
f.do_ncg = 2;
f.do_ista = 1;
f.do_fista = 4;
f.do_mfista = 4;
f.do_salsa = 12;
f.do_alp1 = 3;
f.do_alp2 = 1;
f.do_reeves = 3;

f.ne = 0; % experiment counter
% # of experiments to run
f.exper = f.do_ncg + f.do_ista + f.do_fista + f.do_mfista + f.do_salsa ...
	+ f.do_alp1 + f.do_alp2 + f.do_reeves;

f.etype = cell(f.exper,1);

f.fname = '../results/mat/res_anal_exact_wave_bsnr40_circ_t.mat';
f.prefix = '../results/mat/ind/res_anal_exact_wave_bsnr40_circ_t';

f.snr = 40; % Target BSNR
f.lam = 2^-15*ones(f.exper,1); % regularization parameter
% f.lam = [2^-10 2^-10.5 2^-11 2^-11.5 2^-12];
f.reg_type = 'tv,aniso';
f.model_type = 'non-circ';
f.gen_type = 'non-circ';
% f.chiter = 1*ones(f.exper,1);
f.chiter = [1 1 1 1 4 10 20 1 4 10 20 1 4 10 15 1 4 10 15 1 4 10 15 1 1 1 1 1 1 1];
% f.cgiter = 1*ones(f.exper,1);
% f.cgiter = [ones(1,15) 1000*ones(1,4) 20*ones(1,4) 1 1000 20 1 1 4 10];
f.cgiter = [ones(1,15) 1000*ones(1,4) 20*ones(1,4) 1 1000 20 1 4 4 4];
% f.lniter = 4*ones(f.exper,1);
f.lniter = [2 5 ones(1,f.exper)];
% f.mu = 2^-4*ones(f.exper,1);
f.mu = [ones(1,11) 2^-8*ones(1,12) 2^-10*ones(1,3) 2^-4 2^-10*ones(1,3)];
f.nu = 2^-6*ones(f.exper,1);
% f.nu = [2^-2 2^-3 2^-4 2^-5 2^-6];
% f.precon = [1 1 0 0 0];
% f.epsilon = [1e-6 1e-4 1e-6 1e-4 0];

% f.niter = 1*ones(f.exper,1); % # of iterations
% f.niter = [10*ones(1,5) 20*ones(1,10) 15*ones(1,10)];
% # iter for xinf (TV)
f.niter = [3000 3000 15000 10000 8000 6000 4000 10000 8000 6000 4000 ...
	400*ones(1,15) 5000 400*ones(1,3)];
f.niter = 400*ones(1,30);
% # iter for xtrue (TV)
% f.niter = [200 200 2000 600 500 400 300 600 500 400 300 400 300 200 ...
% 	100 400 300 200 100 400 300 200 100 300 300 300 1000];
f.cgr_thres = [ones(1,27) 0 2 4];
f.maxiter = max(f.niter);
f.plot = 0; % plot results
f.chat = 1; % show printouts


%% Create blur, image, and masks
f.nb = 9; % Size of blur
f.psf = ones(f.nb)/f.nb^2; % Uniform blur
% f.psf = fspecial('gaussian', f.nb, 1.5); % Gaussian blur with sigma=2
% f.psf = padarray(fspecial('motion', f.nb, 30), [2 0]); % Gaussian blur with sigma=1.5

% True image, Cameraman scaled at [0,1]
xtrue = double(imread('../data/fig/house.tif')).';
xtrue = xtrue/max(xtrue(:));

% xtrue = zeros(32);
% xtrue(9:24,9:24) = 1;

[nx ny] = size(xtrue);

if strcmpi(f.gen_type, 'non-circ')
	f.mask = false(nx,ny);
	f.mask((f.nb+1)/2:end-(f.nb-1)/2,(f.nb+1)/2:end-(f.nb-1)/2) = true;
	
	nxr = nx-f.nb+1;
	nyr = ny-f.nb+1;
	f.mask_res = true(nxr,nyr);
else
	nxr = nx-f.nb+1;
	nyr = ny-f.nb+1;
	f.mask_res = true(nx,ny);
	
% 	nx = nxr;
% 	ny = nyr;
	f.mask = true(nx,ny);

% 	xtrue = xtrue((f.nb+1)/2:end-(f.nb-1)/2,(f.nb+1)/2:end-(f.nb-1)/2);
end


%% Create System matrix and data
Ac = Gblur(true(nx,ny), 'psf', f.psf, 'type', 'fft,same');
T = Gdiag2(ones(nxr,nyr),'mask_in',f.mask);

rng(0);
ns = randn(nxr,nyr);
% ns = randn(nx,ny);

yb = T*(Ac*xtrue); % Blurred noiseless data
% yb = Ac*xtrue;

f.sig_pow = sum(abs(yb(:)-mean(yb(:))).^2)/numel(yb);
f.sigma = sqrt(f.sig_pow/(10^(f.snr/10)));

y = yb + f.sigma*ns;

xtrue = T*xtrue;
f.climx = [0 max(xtrue(:))+0.1];

load ../data/mat/xinf_anal_exact_tv_bsnr40.mat

xtrue = mean(xres,3);
xtrue = xres(:,:,2);
f.climx = [0 max(xtrue(:))+0.1];

nrms_lim = nrms(xres(:,:,1),xtrue);
nrms_dblim = 20*log10(nrms_lim);

printm('NRMSE bound =%6.3fx10^-4, in db =%6.2f\n', nrms_lim*1e4, nrms_dblim);

%% Setup and run experiments

if f.do_ncg > 0
	for ii=1:f.do_ncg
		f.etype{f.ne+ii} = 'ncg';
	end
	f.ne = f.ne+f.do_ncg;
end

if f.do_ista > 0
	for ii=1:f.do_ista
		f.etype{f.ne+ii} = 'ista';
	end
	f.ne = f.ne+f.do_ista;
end

if f.do_fista > 0
	for ii=1:f.do_fista
		f.etype{f.ne+ii} = 'fista';
	end
	f.ne = f.ne+f.do_fista;
end

if f.do_mfista > 0
	for ii=1:f.do_mfista
		f.etype{f.ne+ii} = 'mfista';
	end
	f.ne = f.ne+f.do_mfista;
end

if f.do_salsa > 0
	for ii=1:f.do_salsa
		f.etype{f.ne+ii} = 'salsa';
	end
	f.ne = f.ne+f.do_salsa;
end

if f.do_alp1 > 0
	for ii=1:f.do_alp1
		f.etype{f.ne+ii} = 'alp1';
	end
	f.ne = f.ne+f.do_alp1;
end

if f.do_alp2 > 0
	for ii=1:f.do_alp2
		f.etype{f.ne+ii} = 'alp2';
	end
	f.ne = f.ne+f.do_alp2;
end

if f.do_reeves > 0
	for ii=1:f.do_reeves
		f.etype{f.ne+ii} = 'reeves';
	end
	f.ne = f.ne+f.do_reeves;
end

xres = zeros(nxr,nyr,f.exper);
% xres = zeros(nx,ny,f.exper);
isnr = zeros(f.maxiter+1,f.exper);
err = zeros(f.maxiter+1,f.exper);
time = zeros(f.maxiter+1,2,f.exper);
f.cputime_iter = zeros(f.exper,1);
f.cputime_cgiter = zeros(f.exper,1);
f.cputime_chiter = zeros(f.exper,1);
f.cputime_lniter = zeros(f.exper,1);
f.time_iter = zeros(f.exper,1);
f.time_cgiter = zeros(f.exper,1);
f.time_chiter = zeros(f.exper,1);
f.time_lniter = zeros(f.exper,1);

% f.exper = [16:18 25 27 2 3 8:10];
% f.exper = [24 25 28 29];
f.exper = [27 28 29 30];

for jj=f.exper
		
	[xr out] = image_restore(y, f.psf, 'est_type', f.etype{jj}, ...
		'model_type', f.model_type, 'lambda', f.lam(jj), 'reg_type', f.reg_type, ...
		'xtrue', xtrue,  'mu', f.mu(jj), 'nu', f.nu(jj), 'niter', f.niter(jj), ...
		'chiter', f.chiter(jj), 'cgiter', f.cgiter(jj), 'lniter', f.lniter(jj), ...
		'cgr_thres',f.cgr_thres(jj),'plot', f.plot, 'chat', f.chat);
	
	xr = embed(xr(f.mask),f.mask_res);
	xres(:,:,jj) = xr;
	isnr(1:f.niter(jj)+1,jj) = out.isnr;
	err(1:f.niter(jj)+1,jj) = out.err;
	time(1:f.niter(jj)+1,:,jj) = out.time;
	f.cputime_iter(jj) = out.cputime_iter;
	f.cputime_cgiter(jj) = out.cputime_cgiter;
	f.cputime_chiter(jj) = out.cputime_chiter;
	f.cputime_lniter(jj) = out.cputime_lniter;
	f.time_iter(jj) = out.time_iter;
	f.time_cgiter(jj) = out.time_cgiter;
	f.time_chiter(jj) = out.time_chiter;
	f.time_lniter(jj) = out.time_lniter;
	
	if f.niter(jj) < f.maxiter
		isnr(f.niter(jj)+2:end,jj) = isnr(f.niter(jj)+1,jj);
		err(f.niter(jj)+2:end,jj) = err(f.niter(jj)+1,jj);
		time(f.niter(jj)+2:end,1,jj) = time(f.niter(jj)+1,1,jj);
		time(f.niter(jj)+2:end,2,jj) = time(f.niter(jj)+1,2,jj);
	end
	
	figure(1), im(xtrue,'Original image','cbar',f.climx), drawnow;
	figure(2), im(y,'Blured noisy image','cbar'), drawnow;
	figure(3), im(xr,'Restored image','cbar',f.climx), drawnow;	
	figure(4), im(xr-xtrue,'Difference image','cbar'), drawnow;
	
	save([f.prefix num2str(jj)], 'xr', 'out');
	
% 	printm('Image NRMS in db: %6.4f\n',err(end,jj));
% 	printm('Image ISNR in db: %6.4f\n',isnr(end,jj));
end

figure(1), im clf
im(xres(:,:,f.exper),'Restored images','cbar',f.climx);
% print('-dpng','-r300','-f1',[f.name '_img_full.png'])
% figure(2), im(xres-repmat(xtrue,[1 1 f.exper]),'Difference images','cbar');
% print('-dpng','-r300','-f2',[f.name '_diff_full.png'])
% figure(3), im(wres,'Restored coefficients','cbar',f.climw);
% print('-dpng','-r300','-f3',[f.name '_img_cen.png'])
% figure(4), im((wres-repmat(wtrue,[1 1 1 f.exper])),'Difference coefficients','cbar');
% print('-dpng','-r300','-f4',[f.name '_diff_cen.png'])

figure(2), plot(0:f.maxiter,err(:,f.exper)),title 'Error per iteration'
figure(3), plot(squeeze(time(:,1,f.exper)),err(:,f.exper)),title 'Error with CPU time'
figure(4), plot(squeeze(time(:,2,f.exper)),err(:,f.exper)),title 'Error with Actual time'

disp('NRMSE and ISNR')
disp([nrms(xres, xtrue) err(end,:).' isnr(end,:).']);

disp('Iterations CPU time')
disp([f.cputime_iter/1000 f.cputime_cgiter f.cputime_chiter f.cputime_lniter]*1000) 

disp('Iterations Actual time')
disp([f.time_iter f.time_cgiter f.time_chiter f.time_lniter]*1000)

% imwrite(xres(:,:,1).',f.fname,'png');
% imwrite(y.','../results/fig/bsnr50/img_blur_bsnr50.png','png');

% save(f.fname, 'xres', 'isnr', 'err', 'time', 'f');
