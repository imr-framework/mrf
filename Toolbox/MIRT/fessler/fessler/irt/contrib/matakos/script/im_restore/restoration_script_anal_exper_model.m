%% Initialization parameters

% For BSNR30, TV reg: Smoothed image at lam=2^-7
% For BSNR30, Wavelet reg: Smoothed image at lam=2^-7
% For BSNR40, TV reg: Smoothed image at lam=2^-8
% For BSNR40, Wavelet reg: Smoothed image at lam=2^-8
% For BSNR50, TV reg: Smoothed image at lam=2^-8
% For BSNR50, Wavelet reg: Smoothed image at lam=2^-8

f.do_ncg = 0;
f.do_ista = 0;
f.do_fista = 0;
f.do_mfista = 0;
f.do_salsa = 0;
f.do_reeves = 1;
f.do_alp1 = 3;
f.do_alp2 = 0;
f.do_cgq = 0;
f.do_alq = 0;

f.ne = 0; % experiment counter
% # of experiments to run
f.exper = f.do_ncg + f.do_ista + f.do_fista + f.do_mfista + f.do_salsa ...
	+ f.do_alp1 + f.do_alp2 + f.do_cgq + f.do_alq;

f.etype = cell(f.exper,1);

f.fname1 = 'img_anal_tviso_bsnr40_test_ex.png';
f.fname2 = 'img_anal_tviso_bsnr40_test_ap.png';
f.fname3 = 'img_anal_tviso_bsnr40_test_aps.png';
f.fname4 = 'img_anal_tviso_bsnr40_test_ref.png';
f.bname = 'img_blur_bsnr40_test.png';
f.bnamec = '../results/fig/img_blur_circ_bsnr40.png';

f.snr = 50; % Target BSNR
f.lam = 2^-17*ones(f.exper,1); % regularization parameter
% f.lam = [2^-12 2^-15 2^-17 2^-14 2^-13 2^-8 2^-12];
f.reg_type = 'tv,iso';
% f.model_type = {'circ','circ-edge','circ-rpad','circ-rpadedge','circ-spad',...
% 	'circ-spadedge','circ-ppad','circ-ppadedge','circ-zpad','circ-zpadedge'};
% f.model_type = {'circ','circ-edge','circ-zpadedge','circ-rpadedge','circ-spadedge','non-circ'};
% f.model_type = {'non-circ','non-circ','non-circ','non-circ','non-circ'};
f.model_type = {'non-circ','circ','circ-rpadedge'};

% f.model_type = {'circ'};
f.gen_type = 'non-circ';
f.chiter = 4*ones(f.exper,1);
% f.chiter = [1 1 1 1 4 10 1 4 10 1 4 10 1 4 10 1 4 10 1 1 1 1];
f.cgiter = 10*ones(f.exper,1);
% f.cgiter = [1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 4 4 4 1 2 4 1];
f.lniter = 4*ones(f.exper,1);
% f.lniter = [2 5 ones(1,f.exper)];
f.mu = 2^-12*ones(f.exper,1);
% f.mu = 2.^(-6:2:2);
% f.mu = [2^-6 2^-6];
f.nu = 2^-6*ones(f.exper,1);
% f.nu = [2^-2 2^-3 2^-4 2^-5 2^-6];

f.niter = 500*ones(f.exper,1); % # of iterations
f.maxiter = max(f.niter);
f.plot = 0; % plot results
f.chat = 1; % show printouts


%% Create blur, image, and masks
f.nb = 9; % Size of blur
% f.psf = ones(f.nb)/f.nb^2; % Uniform blur
f.psf = fspecial('gaussian',15,3); % Gaussian blur
% f.psf = padarray(fspecial('motion',15,30),[3 0],'both'); % Motion blur

% True image, Cameraman scaled at [0,1]
xtrue = double(imread('cameraman.tif')).';
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
	f.mask_res = true(nxr,nyr);
	
	nx = nxr;
	ny = nyr;
	f.mask = true(nx,ny);

	xtrue = xtrue((f.nb+1)/2:end-(f.nb-1)/2,(f.nb+1)/2:end-(f.nb-1)/2);
end


%% Create System matrix, noise and data
Ac = Gblur(true(nx,ny), 'psf', f.psf, 'type', 'fft,same');
T = Gdiag2(ones(nxr,nyr),'mask_in',f.mask);

rng(0);
ns = randn(nxr,nyr);

yb = T*(Ac*xtrue); % Blurred noiseless data

f.sig_pow = sum(abs(yb(:)-mean(yb(:))).^2)/numel(yb);
f.sigma = sqrt(f.sig_pow/(10^(f.snr/10)));

y = yb + f.sigma*ns;

xtrue = T*xtrue;
f.climx = [0 max(xtrue(:))+0.1];


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

if f.do_reeves > 0
	for ii=1:f.do_reeves
		f.etype{f.ne+ii} = 'reeves';
	end
	f.ne = f.ne+f.do_reeves;
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

if f.do_cgq > 0
	for ii=1:f.do_cgq
		f.etype{f.ne+ii} = 'cgq';
	end
	f.ne = f.ne+f.do_cgq;
end

if f.do_alq > 0
	for ii=1:f.do_alq
		f.etype{f.ne+ii} = 'alq';
	end
	f.ne = f.ne+f.do_alq;
end

xres = zeros(nxr,nyr,f.exper);
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

for jj=1:f.exper
	
	[xr out] = image_restore(y, f.psf, 'est_type', f.etype{jj}, ...
		'model_type', f.model_type{jj}, 'lambda', f.lam(jj), 'reg_type', f.reg_type, ...
		'xtrue', xtrue,  'mu', f.mu(jj), 'nu', f.nu(jj), 'niter', f.niter(jj), ...
		'chiter', f.chiter(jj), 'cgiter', f.cgiter(jj), 'lniter', f.lniter(jj), ...
		'plot', f.plot, 'chat', f.chat);
	
	if strcmpi(f.model_type{jj}, 'non-circ')
		xr = T*xr;
	end
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
		
% 	printm('Image NRMS in db: %6.4f\n',err(end,jj));
% 	printm('Image ISNR in db: %6.4f\n',isnr(end,jj));
end

figure(1)
im(xres,'Restored images','cbar',f.climx);
% print('-dpng','-r300','-f1',[f.name '_img_full.png'])
% figure(2), im(xres-repmat(xtrue,[1 1 f.exper]),'Difference images','cbar');
% print('-dpng','-r300','-f2',[f.name '_diff_full.png'])
% figure(3), im(wres,'Restored coefficients','cbar',f.climw);
% print('-dpng','-r300','-f3',[f.name '_img_cen.png'])
% figure(4), im((wres-repmat(wtrue,[1 1 1 f.exper])),'Difference coefficients','cbar');
% print('-dpng','-r300','-f4',[f.name '_diff_cen.png'])

figure(2), plot(0:f.maxiter,err),title 'Error per iteration'
figure(3), plot(squeeze(time(:,1,:)),err),title 'Error with CPU time'
figure(4), plot(squeeze(time(:,2,:)),err),title 'Error with Actual time'

disp('NRMSE and ISNR')
disp([nrms(xres, xtrue) err(end,:).' isnr(end,:).']);

disp('Iterations CPU time')
disp([f.cputime_iter f.cputime_cgiter f.cputime_chiter f.cputime_lniter]*1000) 

disp('Iterations Actual time')
disp([f.time_iter f.time_cgiter f.time_chiter f.time_lniter]*1000)

imwrite(xres(:,:,1).',f.fname1,'png');
imwrite(xres(:,:,2).',f.fname2,'png');
imwrite(xres(:,:,3).',f.fname3,'png');
imwrite(y.',f.bname,'png');
% 
% imwrite(xres(:,:,3).','img_anal_tviso_bsnr40_aps.png','png');
% imwrite(y.',f.bnamec,'png');

% save(f.fname, 'xres', 'isnr', 'err', 'time', 'f');
