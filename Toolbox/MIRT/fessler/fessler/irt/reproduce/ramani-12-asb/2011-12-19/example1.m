%| example1.m
%| CT 2D Recon Example - simulation using a numerical 2D NCAT phantom
%|
%| Comparison of
%|
%| MFISTA (Beck et al, IEEE TIP, vol. 18, no. 11, 2009)
%|
%| and
%|
%| ADMM (Ramani et al, IEEE TMI, doi: 10.1109/TMI.2011.2175233, in press)
%|
%| algorithm for statistical 2-D CT reconstruction
%|
%| 2011-11-30, Sathish Ramani, University of Michigan
%| 2015-06-17, J Fessler, cosmetic changes, derived from "Example.m"

%% ensure path includes IRT
if isempty(strfind(path, '/fbp'))
	error('irt not in path')
end


%% Image geometry for recon image
if ~isvar('ig')
	fov = 65; % in cms
	down = 8; % downsampling factor for recon image
	nxf = 1024; % NCAT size
	nyf = 1024; % NCAT size

	nx = nxf / down;
	ny = nyf / down;

	N = [nx ny];
	Npix = prod(N); % # of pixels in image

	ig = image_geom('nx', nxf, 'fov', fov, 'down', down); % image geometry structure
	ig_full = ig; % no masking
	ig.mask = logical(ig.circ(fov/2.5, fov/3, 0, 0)); % Mask: image outside mask = 0; also helps speed up computation of Ax
	mask = ig.mask;
	Nmask = length(find(mask>0)); % # of masked samples

	Nc = floor(N / 2) + 1; % Center of image
	idmask = find(mask>0);
	[zr zc] = ind2sub(N, idmask);
	distx = max([abs(min(zr(:)) - Nc(1)) abs(max(zr(:)) - Nc(1))]);
	disty = max([abs(min(zc(:)) - Nc(2)) abs(max(zc(:)) - Nc(2))]);

	zr = Nc(1) - distx : Nc(1) + distx - 1; % Range of pixels to display
	zc = Nc(2) - disty : Nc(2) + disty - 1; % Range of pixels to display
end


%% Sinogram geometry structure (fan-beam)
if ~isvar('sg')
	sg = sino_geom('ge1', 'units', 'cm', 'down', 1, 'na', 984); % Jeff's sinogram geometry structure, fan-beam geometry
	Ndata = sg.nb * sg.na; % # of data samples
	sg.plot(ig);
end


%% Use high res phantom to generate sinogram and reconstruct image on a downsampled grid
% Read noisefree test image and generate more realistic noisefree sinogram from finer image
if ~isvar('data_true')
	imgBig = read_ncat('nx', nxf, 'ny', nyf, 'marrow', true, ...
		'mu', [0 0.05 0.2 0.4 0.4 0.2]/0.2*1000); % True highres NCAT
	imgBig = single(imgBig) / single(max(imgBig(:))) * 0.4; % convert to 1/cm units
	im(imgBig)

	% noisefree phantom on recon grid for RMSE purposes
	img = downsample2(imgBig, down);
	mn = min(img(:));
	mx = max(img(:));
	im(img)
	jf_equal(img, mask .* img) % ensure object within mask!

	if 0 % no stored data
		try
			load('data_true.mat', 'data_true');
			if(~exist('data_true', 'var'))
			 	rundata = 1; else rundata = 0; end
		catch ME
			rundata = 1;
		end
	end

	% image geometry for generating data (no downsampling)
	ig_big = image_geom('nx', nxf, 'fov', fov, 'down', 1);
	Abig = Gtomo2_dscmex(sg, ig_big); % System matrix corresponding to highres phantom

	printm('Generating sinogram from high-res phantom for simulation...')
	cpu etic
	data_true = Abig * imgBig; % use highres phantom for A*x
	cpu etoc
	im(data_true)
%	save data_true.mat data_true;
%	clear ig_big
	clear Abig
end


%% Data Generation - Noisy sino data
if ~isvar('yi')
	I0 = 1e5; % incident photons; decrease this for "low dose" scans
	poissonfactor = 0.4; % for generating poissonn noise using rejection method
	rng(0)
%	seedr = 0; rand('state', seedr);

	yi = poisson(I0 * exp(-data_true), 0, 'factor', poissonfactor); % poisson noise for transmission data:

	if any(yi(:) == 0)
		warn('%d of %d values are 0 in sinogram!', sum(yi(:)==0), length(yi(:)));
	end
	data = log(I0 ./ max(yi,1)); % noisy fan-beam sinogram: form of the data in the quadratic approximation to the actual log-likelihood

	im plc 2 2
	im(1, data_true)
	im(2, data)
	im(3, data - data_true)
	im(4, yi)
end


%% Forward system model for recon as a fatrix object
if ~isvar('wi')
	A = Gtomo2_dscmex(sg, ig); % System matrix used for recon
	Af = Gtomo2_dscmex(sg, ig_full); % System matrix used for generating FFT-preconditioner in ADMM

	%% PWLS weights
	wi = yi / max(yi(:)); % PWLS weights; gives 0 weight to any ray where yi=0!
	mxW = max(wi(:));
	mnW = min(wi(:));

	% Crop weights to avoid "air" regions in raw data
	wistack = sort( wi(:) );
	wistack = wistack( wistack <= 0.9 );
	wistack = wistack( wistack > 0 ); % only nonzero weights less than 0.9 * max are considered for selecting the AL parameters
	medwi = median( wistack(:) ); % Median of weights for use in ADMM
	clear wistack

	im(1, wi)
	pr medwi

	Wy = wi .* data;
	AWy = A' * Wy; % (weighted) back-projection
	im(2, AWy)
end


%% Kappa for space-varying weights
if ~isvar('kappa'), printm('kappas')
	kappa = sqrt( div0(A' * wi, A' * ones(size(wi))) );
	skap = sort(kappa(:));
	skap = skap(skap > 0);
	kappa(kappa==0) = skap(1); % Ensure kappa does not have any zeros
	rw = kappa .^ 2; % spatially varying regularization weights; usually kappa ^ 2, but can be adjusted for quality
	im(kappa)
end


%% FBP Recon
if ~isvar('fbpnoisy'), printm('fbp')
	tmp = fbp2(sg, ig); % Jeff's FBP function

	tic;
	fbpnoisy = fbp2(data, tmp); % FBP with ramp filter
	TFBP = toc;
	errfbp = mask .* ( fbpnoisy - img );
	RMSEFBP = sqrt( sum( abs( errfbp(:) ).^2 ) / Nmask );

	xini = fbpnoisy; % FBP initialization for iterative recon
	xini = xini .* mask;
	errini = mask .* ( xini - img );
	errini = sqrt( sum( abs( errini(:) ).^2 ) / Nmask );

	scale = 1; % adjust this to shift units of recon to Hounsefield
	clim = scale * [mn mx]; % Display limits
	im(xini, clim)
end


%% Params structure used globally in all algos
if ~isvar('params'), printm('params')
	params.zoomr = zr;
	params.zoomc = zc;
	params.clim = clim;
	params.scale = scale;
	params.fov = fov;
	params.ig = ig;
	params.sg = sg;
	params.img = img;
	params.rs = nx;
	params.cs = ny;
	params.mn = mn;
	params.mx = mx;
	params.N = N;
	params.Npix = Npix;
	params.Ndata = Ndata;
	params.Nmask = Nmask;

	params.kappa = kappa;
	params.rw = rw;
	params.fbpnoisy = fbpnoisy;
	params.xini = xini;
	params.AWy = AWy;
	params.Wy = Wy;

	params.A = A;
	params.W = wi;
	params.mxW = mxW;
	params.mnW = mnW;

	% Regularization parameter for statistical recon
	lambda = medwi * 20; % heuristic!
	params.lambda = lambda;

	% Parameters for Preconditioned CG-solver inside ADMM for "inverting" A'A + nu * R'R
	params.CG.precon = 1; % Do Precondition CG-solver inverting (A'A + mu*R'R) in ADMM
	params.CG.restol = 1e-8; % Tol for residue
	params.CG.nCG = 2; % # of CGs
	params.MFISTA.nCG = 5; % maximum # of CG iterations inside MFISTA

	% x_infinity solution
	params.xinf = zeros(nx,ny); % to be populated when x_infinity ( minimizer of the cost ) is available
	params.xinfnorm = 1;

	% Parameters for cost decrease, etc
	params.dcosttol = 1e18;
	params.dxtol = 0;

	% Display parameters
	params.dispfig = 0; % turn on to display recon at every iteration
	params.dispitr = 1; % turn on to display cost, etc., at every iteration
	params.dispitrnum = 3; % display info once these many iterations

	% ADMM parameters
	params.AL.mu = 1; % To be populated later
	params.AL.nu1 = 1; % To be populated later
	params.AL.nu1AL1factor = 100; % divide nu in (CAA + nu * R'R) by this factor

	params.AL.iWmu = [];
	params.AL.iRnu2nu1 = [];

	params.formatstringC = '%0.5E';
	params.formatstringT = '%0.3E';
	params.formatstringO = '%0.2E';

	params.Constants.EPSL = 1.0E-10; % Some constants needed for Brent iterations for minimization without derivatives; used for condition number minimization
	params.Constants.CGOLD = 0.3819660;
	params.Constants.GOLD = 1.618034;
	params.Constants.GLIMIT = 100.0;
	params.Constants.Brenttol = 0.01;
	params.Constants.maxBitr = 100;

	%% Recon Setup
	% Potential Function
	params.PriorType = 'l1'; % Absolute function = |t|
	% params.PriorType = 'FP'; % Fair potential = alpha^2 * ( |t| / alpha - log(1 + |t| / alpha ) )
	% params.Prior.alpha = 1e-4; % smoothing parameter for Fair potential

	% Type of reg op
	% params.Operator = 'AFD'; % plain finite differences (bilateral form)
	% params.Operator = 'FD'; % gradient-norm of finite differences; **** for TV use l1-FD ****
	params.Operator = 'W'; % Undecimanted (shift-invariant) wavelet transform

	% Wavelet Options
	params.dwtmode = 'per'; % Period boundaries for wavelet implementation
	params.redundancy = 'undecimated'; % Undecimated wavelet transform using wavelet filters corresponding to standard orthonormal wavelets
	params.wname = 'haar'; % wavelet name, see waveinfo.m
	params.nlev = 2; % # of levels in wavelet transform
	params.includeApprox = false; % do not include wavelet approximation level in the reg op

	% Wavelet filters
	[lod, hid, lor, hir] = wfilters(params.wname);

	% Normalize filters so as to avoid a product with 0.5 during inverse undecimated wavelet transform
	params.lod = lod/sqrt(2);
	params.hid = hid/sqrt(2);
	params.lor = lor/sqrt(2);
	params.hir = hir/sqrt(2);
end


% Regularization object
if ~isvar('params.R'), printm('R')
	switch params.Operator
	case {'W'}
		% With mask supporting only the object, used in the recon algo
		R = Gwave2('mask', mask, 'dwtmode', 'per', ...
			'redundancy', 'undecimated', 'wname', 'haar', ...
			'nlevel', params.nlev, 'includeApprox', false, ...
			'scale_filters', 1/sqrt(2));

		% Without mask, used only to obtain the freq.resp of Rf'Rf for preconditioner
		Rf = Gwave2('mask', true(N), 'dwtmode', 'per', ...
			'redundancy', 'undecimated', 'wname', 'haar', ...
			'nlevel', params.nlev, 'includeApprox', false, ...
			'scale_filters', 1/sqrt(2));

	case {'FD', 'AFD'}
		fail('Fatrix for finite differences should be inserted here...')

	otherwise
		fail('Wrong choice for regularization operator')
	end

	params.R = R;
end


%% preconditioner related
if ~isvar('params.mEAWA')
	ei = zeros(nx,ny); ei(nx/2+1,ny/2+1) = 1; RRei = Rf' * ( Rf * ei );
	RR = abs(fft2(RRei)); % freq. response corresponding to circulant R'R (uses Rf without masks)
	mxRR = max(RR(:));
	mnRR = min(RR(:));

	params.RR = RR;
	params.mxRR = mxRR;
	params.mnRR = mnRR;

	%% Calculate (or load) max eigenvalue of A'WA for MFISTA
	fname_EvalsAWA = ['Z_Nx' int2str(nx) '_FOV' num2str(fov) '_DS' int2str(down) '_Sino' int2str(sg.nb) 'x' int2str(sg.na) '_Int' num2str(I0) '_mEAWA.mat'];
	params.fname_EvalsAWA = fname_EvalsAWA;

	% Power method parameters; Power method is used to find the max.eigenvalue of a matrix
	if exist(fname_EvalsAWA, 'file')
		load(fname_EvalsAWA, 'mEAWA');
	else
		warn 'no file?'
		params.eigtol = eps; % Matlab epsilon
		params.eigpowermaxitr = 10000;
		params.eigdispitr = 10;

		printm('Computing max eigvalue of AWA using Power method for MFISTA...');
		tic
		mEAWA = get_MaxEigenValue_1by1(params, 'AWA'); % Find maximum eigenvalues of the forward system and the regularization
		toc
		printm(['Max eigvalue of AWA =' num2str(mEAWA)]);
		% save(fname_EvalsAWA, 'mEAWA');
	end

	params.mEAWA = mEAWA;
	params.eiginflate = 1.0025;
	params.mEval_adjust = mEAWA * params.eiginflate;
end


%% Circulant approximation to A'A by taking response of A'A to an impulse at the center
if ~isvar('CAA'), printm('ADMM: Minimizing the condition number of CAA + nu1*RR');
	ei = zeros(nx,ny); ei(nx/2+1,ny/2+1) = 1; AAei = Af' * Af * ei;

	CAA = abs( real( fft2( fftshift( AAei ) ) ) ); % freq. resp corresponding to the circulant approximation to A'A
	params.CAA = CAA;
	params.BRbx = 1e-2; %% Center of first bracket for Brent optimization of condition number

	[minCondCAARRapprox, nuoptapprox] = Brent_Linmin(params); %% Input arguments unused except for "params"; Need to simplify code using varargin!

	params.nuoptapprox = nuoptapprox; % approximate nu that minimizes condition number of original system (A'A + nu * R'R)
	params.CondCAARRapprox = minCondCAARRapprox; % approximate min.cond.num
	printm(['ADMM: Approx min cond num = ' num2str(minCondCAARRapprox, '%0.4E'), '; approx optimal nu = ' num2str(nuoptapprox, '%0.4E')]);
end

%% Compare algos
% Run MFISTA
if ~isvar('xMFIS'), printm('Doing MFISTA ...')
	params.figno = 2;
	params.doMFISTA = 1;
	params.maxitr = 39; % Max # of outer iterations

	[xMFIS CMFIS TFIS l2DFIS RMSEFIS] = runMFISTA(data, AWy, xini, params);
	TFIS = cumsum(TFIS);
	im(xMFIS)
end


% Run ADMM
if ~isvar('xADMM'), printm('Doing ADMM ...')
	params.figno = 3;
	params.AL.mu = medwi;
	params.AL.nu1 = nuoptapprox / params.AL.nu1AL1factor;
	kWAL1 = (mxW + params.AL.mu)/(mnW + params.AL.mu);
	params.AL.iWmu = 1 ./ (wi + params.AL.mu);
	params.maxitr = 18; % Max # of outer iterations

	printm(['Doing ADMM with mu = ' num2str(params.AL.mu) '; nu1 = ' num2str(params.AL.nu1) '...'])
	[xADMM CADMM TADMM l2DADMM RMSEADMM] = runADMM(data, xini, params);
	TADMM = cumsum(TADMM);
	im(xADMM)
end


%% Plot results and show recon
if 1
	im plc 2 3
	subplot(133)
	plot(TFIS, RMSEFIS, 'b', TADMM, RMSEADMM, 'r')
	xlabel('Time (seconds)')
	ylabel('RMSE in cm^{-1}')
	legend(['MFISTA - ' int2str(params.MFISTA.nCG) ' inner CG itrs.'], ['ADMM - ' int2str(params.CG.nCG) ' inner CG itrs.'])
	title('RMSE versus algorithm run-time')

	im(1, scale * img(zr,zc), clim), cbar
	title('Noisefree')

	im(2, scale * fbpnoisy(zr, zc), clim), cbar
	title(['FBP recon after ' num2str(TFBP, '%0.1f') ' seconds'])
	xlabel(['RMSE = ' num2str(RMSEFBP(end), '%0.4f') ' cm^{-1}'])

	im(4, scale * xMFIS(zr, zc), clim), cbar
	title(['MFISTA recon after ' num2str(TFIS(end), '%0.1f') ' seconds'])
	xlabel(['RMSE = ' num2str(RMSEFIS(end), '%0.4f') ' cm^{-1}'])

	im(5, scale * xADMM(zr, zc), clim), cbar
	title(['ADMM recon after ' num2str(TADMM(end), '%0.1f') ' seconds'])
	xlabel(['RMSE = ' num2str(RMSEADMM(end), '%0.4f') ' cm^{-1}'])
end
