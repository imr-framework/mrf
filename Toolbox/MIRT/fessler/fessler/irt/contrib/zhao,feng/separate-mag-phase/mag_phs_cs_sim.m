% mag_phs_cs_sim.m
%
% This script is an example of MRI reconstruction with separate magnitude and
% phase regularization via Compressed Sensing (CS) for simulation experiments.
% 2D MRI thermometry data is simulated based on data from "ISMRM Reconstruction
% Challenge 2010". For details of the data, please refer to our paper.
%
% Copyright 2012-06-15, Feng Zhao, University of Michigan
% 2013-04-11 modified to use fatrix and to enable octave by JF

%% Load data
if ~isvar('xt')
	more off
	load abd_therm_sim
	mt = abs(im2_l); % true magnitude image
	xt = xt2; % true phase image
end

%% Setup reconstruction parameters
if ~isvar('clim')
	Nd = size(mask1); % high res image size
	Nd1 = size(mt); % low res image size
	curv = prod(Nd1); % spectral radius of A'*A

	sig = 30/1.41; % control the Gaussian noise level
	sr = 0.38; % sampling rate without the fully sampled k-space center
	calirt = 0.03; % ratio of sampling for the fully sampled k-space center
	edgepr = 1; % use edge preserving smoothing regularizer or not
	mask = msk2; % loose low-res mask

	piter1 = 15; % # of iterations step 3 of initialization
	piter2 = 65; % # of iterations for updating phase
	nsubiter = [4 1]; % # of subiterations in each iteration
	beta1 = 1e4; % regularization parameter for phase (rg1/rg3)
	typ_init = 1; % 0: initialize by rg3; 1: initialize by rg2
	beta11 = 7e6; % regularization parameter for phase (rg2/rg4)
	beta22 = 0.5*2^-4 * curv; % regularization parameter for magnitude
	if edgepr == 1;
		del = 0.0005; % parameter for edge-preserving
	else
		piter2 = 0;
	end
	preiter = 2; % # of iterations for intialization by conventional CS
end

%% Setup system objects
if ~isvar('A')
%	rand('state', 1)
	rng(0) % jf - may differ from seed used in paper
	samp = rand(Nd1) > 1 - sr;
	rc = calirt;
	i1 = round(Nd1(1)*(1-sqrt(rc))/2);
	i2 = round(Nd1(2)*(1-sqrt(rc))/2);
	samp(i1+1:Nd1(1)-i1,i2+1:Nd1(2)-i2) = true; % random sampling for the low res data
	unsamp = sum(samp(:))/prod(Nd1); % calculate the actual sampling rate;

	samp2 = logical(embed(samp(:),samp1)); % sampling pattern on the high res data

	% system fatrix for image reconstruction
	A = Gdft('mask', mask1, 'samp', samp2, 'ifftshift', 1, 'fftshift', 1);
	% system fatrix for data simulation
	A1 = Gdft('mask', mask, 'samp', samp, 'ifftshift', 1, 'fftshift', 1);


	%% Setup the wavelet transform object
%	U = Gwave1('mask', mask, 'noCA', 0, 'level', 3); % FZ
	U = Godwt1(mask, 'level', 3); % jf, for octave
	U = U';
end


%% Simulate the MRI data acquistion
if ~isvar('yi')
	rng(0) % jf, this may differ from seed used in paper
	rt = prod(Nd)/prod(Nd1); % ratio between pixel # of high-res image and low-res image
	yb = A * img(mask1)/rt;
	yi = yb + sig * (randn(size(yb)) + 1i * randn(size(yb))); % add complex Gaussian noise
	snr = 20 * log10( norm(yb(:)) / norm(yi(:) - yb(:)) ); pr snr
end


%% The 1st step initialization by inverse FFT
if ~isvar('xi')
	img_i = A1' * yi / prod(Nd1);
	mi = abs(img_i);
	xi = angle(img_i);
end


%% Regularized reconstruction: alternate updating for magnitude mi and phase xi
if ~isvar('niter')
	niter = piter1 + piter2; % # of total iterations

	% 1st-order finite difference fatrix for regularization
	C = Cdiffs([Nd1(1) Nd1(2)], 'mask', mask);
	soft = @(t,a) (t - a * sign(t)) .* (abs(t) > a); % soft thresholding fun

  for iter = 1:niter
	ticker(mfilename, iter, niter)

	% update xi
	if iter > preiter
        	if iter > piter1 % minimize the cost function
			xi = pcg_bls_exp_ep(A1, C, yi, mi, xi, beta11, del, ...
				nsubiter(1), 'mask', mask);
			% with regularizer 4
%			xi = pcg_bls(A1, C, yi, mi, xi, beta1, nsubiter(1));
			% with regularizer 1
		else % optional 3rd step initialization
			if typ_init == 1 % with regularizer 2
		                xi = pcg_bls_exp(A1, C, yi, mi, xi, beta1, ...
					nsubiter(1), 'mask', mask);
			else % with regularizer 3
                		xi = pcg_bls_ep(A1, C, yi, mi, xi, beta11, ...
					del, nsubiter, 'mask', mask);
            		end
	        end
	end

	% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

	% update mi or img_i by soft thresholding
%	Ax0 = A1*diag_sp(exp(1i*xi)); % plug in the newest phase terms
	Ax0 = A1 * Gdiag(exp(1i*xi), 'mask', mask); % jf

	if iter <= preiter % 2nd step initialization by conventional CS
		tmp = img_i + 1/curv * A1' * (yi - A1 * img_i);
		img_i = U * soft(U' * tmp, beta22 / curv);
	else % CS update for magnitude image
		for subiter = 1:nsubiter(2)
			tmp = U' * mi + 1/curv * real(U'*Ax0'*(yi-Ax0*mi));
			mi = U * soft(tmp, beta22 / (2 * curv));
		end
	end

	if iter == preiter
		mi = abs(img_i);
		xi = angle(img_i);
	end

	if mod(iter,10) == 1
		printm('iteration %d/%d', iter, niter)
	end
 end % iter

%% Mask out the unreliable phase by mask
	xi = embed(xi,mask) .* msk1;
	mi = embed(mi,mask) .* msk1;

%% Phase unwrapping (not necessary, just for evaluation)
	mk = abs(xi-xt) > pi; % find the wrapped pixels
	xi_wp = xi(mk);
	xt_wp = xt(mk);
	del_wp = round((xt_wp-xi_wp)/(2*pi)); % # of periods needed to be compensated
	xi(mk) = xi(mk) + del_wp*2*pi; % unwrapping
end


%% Display the results
ir_fontsize im_axes 5
clim1 = [0 9]; % limits for color bar of magnitude results
clim2 = [-1 1] * 6; % limits for color bar of phase results
im plc 2 3
im(1, mt.', 'true',clim1), cbar
ylabelf 'magnitude'
im(4, xt.', ' ',clim2), cbar
ylabelf 'phase'

im(2, mi.', 'proposed', clim1), cbar
im(5, xi.', ' ', clim2),cbar

im(3, (mi-mt).', 'error', [-1 1]), cbar
im(6, (xi-xt).', ' ', [-1.5 1.5]), cbar
