% mag_phs_cs_real.m
% This script is an example of MRI reconstruction with separate magnitude and
% phase regularization via Compressed Sensing (CS) for real experiments.
% A human abdominal aorta velocity mapping data is used in this script. For
% details of the data, please refer to our paper.
%
% Copyright 2012-06-15, Feng Zhao, The University of Michigan
% 2013-04-11 modify for octave by JF

%% Load data
if ~isvar('yb')
	load abd_vl_real;
	img = imt26; % complex image reconstructed by DFT of fully sampled k-space data, frame 6
	yb = kdata2(:,:,6); % raw fully sampled k-space data, frame 6
	mask = msk2; % loose mask for image domain, frame 6
end


%% Setup reconstruction parameters
if ~isvar('mt')
	mt = abs(img); xt = angle(img);
	Nd = size(mt); % image size
	Np = prod(Nd); % # of pixels
	curv = prod(Nd); % spectral radius of A'*A

	sr = 0.3; % sampling rate without the fully sampled k-space center
	calirt = 0.04; % ratio of sampling for the fully sampled k-space center
	edgepr = 1; % use edge preserving smoothing regularizer or not

	piter1 = 80; % # of iterations step 3 of initialization
	piter2 = 80; % # of iterations for updating phase
	nsubiter = [2 1]; % number of subiterations in each iteration
	beta1 = 5e4*1; % regularization paramter for phase (rg1/rg3)
	typ_init = 0; % 0: initialize by rg3; 1: initialize by rg2
	beta11 = 6e4; % regularization paramter for phase (rg2/rg4)
	beta22 = 0.5*2^-4 * curv; % regularization paramter for magnitude
	if edgepr == 1;
		del = 0.005/1; % parameter for edge-preserving
	else
		piter2 = 0;
	end
	preiter = 2; % number of iterations for intialization
end


%% Setup system objects
if ~isvar('A')
	%rand('state', 1)
	rng(0)
	samp = rand(Nd) > 1 - sr;
	rc = calirt;
	i1 = round(Nd(1)*(1-sqrt(rc))/2);
	i2 = round(Nd(2)*(1-sqrt(rc))/2);
	samp(i1+1:Nd(1)-i1,i2+1:Nd(2)-i2) = true; % random sampling
	unsamp = sum(samp(:))/prod(Nd); % calculate the actual sampling rate;
	% A = Gdft('mask', mask, 'samp', samp, 'ifftshift',1,'fftshift',1,'class','Fatrix');
	A = Gdft('mask', mask, 'samp', samp, 'ifftshift', 1, 'fftshift', 1); % system fatrix


	%% Setup the wavlet transform object
	% U = Gwave1('mask', mask, 'noCA', 0, 'level', 3);
	U = Godwt1(mask, 'level', 3); % jf
	U = U';

	%% Simulate the MRI data acquistion
	yi = yb(samp);

	%% incorporate the reference phase map
	ref = im20;

	ref = ref(mask);
	%rfph = diag_sp(exp(1i*angle(ref)));
	rfph = Gdiag(exp(1i*angle(ref)), 'mask', mask); % jf
	A = A * rfph;
end


%% The 1st step initialization by inverse FFT
if ~isvar('xi')
	img_i = A' * yi / prod(Nd);
	mi = abs(img_i);
	xi = angle(img_i);
end

%% Regularized reconstruction: alternate updating for mt and xt
if ~isvar('iter')

niter = piter1+piter2; % number of total iterations
% C = Cdiffs([Nd(1),Nd(2)], 'mask', mask, 'order', 1, 'class','Fatrix'); % 1st order finite difference fatrix
C = Cdiffs([Nd(1),Nd(2)], 'mask', mask, 'order', 1); % jf
soft = @(t,a) (t - a * sign(t)) .* (abs(t) > a); % soft thresholding function

for iter = 1:niter
	% update xi
	if iter > preiter
		if iter > piter1 % minimize the cost function
			xi = pcg_bls_exp_ep(A, C, yi, mi, xi, beta11, del, ...
				nsubiter(1), 'mask', mask); % with regularizer 4
		else % optional 3rd step initialization
			if typ_init == 1 % with regularizer 2
				xi = pcg_bls_exp(A, C, yi, mi, xi, beta1, ...
					nsubiter(1), 'mask', mask);
			else % with regularizer 3
				xi = pcg_bls_ep(A, C, yi, mi, xi, beta11, ...
					del, nsubiter, 'mask', mask);
			end
		end
	end

	% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

	% update mi or img_i by soft thresholding
%	Ax0 = A*diag_sp(exp(1i*xi)); % plug in the newest phase terms
	Ax0 = A * Gdiag(exp(1i*xi), 'mask', mask); % jf

	if iter <= preiter
		% 2nd step initialization by conventional CS
		tmp = img_i + 1/curv * A' * (yi - A * img_i);
		img_i = U * soft(U' * tmp, beta22 / curv);
	else
		% CS update for magnitude image
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
		printf('iteration %d/%d', iter, niter)
	end
end

	xi = embed(xi,mask);
	mi = embed(mi,mask);

	mi = rot90(mi,2); mt = rot90(mt,2);
	im20 = rot90(im20,2); imt20 = rot90(imt20,2);
	xi = rot90(xi,2); xt = rot90(xt,2); dx26 = rot90(dx26,2);

end % iter

%% Display the results

clim1 = [0 max(mt(:))]; % limits for color bar of magnitude results
clim2 = [-pi pi]; % limits for color bar of phase results
vlim = [-3 3];
im plc 2 3
im(1, mt.', 'magnitude', clim1), cbar
ylabelf 'DFT'
im(2, angle(imt20).', 'ref. phase', clim2), cbar
im(3, (dx26/fc).', 'velocity', vlim), cbar

im(4, mi.', 'magnitude', clim1), cbar
ylabelf 'proposed'
im(5, angle(im20).', 'ref. phase', clim2), cbar
im(6, (xi/fc).', 'velocity', vlim),cbar
