function [kspace, y, dkf, phi, pmap, out] = ghost_correct(yi,ks,smap,varargin)
%| function [kspace, y, dkf, phi, pmap, out] = 
%|				ghost_correct(yi,ks,N,smap,varargin)
%|
%| Function to apply linear phase correction to EPI data for N/2 ghost
%| correction. Modifies kspace with linear shift (if model based
%| reconstruciton is used) and data by applying linear phase correction in
%| the image domain and the FFT (used for FFT reconstructions or to avoid
%| slower NUFFT)
%|
%| in
%|	yi		[nd ns nc]	Original dataset, nd=#samples per shot,
%|						ns=# of shots, and nc=# of coils
%|	ks		[nd*ns 2]	Original kspace sampling pattern
%|	smap	[nx ny nc]	Sensitivity maps
%|
%| options
%|	kshift	[2 1]		linear shift in kspace to apply correction in kx and ky
%|	phi		[3 1]		linear phase correction terms, used to create phase map
%|						as phi(1) + phi(2)*x + phi(3)*y
%|	phiinit [3 1]		Initial estimate of phase map parameters
%|	pinit	[nx ny]		Initial estimate of 2D phase map
%|	mask	[nx ny]		reconstruction mask
%|	emask	[nx ny]		tight object mask used for phase map estimation
%|	fmap	[nx ny]		fieldmap to use in image reconstruction
%|	ti		[nd 1]		time samples when fmap is used
%|	tseg	scalar		time segments for time segmentation when fmap used
%|	fov		[2 1]		Field of view in cm
%|	niter	scalar		# of CG iterations
%|	updates scalar		# of alternating updates for joint estimation
%|	chat	bool		Show printouts
%|	plot	bool		Plot figures
%|	basis	string		Basis function to use
%|		'dirac'			Dirac impulses
%|		'rect'			Rectangular basis
%|	etype	string		Phase map estimation type. All options from "fit" 
%|						function of "curve fitting" toolbox are applicable. 
%|						Recommended options are polynomial fit or local linear 
%|						regression. New option for joint estimation of phase map
%|						can be more effective in estimating linear 2D map. 
%|						Higher order estimation available up to 3rd oder 
%|						polynomial.
%|		'polyXY'		2D map from polynomial interpolation. Use X=Y=1 for 
%|						plane. Use poly10 for standard 1D correction
%|		'linearinterp'	2D map by piecewise linear interpolation
%|		'cubicinterp'	2D map by piecewise cubic interpolation
%|		'lowess'		2D map by local linear regression
%|		'regfmap'		2D map by using regularized fieldmap estimation
%|						method
%|		'joint'			Joint estimation of image and 2D planar phase map
%|		'vproj2d'		Estimation of 2D phase map using variable projection
%|		'vproj1d'		Estimation of 1D phase map using variable projection
%|		'approx2d'		Approximate FFT based 2D estimation
%|		'approx1d'		Approximate FFT based 1D estimation
%|	morder	scalar		Polynomial order for phase map parameter estimation
%|						Default is 1 for linear. Up to order 3 implemented
%|	indcoil bool		Turns on individual phase map calculation for each coil
%|						The output vector phi is [3 nc]
%|	indshot bool		Turns on individual phase map calculation for each shot
%|						The output vector phi is [3 ns]. If used in conjunction
%|						with the indcoil option then output is [3 ns nc].
%|	rows	[nrow 1]	Rows to use for initial phase map estimation
%|	cols	[ncol 1]	Columns to use for initial phase map estimation
%|	itnrms	scalar		Nrms difference between iterations used as stopping
%|						criterion for qs and joint estimation
%|
%| out
%|	kspace	[nd*ns 2]	Modified kspace with linear shift
%|	y		[nd ns nc]	Modified dataset after inage domain correction
%|						and FFT
%|	dkf		[2 1]		Linear shifts it kspace for subsequent correction
%|	phi		[3 1]		Linear terms for 2D phase map used for subsequent
%|						correction
%|	pmap	[nx ny ns nc]	2D phase map given as phi(1) + phi(2)*x + phi(3)*y
%|							possibly different map for each shot and coil
%|	out		struct		Miscelaneous output parameters
%|
%| Copyright 2010-11-18, Antonis Matakos, University of Michigan
%| Updated 2013-03-14, Antonis Matakos, University of Michigan

out = [];
[arg.nx, arg.ny, arg.nc] = size(smap);

arg.kshift = [];
arg.phi = [];
arg.phiinit = [];
arg.mask = true(arg.nx,arg.ny);
arg.emask = [];
arg.fmap = [];
arg.ti = [];
arg.tseg = [];
arg.fov = [24 24];
arg.niter = 20;
arg.riter = 100;
arg.updates = 20;
arg.chat = 0;
arg.plot = 0;
arg.basis = 'rect';
arg.etype = 'joint';
arg.morder = 1;
arg.indcoil = 0;
arg.indshot = 0;
arg.cols = 29:39;
arg.rows = 10:46;
arg.fft_osmp = 4;
arg.itnrms = 1e-5;

arg = vararg_pair(arg, varargin);

if length(arg.fov) == 1
	arg.fov(2) = arg.fov(1);
end


%% Split dataset
arg.N = [arg.nx arg.ny];
arg.nd = size(yi,1);
arg.ns = size(yi,2);

if arg.indcoil && arg.indshot
	[Gp, Gn, yp, yn, ksp, ksn, arg] = create_system_shot_coil(yi, ks, smap, arg);
elseif arg.indshot
	[Gp, Gn, yp, yn, ksp, ksn, arg] = create_system_shot(yi, ks, smap, arg);
elseif arg.indcoil
	[Gp, Gn, yp, yn, ksp, ksn, arg] = create_system_coil(yi, ks, smap, arg);
else
	[Gp, Gn, yp, yn, ksp, ksn, arg] = create_system(yi, ks, smap, arg);
end


arg.xx = (1:arg.nx)'-arg.nx/2-1;
arg.yy = (1:arg.ny)'-arg.ny/2-1;
[arg.X,arg.Y] = ndgrid(arg.xx,arg.yy);

if arg.morder == 1
	T = [ones(sum(arg.mask(:)),1) arg.X(arg.mask) arg.Y(arg.mask)];
	Tr = [ones(numel(arg.mask),1) arg.X(:) arg.Y(:)];
elseif arg.morder == 2
	T = [ones(sum(arg.mask(:)),1) arg.X(arg.mask) arg.Y(arg.mask) ...
		arg.X(arg.mask).^2 arg.X(arg.mask).*arg.Y(arg.mask) arg.Y(arg.mask).^2];
	Tr = [ones(numel(arg.mask),1) arg.X(:) arg.Y(:) arg.X(:).^2 ...
		arg.X(:).*arg.Y(:) arg.Y(:).^2];
elseif arg.morder == 3
	T = [ones(sum(arg.mask(:)),1) arg.X(arg.mask) arg.Y(arg.mask) ...
		arg.X(arg.mask).^2 arg.X(arg.mask).*arg.Y(arg.mask) arg.Y(arg.mask).^2 ...
		arg.X(arg.mask).^3 arg.X(arg.mask).^2.*arg.Y(arg.mask) ...
		arg.X(arg.mask).*arg.Y(arg.mask).^2 arg.Y(arg.mask).^3];
	Tr = [ones(numel(arg.mask),1) arg.X(:) arg.Y(:) arg.X(:).^2 ...
		arg.X(:).*arg.Y(:) arg.Y(:).^2 arg.X(:).^3 arg.X(:).^2.*arg.Y(:) ...
		arg.X(:).*arg.Y(:).^2 arg.Y(:).^3];
else
	warn('Higer order than 3 not implemeted. Reverting to default order=1.')
	arg.morder = 1;
	T = [ones(sum(arg.mask(:)),1) arg.X(arg.mask) arg.Y(arg.mask)];
	Tr = [ones(numel(arg.mask),1) arg.X(:) arg.Y(:)];
end

if isempty(arg.emask)
	arg.emask = false(arg.nx, arg.ny);
	arg.emask(arg.rows, arg.cols) = true;
end

if ~isempty(arg.phi)
	
	phi = arg.phi;
	if size(phi,3) > 1
		phi = reshape(phi,size(phi,1),[]);
% 		pmap = embed(T*phi, arg.mask);
% 		pmap = reshape(pmap,arg.nx,arg.ny,[],arg.nc);
		phi = reshape(phi,size(phi,1),[],arg.nc);
	end
	dkf = phi(2:3,:)*arg.ny/2/pi/arg.fov(2);
	
elseif ~isempty(arg.kshift)
	
	dkf = arg.kshift(:);
	phi = [0; 2*pi*arg.fov(2)/arg.ny*dkf];
	
% 	pmap = embed(T*phi, arg.mask);
	
else
	if strcmpi(arg.etype, 'regfmap')
		[pinit, phi, out] = fit_est('poly11', Gp, Gn, yp, yn, arg);
		xg(:,:,1) = out.xtn;
		xg(:,:,2) = out.xtp;
		pmap = mri_field_map_reg(xg, [0 1], 'winit', pinit, 'mask', arg.mask);
	elseif strcmpi(arg.etype, 'joint')
		if isempty(arg.phiinit)
			[Gpi, Gni, ypi, yni] = create_system(yi, ks, smap, arg);
			[~, phi] = fit_est('poly11', Gpi, Gni, ypi, yni, arg);
		else
			phi = arg.phiinit;
		end
		
		if arg.indcoil && size(phi,3) == 1
			phi = repmat(phi, [1 1 arg.nc]);
		end
		
		if arg.indshot && size(phi,2) == 1
			phi = repmat(phi, [1 2*arg.ns 1]);
			phi(:,1:arg.ns,:) = 0;
		end
		
		phi = squeeze(phi);
		
		if arg.indcoil && arg.indshot
			[~, phi] = joint_est_shot_coil(T, phi, Gp, Gn, yp, yn, arg);
		elseif arg.indshot
			[~, phi] = joint_est_shot(T, phi, Gp, Gn, yp, yn, arg);
		elseif arg.indcoil
			[~, phi] = joint_est_coil(T, phi, Gp, Gn, yp, yn, arg);
		else
			[~, phi] = joint_est(T, phi, Gp, Gn, yp, yn, arg);
		end
	elseif strcmpi(arg.etype, 'vproj2d')
		if isempty(arg.phiinit)
			[~, phi] = approx2d_est(T, smap, Gp, Gn, yp, yn, arg);
		else
			phi = arg.phiinit;
		end
		
		[~, phi] = vproj2d_est(T, phi, smap, Gp, Gn, yp, yn, arg);
		
	elseif strcmpi(arg.etype, 'vproj1d')
		if isempty(arg.phiinit)
			[~, phi] = approx1d_est(T, smap, Gp, Gn, yp, yn, arg);
		else
			phi = arg.phiinit;
		end
		
		if arg.indshot && size(phi,2) == 1
			phi = repmat(phi, [1 2*arg.ns 1]);
			phi(:,1:arg.ns,:) = 0;
		end
		
		phi = squeeze(phi);
		
		if arg.indshot
			[~, phi] = vproj1d_est_shot(T, phi, smap, Gp, Gn, yp, yn, arg);
		else
			[~, phi] = vproj1d_est(T, phi, smap, Gp, Gn, yp, yn, arg);
		end
		
	elseif strcmpi(arg.etype, 'approx2d')
		[~, phi] = approx2d_est(T, smap, Gp, Gn, yp, yn, arg);
	elseif strcmpi(arg.etype, 'approx1d')
		[~, phi] = approx1d_est(T, smap, Gp, Gn, yp, yn, arg);
	else
		if arg.indcoil
			[~, phi] = fit_est_coil(arg.etype, Gp, Gn, yp, yn, arg);
		elseif arg.indshot
			[~, phi] = fit_est_shot(arg.etype, Gp, Gn, yp, yn, arg);
		else
			[~, phi] = fit_est(arg.etype, Gp, Gn, yp, yn, arg);
		end
	end
	
	while phi(1) < -pi
		phi(1) = phi(1) + 2*pi;
	end
	
	while phi(1) >= pi
		phi(1) = phi(1) - 2*pi;
	end
	
	% Calculate offset using linear fit
	dkf = phi(2:3,:,:)*arg.ny/2/pi/arg.fov(2);
end


%% Apply correction to kspace
[kspace, yk] = kspace_correction(ksp, ksn, yp, yn, dkf, phi, arg);

%% Apply correction to data
y = data_correction(yp, yn, phi, arg);


%% Test correction by reconstructing image
arg.mask = true(arg.N);
if ~strcmpi(arg.etype, 'regfmap')
	pmap = embed(Tr*reshape(phi,size(phi,1),[]), arg.mask);
end

if arg.plot
	if strcmpi(arg.etype, 'joint');
		
		if arg.indshot && arg.indcoil
			pmap = reshape(pmap, arg.nx, arg.ny, [], arg.nc);
		end
	end
	if arg.indshot
		out = plot_results_shot(yi, yp, yn, y, yk, ks, kspace, ...
			smap, pmap, arg);
	else
		out = plot_results(yi, yp, yn, y, yk, ks, kspace, ...
			ksp, ksn, smap, pmap, arg);
	end
end

out.yk = yk;

end

function [Gp, Gn, yp, yn, ksp, ksn, arg] = create_system(yi, ks, smap, arg)

mask = arg.mask;
y = reshape(yi,arg.nx,[],arg.ns,arg.nc);
yp = y(:,1:2:end,:,:); % data going right
yn = y(:,2:2:end,:,:); % data goint left

ksx = reshape(ks(:,1),arg.nx,[]);
ksxp = ksx(:,1:2:end); % data going right
ksxn = ksx(:,2:2:end); % data goint left

ksy = reshape(ks(:,2),arg.nx,[]);
ksyp = ksy(:,1:2:end); % data going right
ksyn = ksy(:,2:2:end); % data goint left

ksp = [ksxp(:) ksyp(:)];
ksn = [ksxn(:) ksyn(:)];

if ~isempty(arg.ti)
	arg.ti = reshape(arg.ti,arg.nx,[]);
	arg.tip = vec(arg.ti(:,1:2:end));
	arg.tin = vec(arg.ti(:,2:2:end));
end

nufft = {arg.N, [1 1], arg.N, arg.N/2, 'linear'};

Gmp = Gmri(ksp, arg.mask, 'fov', arg.fov, 'basis', {arg.basis},'nufft', nufft);
Gmn = Gmri(ksn, arg.mask, 'fov', arg.fov, 'basis', {arg.basis},'nufft', nufft);

if ~isempty(arg.fmap)
	Gmp = feval(Gmp.arg.new_zmap, Gmp, arg.tip(:), arg.fmap, arg.tseg);
	Gmn = feval(Gmn.arg.new_zmap, Gmn, arg.tin(:), arg.fmap, arg.tseg);
end

Gcp = cell(1,arg.nc);
Gcn = cell(1,arg.nc);

parfor ic = 1:arg.nc
	tmp = smap(:,:,ic);
	Gcp{ic} = Gmp * Gdiag(tmp(mask),'mask',mask); % cascade
	Gcn{ic} = Gmn * Gdiag(tmp(mask),'mask',mask); % cascade
end
Gp = block_fatrix(Gcp, 'type', 'col'); % [G1; G2; ... ]
Gn = block_fatrix(Gcn, 'type', 'col'); % [G1; G2; ... ]

end

function [Gp, Gn, yp, yn, ksp, ksn, arg] = create_system_coil(yi, ks, smap, arg)

mask = arg.mask;
y = reshape(yi,arg.nx,[],arg.ns,arg.nc);
yp = y(:,1:2:end,:,:); % data going right
yn = y(:,2:2:end,:,:); % data goint left

ksx = reshape(ks(:,1),arg.nx,[]);
ksxp = ksx(:,1:2:end); % data going right
ksxn = ksx(:,2:2:end); % data goint left

ksy = reshape(ks(:,2),arg.nx,[]);
ksyp = ksy(:,1:2:end); % data going right
ksyn = ksy(:,2:2:end); % data goint left

ksp = [ksxp(:) ksyp(:)];
ksn = [ksxn(:) ksyn(:)];

if ~isempty(arg.ti)
	arg.ti = reshape(arg.ti,arg.nx,[]);
	arg.tip = vec(arg.ti(:,1:2:end));
	arg.tin = vec(arg.ti(:,2:2:end));
end

nufft = {arg.N, [1 1], arg.N, arg.N/2, 'linear'};

Gmp = Gmri(ksp, arg.mask, 'fov', arg.fov, 'basis', {arg.basis},'nufft', nufft);
Gmn = Gmri(ksn, arg.mask, 'fov', arg.fov, 'basis', {arg.basis},'nufft', nufft);

if ~isempty(arg.fmap)
	Gmp = feval(Gmp.arg.new_zmap, Gmp, arg.tip(:), arg.fmap, arg.tseg);
	Gmn = feval(Gmn.arg.new_zmap, Gmn, arg.tin(:), arg.fmap, arg.tseg);
end

Gp = cell(1,arg.nc);
Gn = cell(1,arg.nc);

parfor ic = 1:arg.nc
	tmp = smap(:,:,ic);
	Gp{ic} = Gmp * Gdiag(tmp(mask),'mask',mask); % cascade
	Gn{ic} = Gmn * Gdiag(tmp(mask),'mask',mask); % cascade
end

end

function [Gp, Gn, yp, yn, ksp, ksn, arg] = create_system_shot(yi, ks, smap, arg)

mask = arg.mask;
fov = arg.fov;
basis = arg.basis;
fmap = arg.fmap;

tseg = arg.tseg;
nc = arg.nc;

y = reshape(yi,arg.nx,[],arg.ns,arg.nc);
yp = y(:,1:2:end,:,:); % data going right
yn = y(:,2:2:end,:,:); % data goint left

ksx = reshape(ks(:,1),arg.nx,[],arg.ns);
ksxp = ksx(:,1:2:end,:); % data going right
ksxn = ksx(:,2:2:end,:); % data goint left

ksy = reshape(ks(:,2),arg.nx,[],arg.ns);
ksyp = ksy(:,1:2:end,:); % data going right
ksyn = ksy(:,2:2:end,:); % data goint left

if ~isempty(arg.ti)
	arg.ti = reshape(arg.ti,arg.nx,[],arg.ns);
	arg.tip = arg.ti(:,1:2:end,:);
	arg.tin = arg.ti(:,2:2:end,:);
	tip = arg.tip;
	tin = arg.tin;
end

nufft = {arg.N, [1 1], arg.N, arg.N/2, 'linear'};

Gp = cell(1,arg.ns);
Gn = cell(1,arg.ns);

for ii=1:arg.ns
	
	ksp = [vec(ksxp(:,:,ii)) vec(ksyp(:,:,ii))];
	ksn = [vec(ksxn(:,:,ii)) vec(ksyn(:,:,ii))];
	
	Gmp = Gmri(ksp, mask, 'fov', fov, 'basis', {basis}, 'nufft', nufft);
	Gmn = Gmri(ksn, mask, 'fov', fov, 'basis', {basis}, 'nufft', nufft);
	
	if ~isempty(fmap)
		Gmp = feval(Gmp.arg.new_zmap, Gmp, vec(tip(:,:,ii)), fmap, tseg);
		Gmn = feval(Gmn.arg.new_zmap, Gmn, vec(tin(:,:,ii)), fmap, tseg);
	end
	
	Gcp = cell(1,nc);
	Gcn = cell(1,nc);
	
	for ic = 1:nc
		tmp = smap(:,:,ic);
		Gcp{ic} = Gmp * Gdiag(tmp(mask),'mask',mask); % cascade
		Gcn{ic} = Gmn * Gdiag(tmp(mask),'mask',mask); % cascade
	end
	Gp{ii} = block_fatrix(Gcp, 'type', 'col'); % [G1; G2; ... ]
	Gn{ii} = block_fatrix(Gcn, 'type', 'col'); % [G1; G2; ... ]
end

ksp = [ksxp(:) ksyp(:)];
ksn = [ksxn(:) ksyn(:)];

end

function [Gp, Gn, yp, yn, ksp, ksn, arg] = create_system_shot_coil(yi, ks, smap, arg)

mask = arg.mask;
fov = arg.fov;
basis = arg.basis;
fmap = arg.fmap;
tseg = arg.tseg;
nc = arg.nc;

y = reshape(yi,arg.nx,[],arg.ns,arg.nc);
yp = y(:,1:2:end,:,:); % data going right
yn = y(:,2:2:end,:,:); % data goint left

ksx = reshape(ks(:,1),arg.nx,[],arg.ns);
ksxp = ksx(:,1:2:end,:); % data going right
ksxn = ksx(:,2:2:end,:); % data goint left

ksy = reshape(ks(:,2),arg.nx,[],arg.ns);
ksyp = ksy(:,1:2:end,:); % data going right
ksyn = ksy(:,2:2:end,:); % data goint left

if ~isempty(arg.ti)
	arg.ti = reshape(arg.ti,arg.nx,[],arg.ns);
	arg.tip = arg.ti(:,1:2:end,:);
	arg.tin = arg.ti(:,2:2:end,:);
	tip = arg.tip;
	tin = arg.tin;
end

nufft = {arg.N, [1 1], arg.N, arg.N/2, 'linear'};

Gp = cell(arg.ns, arg.nc);
Gn = cell(arg.ns, arg.nc);

for ii=1:arg.ns
	
	ksp = [vec(ksxp(:,:,ii)) vec(ksyp(:,:,ii))];
	ksn = [vec(ksxn(:,:,ii)) vec(ksyn(:,:,ii))];
	
	Gmp = Gmri(ksp, mask, 'fov', fov, 'basis', {basis},'nufft', nufft);
	Gmn = Gmri(ksn, mask, 'fov', fov, 'basis', {basis},'nufft', nufft);
	
	if ~isempty(fmap)
		Gmp = feval(Gmp.arg.new_zmap, Gmp, vec(tip(:,:,ii)), fmap, tseg);
		Gmn = feval(Gmn.arg.new_zmap, Gmn, vec(tin(:,:,ii)), fmap, tseg);
	end
	
	parfor ic = 1:nc
		tmp = smap(:,:,ic);
		Gp{ii,ic} = Gmp * Gdiag(tmp(mask),'mask',mask); % cascade
		Gn{ii,ic} = Gmn * Gdiag(tmp(mask),'mask',mask); % cascade
	end
end

ksp = [ksxp(:) ksyp(:)];
ksn = [ksxn(:) ksyn(:)];

end

function [pmap, phi] = joint_est(T, phi, Gp, Gn, yp, yn, arg)

xj = zeros(arg.N);
pmap = embed(T*phi, arg.mask);
phiold = phi;

% figure(1), im clf
for ii=1:arg.updates
	printm('Update %d',ii);
	Dm = Gdiag(exp(-1i*pmap(arg.mask)), 'mask', arg.mask);
	Gj = block_fatrix({Gp, Gn*Dm}, 'type', 'col');
	
	ticker reset, ticker print, ticker(4);
	xi = qpwls_pcg1(xj(arg.mask), Gj, 1, [yp(:);yn(:)], 0,'niter', arg.niter);
	xj = embed(xi, arg.mask);
	Dx = Gdiag(xj(arg.mask), 'mask', arg.mask);
	Gx = Gn * (Dm * (Dx*T));
	
	ynt = yn(:) - Gn*(Dm*xj(arg.mask)) - 1i*Gn*(Dm*(Dx*pmap(arg.mask)));
	phi = real(1i*(Gx\ynt));
	pmap = embed(T*phi, arg.mask);
% 	figure(1),im(pmap),drawnow;
% 	figure(2),im(xj),drawnow;
	diff = nrms(phi,phiold);
	if diff < arg.itnrms
		printm('Stopped at iteration %d',ii);
		break;
	end
	phiold = phi;
end


end

function [pmap, phi] = joint_est_coil(T, phi, Gp, Gn, yp, yn, arg)

mask = arg.mask;
ym = [yp(:);yn(:)];

xj = zeros(arg.N);
pmap = embed(T*phi, arg.mask);
phiold = phi;

Gcn = cell(1,arg.nc);

Gbp = block_fatrix(Gp, 'type', 'col');

% figure(1), im clf
for ii=1:arg.updates
	printm('Update %d',ii);
	parfor ic=1:arg.nc
		tmp = exp(-1i*pmap(:,:,ic));
		Gcn{ic} = Gn{ic} * Gdiag(tmp(mask), 'mask', mask);
	end
	Gbn = block_fatrix(Gcn, 'type', 'col');
	Gj = block_fatrix({Gbp, Gbn}, 'type', 'col');
	
	ticker reset, ticker print, ticker(4);
	xi = qpwls_pcg1(xj(arg.mask), Gj, 1, ym, 0,'niter', arg.niter);
	xj = embed(xi, arg.mask);
	Dx = Gdiag(xj(arg.mask), 'mask', arg.mask);
	
	parfor ic=1:arg.nc
		Gx = Gcn{ic} * (Dx * T);
		tmp = pmap(:,:,ic);
		ynt = vec(yn(:,:,:,ic))-Gcn{ic}*xi-1i*Gcn{ic}*(Dx*tmp(mask));
		phi(:,ic) = real(1i*(Gx\ynt));
	end
	pmap = embed(T*phi, arg.mask);
	diff = nrms(phi,phiold);
	if diff < arg.itnrms
		printm('Stopped at iteration %d',ii);
		break;
	end
	phiold = phi;
% 	figure(1),im(pmap),drawnow;
% 	figure(2),im(xj),drawnow;
end

end

function [pmap, phi] = joint_est_shot(T, phi, Gp, Gn, yp, yn, arg)

mask = arg.mask;
ym = [vec(permute(yp,[1 2 4 3]));vec(permute(yn,[1 2 4 3]))];
ns = arg.ns;

xj = zeros(arg.N);
pmap = embed(T*phi, arg.mask);
phiold = phi;
xold = xj;

Gcn = cell(1,arg.ns);
Gcp = cell(1,arg.ns);

% figure(1), im clf
for ii=1:arg.updates
	printm('Update %d',ii);
	Gcp{1} = Gp{1};
	parfor ic=1:arg.ns
		tmp = exp(-1i*pmap(:,:,ic));
		Gcp{ic} = Gp{ic} * Gdiag(tmp(mask), 'mask', mask);
	end
	Gbp = block_fatrix(Gcp, 'type', 'col');
	
	parfor ic=1:arg.ns
		tmp = exp(-1i*pmap(:,:,ns+ic));
		Gcn{ic} = Gn{ic} * Gdiag(tmp(mask), 'mask', mask);
	end
	Gbn = block_fatrix(Gcn, 'type', 'col');
	
	Gj = block_fatrix({Gbp, Gbn}, 'type', 'col');
	
	ticker reset, ticker print, ticker(4);
	xi = qpwls_pcg1(xj(arg.mask), Gj, 1, ym, 0,'niter', arg.niter);
	xj = embed(xi, arg.mask);
	Dx = Gdiag(xj(arg.mask), 'mask', arg.mask);
	
	parfor ic=2:arg.ns
		Gx = Gcp{ic} * (Dx * T);
		tmp = pmap(:,:,ic);
		ypt = vec(yp(:,:,ic,:))-Gcp{ic}*xi-1i*Gcp{ic}*(Dx*tmp(mask));
		phi(:,ic) = real(1i*(Gx\ypt));
	end
	
	parfor ic=1:arg.ns
		Gx = Gcn{ic} * (Dx * T);
		tmp = pmap(:,:,ns + ic);
		ynt = vec(yn(:,:,ic,:))-Gcn{ic}*xi-1i*Gcn{ic}*(Dx*tmp(mask));
		phi(:,ns+ic) = real(1i*(Gx\ynt));
	end
	
	pmap = embed(T*phi, arg.mask);
% 	diff = nrms(phi,phiold);
	diff = nrms(xold,xj);
	printm('Iteration nrms: %2.4f%%', diff*100);
	if diff < arg.itnrms
		printm('Stopped at iteration %d',ii);
		break;
	end
	phiold = phi;
	xold = xj;
% 	figure(1),im(pmap),drawnow;
% 	figure(2),im(xj),drawnow;
end

end

function [pmap, phi] = joint_est_shot_coil(T, phi, Gp, Gn, yp, yn, arg)

ns = arg.ns;
mask = arg.mask;
ym = [yp(:);yn(:)];
% ym = [vec(permute(yp,[1 2 4 3]));vec(permute(yn,[1 2 4 3]))];

phip = zeros(size(phi,1), ns,arg.nc);
phip(:,2:end,:) = phi(:,1:ns-1,:);
phip = reshape(phip, size(phi,1), []);

phin = reshape(phi(:,ns+1:end,:), size(phi,1), []);

xj = zeros(arg.N);
pmapp = embed(T*phip, arg.mask);
pmapn = embed(T*phin, arg.mask);

Gcn = cell(1,arg.ns*arg.nc);
Gcp = cell(1,arg.ns*arg.nc);

figure(1), im clf
for ii=1:arg.updates
	printm('Update %d',ii);
	parfor ic=1:arg.ns*arg.nc
		tmp = exp(-1i*pmapp(:,:,ic));
		Gcp{ic} = Gp{ic} * Gdiag(tmp(mask), 'mask', mask);
	end
	Gbp = block_fatrix(Gcp, 'type', 'col');
	
	parfor ic=1:arg.ns*arg.nc
		tmp = exp(-1i*pmapn(:,:,ic));
		Gcn{ic} = Gn{ic} * Gdiag(tmp(mask), 'mask', mask);
	end
	Gbn = block_fatrix(Gcn, 'type', 'col');
	
	Gj = block_fatrix({Gbp, Gbn}, 'type', 'col');
	
	ticker reset, ticker print, ticker(4);
	xi = qpwls_pcg1(xj(arg.mask), Gj, 1, ym, 0,'niter', arg.niter);
	xj = embed(xi, arg.mask);
	Dx = Gdiag(xj(arg.mask), 'mask', arg.mask);
	
	parfor ic=1:arg.ns*arg.nc
		Gx = Gcp{ic} * (Dx * T);
		tmp = pmapp(:,:,ic);
		ypt = vec(yp(:,:,ic))-Gcp{ic}*xi-1i*Gcp{ic}*(Dx*tmp(mask));
		phip(:,ic) = real(1i*(Gx\ypt));
	end
	phip = reshape(phip, size(phi,1), arg.ns, arg.nc);
% 	phip(:,1,:) = 0;
	phip = reshape(phip, size(phi,1), []);
	
	parfor ic=1:arg.ns*arg.nc
		Gx = Gcn{ic} * (Dx * T);
		tmp = pmapn(:,:,ic);
		ynt = vec(yn(:,:,ic))-Gcn{ic}*xi-1i*Gcn{ic}*(Dx*tmp(mask));
		phin(:,ic) = real(1i*(Gx\ynt));
	end
	pmapp = embed(T*phip, arg.mask);
	pmapn = embed(T*phin, arg.mask);
	figure(1),im(pmapp),drawnow;
	figure(2),im(pmapn),drawnow;
	figure(3),im(xj),drawnow;
end

pmap = zeros(arg.nx, arg.ny, 2*arg.ns, arg.nc);
pmap(:,:,1:arg.ns,:) = reshape(pmapp, arg.nx, arg.ny, arg.ns, arg.nc);
pmap(:,:,arg.ns+1:end,:) = reshape(pmapn, arg.nx, arg.ny, arg.ns, arg.nc);

phi = zeros(size(phi,1), 2*arg.ns, arg.nc);
phi(:,1:arg.ns,:) = reshape(phip, size(phi,1), arg.ns, arg.nc);
phi(:,arg.ns+1:end,:) = reshape(phin, size(phi,1), arg.ns, arg.nc);

end

function [pmap, phi] = vproj2d_est(T, phi, smap, Gp, Gn, yp, yn, arg)

tm = tic;
vp = Gp'*yp(:);
vn = Gn'*yn(:);
sw = sum(abs(smap).^2,3);
sw = sw(arg.mask);

wm = (abs(vp).*abs(vn))./sw;
delta = angle(vn) - angle(vp);

phiold = [0;0;0];

ts = toc(tm);
printm('Setup time: %5.2fms', ts*1000);
for ii=1:arg.updates
	sm = T*phi + delta;
	grad = T'*(wm.*sin(sm));
	srm = mod(sm + pi, 2*pi) - pi;
	denom = wm.*nufft_sinc_nopi(srm);
	H = T' * (diag_sp(denom) * T);
	phi = phi - (H \ grad);
	pmap = embed(T*phi, arg.mask);
% 	figure(1),im(pmap),drawnow;
	diff = nrms(phi,phiold);
	if diff < arg.itnrms
		printm('Stopped at iteration %d',ii);
		break;
	end
	phiold = phi;
end
time = toc(tm);
printm('Total time: %5.2fms', time*1000);

end

function [pmap, phi] = vproj1d_est(T, phi, smap, Gp, Gn, yp, yn, arg)

phi(3,:) = [];
T = T(:,1:2);
vp = Gp'*yp(:);
vn = Gn'*yn(:);
sw = sum(abs(smap).^2,3);
sw = sw(arg.mask);

wm = (abs(vp).*abs(vn))./sw;
delta = angle(vn) - angle(vp);

for ii=1:arg.updates
	sm = T*phi + delta;
	grad = T'*(wm.*sin(sm));
	srm = mod(sm + pi, 2*pi) - pi;
	denom = wm.*nufft_sinc_nopi(srm);
	H = T' * (diag_sp(denom) * T);
	phi = phi - (H \ grad);
	pmap = embed(T*phi, arg.mask);
	figure(1),im(pmap),drawnow;
end

phi = [phi;0];

end

function [pmap, phi] = vproj1d_est_shot(T, phi, Gp, Gn, yp, yn, arg)

phi(3,:) = [];
phi = phi(:);
T = T(:,1:2);
q = 2*arg.ns;
v = zeros(sum(arg.mask),q);

for ii=1:arg.ns
	v(:,ii) = Gp{ii}'*vec(yp(:,:,ii,:))
	v(:,ii+arg.ns) = Gn{ii}'*vec(yn(:,:,ii,:))
end

sw = sum(abs(smap).^2,3);
sw = sw(arg.mask);

wm = zeros(sum(arg.mask),q*(q-1)/2);
delta = zeros(sum(arg.mask),q*(q-1)/2);
Ts = cell(q*(q-1)/2,q);
idx = 1;
for ii=1:q-1
	for jj=(ii+1):q
		wm(:,idx) = (abs(v(:,ii)).*abs(v(:,jj)))./sw;
		delta(:,idx) = angle(v(:,ii)) - angle(v(:,jj));
		Ts{idx,ii} = T;
		Ts{idx,jj} = -T;
		idx = idx+1;
	end
end

wm = wm(:);
delta = delta(:);
Ts = cell2mat(Ts);

for ii=1:arg.updates
	sm = Ts*phi + delta;
	grad = Ts'*(wm.*sin(sm));
	srm = mod(sm + pi, 2*pi) - pi;
	denom = wm.*nufft_sinc_nopi(srm);
	H = Ts' * (diag_sp(denom) * Ts);
	phi = phi - (H \ grad);
end

phi = reshape(phi,[],q);
pmap = embed(T*phi, arg.mask);
figure(1),im(pmap),drawnow;
phi(3,:) = 0;

end

function [pmap, phi] = approx2d_est(T, smap, Gp, Gn, yp, yn, arg)

tm = tic;
vp = Gp'*yp(:);
vn = Gn'*yn(:);
sw = sum(abs(smap).^2,3);
sw = sw(arg.mask);
% sw = ones(sum(arg.mask(:)),1);

ts = toc(tm);
printm('Setup time: %5.2fms', ts*1000);

mp = arg.fft_osmp;
M = ((mp-1)/2)*arg.N;

u = embed((conj(vp).*vn)./sw,arg.mask);

u = padarray(u,M,'both');
U = fftshift(ifft2(ifftshift(u)));

[~,idx] = max(abs(U(:)));
[nx,ny] = ind2sub(mp*arg.N,idx);

px = 2*pi*(nx - mp*arg.N(1)/2-1)/mp/arg.N(1);
py = 2*pi*(ny - mp*arg.N(2)/2-1)/mp/arg.N(2);
pc = mod(-angle(U(nx,ny)),2*pi);

phi = [pc;px;py];
pmap = embed(T*phi, arg.mask);

time = toc(tm);
printm('Total time: %5.2fms', time*1000);

end

function [pmap, phi] = approx1d_est(T, smap, Gp, Gn, yp, yn, arg)

tm = tic;
vp = Gp'*yp(:);
vn = Gn'*yn(:);
sw = sum(abs(smap).^2,3);
sw = sw(arg.mask);

ts = toc(tm);
printm('Setup time: %5.2fms', ts*1000);

mp = arg.fft_osmp;
M = ((mp-1)/2)*arg.N;

u = embed((conj(vp).*vn)./sw,arg.mask);
u = sum(u,2);
u = padarray(u,[M(1) 0],'both');
U = fftshift(ifft(ifftshift(u)));
[~,nx] = max(abs(U(:)));
px = 2*pi*(nx - mp*arg.N(1)/2-1)/mp/arg.N(1);
pc = mod(-angle(U(nx)),2*pi);

phi = [pc;px;0];
pmap = embed(T*phi, arg.mask);

time = toc(tm);
printm('Total time: %5.2fms', time*1000);


end

function [pmap, phi] = fit_est_coil(type, Gp, Gn, yp, yn, arg)

pmap = zeros(arg.nx, arg.ny, arg.nc);
phi = zeros(3,arg.nc);

parfor ic=1:arg.nc
	[pm, ph] = fit_est(type, Gp{ic}, Gn{ic}, yp(:,:,:,ic), yn(:,:,:,ic), arg);
	pmap(:,:,ic) = pm;
	phi(:,ic) = ph;
end

end

function [pmap, phi] = fit_est_shot(type, Gp, Gn, yp, yn, arg)

pmap = zeros(arg.nx, arg.ny, 2*arg.ns);
phi = zeros(3,2*arg.ns);

Gr = Gp{1};
yr = yp(:,:,1,:);

pmapn = zeros(arg.nx, arg.ny, arg.ns);
phin = zeros(3,arg.ns);

parfor ic=1:arg.ns
	[pm, ph] = fit_est(type, Gr, Gn{ic}, yr, yn(:,:,ic,:), arg);
	pmapn(:,:,ic) = pm;
	phin(:,ic) = ph;
end

pmap(:,:,arg.ns+1:end) = pmapn;
phi(:,arg.ns+1:end) = phin;

parfor ic=2:arg.ns
	[pm, ph] = fit_est(type, Gr, Gp{ic}, yr, yp(:,:,ic,:), arg);
	pmap(:,:,ic) = pm;
	phi(:,ic) = ph;
end

end

function [pmap, phi, out] = fit_est(type, Gp, Gn, yp, yn, arg)

xinit = zeros(arg.N);
tm = tic;
ticker reset, ticker print, ticker(4);
% xi = qpwls_pcg1(xinit(arg.mask), Gp, 1, yp(:), 0,'niter', arg.niter);
xi = Gp'*yp(:);
xtp = embed(xi, arg.mask);

ticker reset, ticker print, ticker(4);
% xi = qpwls_pcg1(xinit(arg.mask), Gn, 1, yn(:), 0,'niter', arg.niter);
xi = Gn'*yn(:);
xtn = embed(xi, arg.mask);
ts = toc(tm);
printm('Setup time: %5.2fms', ts*1000);
out.xtp = xtp;
out.xtn = xtn;

% pmap = angle(xtp(arg.mask))-angle(xtn(arg.mask));
% pmap = embed(pmap, arg.mask);
% phi = [0;0;0];


pmap = mod(angle(xtp(arg.emask))-angle(xtn(arg.emask)),2*pi);
pmap = embed(pmap, arg.emask);
pmap = unwrap(pmap);

if strcmpi(type, 'poly10')
	T = [ones(sum(arg.emask(:)),1) arg.X(arg.emask)];
	Tr = [ones(numel(arg.mask),1) arg.X(:)];
	phi = T\pmap(arg.emask);
	pmap = reshape(Tr*phi, arg.nx, arg.ny);
	phi = [phi;0];
elseif strcmpi(type, 'poly11')
	T = [ones(sum(arg.emask(:)),1) arg.X(arg.emask) arg.Y(arg.emask)];
	Tr = [ones(numel(arg.mask),1) arg.X(:) arg.Y(:)];
	phi = T\pmap(arg.emask);
	pmap = reshape(Tr*phi, arg.nx, arg.ny);
else
	sf = fit([arg.X(arg.emask) arg.Y(arg.emask)], pmap(arg.emask), type);
	pmap = sf(arg.X, arg.Y);
	sf = fit([arg.X(:) arg.Y(:)], pmap(:), 'poly11');
	phi = [sf.p00;sf.p10;sf.p01];
end
ts = toc(tm);
printm('Total time: %5.2fms', ts*1000);


while phi(1,1) < -pi
	phi(1,:) = phi(1,1) + 2*pi;
end

while phi(1) >= pi
	phi(1,:) = phi(1,1) - 2*pi;
end

end

function [ks, yk] = kspace_correction(ksp, ksn, yp, yn, dkf, phi, arg)

yk = zeros(arg.nx, arg.ny/arg.ns, arg.ns, arg.nc);

ksx = zeros(arg.nx, arg.ny/arg.ns, arg.ns);
ksy = zeros(arg.nx, arg.ny/arg.ns, arg.ns);

ksxp = reshape(ksp(:,1),arg.nx,[],arg.ns);
ksyp = reshape(ksp(:,2),arg.nx,[],arg.ns);
ksxn = reshape(ksn(:,1),arg.nx,[],arg.ns);
ksyn = reshape(ksn(:,2),arg.nx,[],arg.ns);

if arg.indshot && arg.indcoil
	yk(:,1:2:end,1,:) = yp(:,:,1,:);
	
	ksx(:,1:2:end,1) = ksxp(:,:,1);
	ksy(:,1:2:end,1) = ksyp(:,:,1);
	
	for ii=1:arg.ns
		for ic=1:arg.nc
			y3 = reshape(yp(:,:,ii,ic),arg.nx,[]);
			f2 = fftshift(ifft(ifftshift(y3,1)),1);
			f2 = f2.*exp(1i*phi(1,ii,ic));
			y3 = fftshift(fft(ifftshift(f2,1)),1);
			yk(:,1:2:end,ii,ic) = reshape(y3,arg.nx,[]);
		end
		
		ksx(:,1:2:end,ii) = ksxp(:,:,ii) + dkf(1,ii,1);
		ksy(:,1:2:end,ii) = ksyp(:,:,ii) + dkf(2,ii,1);
	end
	
	for ii=1:arg.ns
		for ic=1:arg.nc
			y3 = flipud(reshape(yn(:,:,ii,ic),arg.nx,[]));
			f2 = fftshift(ifft(ifftshift(y3,1)),1);
			f2 = f2.*exp(1i*phi(1,arg.ns + ii,ic));
			y3 = fftshift(fft(ifftshift(f2,1)),1);
			yk(:,2:2:end,ii,ic) = reshape(flipud(y3),arg.nx,[]);
		end
		
		ksx(:,2:2:end,ii) = ksxn(:,:,ii) + dkf(1,arg.ns + ii,1);
		ksy(:,2:2:end,ii) = ksyn(:,:,ii) + dkf(2,arg.ns + ii,1);
	end
elseif arg.indshot
	yk(:,1:2:end,1,:) = yp(:,:,1,:);
	
	ksx(:,1:2:end,1) = ksxp(:,:,1);
	ksy(:,1:2:end,1) = ksyp(:,:,1);
	
	for ii=1:arg.ns
		y3 = reshape(yp(:,:,ii,:),arg.nx,[]);
		f2 = fftshift(ifft(ifftshift(y3,1)),1);
		f2 = f2.*exp(1i*phi(1,ii));
		y3 = fftshift(fft(ifftshift(f2,1)),1);
		yk(:,1:2:end,ii,:) = reshape(y3,arg.nx,[],1,arg.nc);
		
		ksx(:,1:2:end,ii) = ksxp(:,:,ii) + dkf(1,ii);
		ksy(:,1:2:end,ii) = ksyp(:,:,ii) + dkf(2,ii);
	end
	
	for ii=1:arg.ns
		y3 = flipud(reshape(yn(:,:,ii,:),arg.nx,[]));
		f2 = fftshift(ifft(ifftshift(y3,1)),1);
		f2 = f2.*exp(1i*phi(1,arg.ns + ii));
		y3 = fftshift(fft(ifftshift(f2,1)),1);
		yk(:,2:2:end,ii,:) = reshape(flipud(y3),arg.nx,[],1,arg.nc);
		
		ksx(:,2:2:end,ii) = ksxn(:,:,ii) + dkf(1,arg.ns + ii);
		ksy(:,2:2:end,ii) = ksyn(:,:,ii) + dkf(2,arg.ns + ii);
	end
	
elseif arg.indcoil
	yk(:,1:2:end,:,:) = yp;
	for ic=1:arg.nc
		y3 = flipud(reshape(yn(:,:,:,ic),arg.nx,[]));
		f2 = fftshift(ifft(ifftshift(y3,1)),1);
		f2 = f2.*exp(1i*phi(1,ic));
		y3 = fftshift(fft(ifftshift(f2,1)),1);
		yk(:,2:2:end,:,ic) = reshape(flipud(y3),arg.nx,[],arg.ns);
	end
	
	ksx(:,1:2:end,:) = ksxp;
	ksy(:,1:2:end,:) = ksyp;
	
	ksx(:,2:2:end,:) = ksxn + dkf(1);
	ksy(:,2:2:end,:) = ksyn + dkf(2);
else
	yk(:,1:2:end,:,:) = yp;
	y3 = flipud(reshape(yn,arg.nx,[]));
	f2 = fftshift(ifft(ifftshift(y3,1)),1);
	f2 = f2.*exp(1i*phi(1));
	y3 = fftshift(fft(ifftshift(f2,1)),1);
	yk(:,2:2:end,:,:) = reshape(flipud(y3),arg.nx,[],arg.ns,arg.nc);
	
	ksx(:,1:2:end,:) = ksxp;
	ksy(:,1:2:end,:) = ksyp;
	
	ksx(:,2:2:end,:) = ksxn + dkf(1);
	ksy(:,2:2:end,:) = ksyn + dkf(2);
end

ks = [ksx(:) ksy(:)];

end

function yc = data_correction(yp, yn, phi, arg)


yc = zeros(arg.nx, arg.ny/arg.ns, arg.ns, arg.nc);

if arg.indshot && arg.indcoil
	yc(:,1:2:end,1,:) = yp(:,:,1,:);
	
	for ii=1:arg.ns
		for ic=1:arg.nc
			fitvec = phi(1,ii,ic) + phi(2,ii,ic)*arg.xx;
			
			yr2 = reshape(yp(:,:,ii,ic),arg.nx,[]);
			f2 = fftshift(ifft(ifftshift(yr2,1)),1);
			f2 = f2.*repmat(exp(1i*fitvec),[1 arg.ny/arg.ns/2]);
			yr2 = fftshift(fft(ifftshift(f2,1)),1);
			yc(:,1:2:end,ii,ic) = reshape(yr2,arg.nx,[]);
		end
	end
	
	for ii=1:arg.ns
		for ic=1:arg.nc
			fitvec = phi(1,arg.ns+ii,ic) + phi(2,arg.ns+ii,ic)*arg.xx;
			
			yr2 = flipud(reshape(yn(:,:,ii,ic),arg.nx,[]));
			f2 = fftshift(ifft(ifftshift(yr2,1)),1);
			f2 = f2.*repmat(exp(1i*fitvec),[1 arg.ny/arg.ns/2]);
			yr2 = fftshift(fft(ifftshift(f2,1)),1);
			yc(:,2:2:end,ii,ic) = reshape(flipud(yr2),arg.nx,[]);
		end
	end
elseif arg.indshot
	yc(:,1:2:end,1,:) = yp(:,:,1,:);
	
	for ii=1:arg.ns
		fitvec = phi(1,ii) + phi(2,ii)*arg.xx;
		
		yr2 = reshape(yp(:,:,ii,:),arg.nx,[]);
		f2 = fftshift(ifft(ifftshift(yr2,1)),1);
		f2 = f2.*repmat(exp(1i*fitvec),[1 arg.ny/arg.ns/2*arg.nc]);
		yr2 = fftshift(fft(ifftshift(f2,1)),1);
		yc(:,1:2:end,ii,:) = reshape(yr2,arg.nx,[],1,arg.nc);
	end
	
	for ii=1:arg.ns
		fitvec = phi(1,arg.ns+ii) + phi(2,arg.ns+ii)*arg.xx;
		
		yr2 = flipud(reshape(yn(:,:,ii,:),arg.nx,[]));
		f2 = fftshift(ifft(ifftshift(yr2,1)),1);
		f2 = f2.*repmat(exp(1i*fitvec),[1 arg.ny/arg.ns/2*arg.nc]);
		yr2 = fftshift(fft(ifftshift(f2,1)),1);
		yc(:,2:2:end,ii,:) = reshape(flipud(yr2),arg.nx,[],1,arg.nc);
	end
elseif arg.indcoil
	yc(:,1:2:end,:,:) = yp;
	
	for ic=1:arg.nc
		fitvec = phi(1,ic) + phi(2,ic)*arg.xx;
		
		yr2 = flipud(reshape(yn(:,:,:,ic),arg.nx,[]));
		f2 = fftshift(ifft(ifftshift(yr2,1)),1);
		f2 = f2.*repmat(exp(1i*fitvec),[1 arg.ny/2]);
		yr2 = fftshift(fft(ifftshift(f2,1)),1);
		yc(:,2:2:end,:,ic) = reshape(flipud(yr2),arg.nx,[],arg.ns);
	end
else
	yc(:,1:2:end,:,:) = yp;
	
	fitvec = phi(1) + phi(2)*arg.xx;
	
	yr2 = flipud(reshape(yn,arg.nx,[]));
	f2 = fftshift(ifft(ifftshift(yr2,1)),1);
	f2 = f2.*repmat(exp(1i*fitvec),[1 arg.ny/2*arg.nc]);
	yr2 = fftshift(fft(ifftshift(f2,1)),1);
	yc(:,2:2:end,:,:) = reshape(flipud(yr2),arg.nx,[],arg.ns,arg.nc);
end



end

function out = plot_results(yi, yp, yn, yc, yk, ks, ksc, ksp, ksn, smap, pmap, arg)

arg.niter = 50;
nufft = {arg.N, [6 6], 2*arg.N, arg.N/2, 'table',2^10,'minmax:kb'};

Gmk = Gmri(ksc, arg.mask, 'fov', arg.fov, 'basis', {arg.basis},'nufft', nufft);

nufft = {arg.N, [1 1], arg.N, arg.N/2, 'linear'};

Gm = Gmri(ks, arg.mask, 'fov', arg.fov, 'basis', {arg.basis},'nufft', nufft);
Gmp = Gmri(ksp, arg.mask, 'fov', arg.fov, 'basis', {arg.basis},'nufft', nufft);
Gmn = Gmri(ksn, arg.mask, 'fov', arg.fov, 'basis', {arg.basis},'nufft', nufft);

if ~isempty(arg.fmap)
	Gm = feval(Gm.arg.new_zmap, Gm, arg.ti(:), arg.fmap, arg.tseg);
	Gmp = feval(Gmp.arg.new_zmap, Gmp, arg.tip(:), arg.fmap, arg.tseg);
	Gmn = feval(Gmp.arg.new_zmap, Gmn, arg.tin(:), arg.fmap, arg.tseg);
	Gmk = feval(Gmk.arg.new_zmap, Gmk, arg.ti(:), arg.fmap, arg.tseg);
end

xxk = zeros(arg.nx,arg.ny,arg.nc);
xxm = zeros(arg.nx,arg.ny,arg.nc);
xxi = zeros(arg.nx,arg.ny,arg.nc);
xinit = zeros(arg.N);

Gc = cell(1,arg.nc);
Gck = cell(1,arg.nc);
Gcn = cell(arg.nc,1);
Gcp = cell(arg.nc,1);
for ii=1:arg.nc
	if arg.indcoil
		tmp = exp(-1i*pmap(:,:,ii));
	else
		tmp = exp(-1i*pmap);
	end
	
	smp = smap(:,:,ii);
	
	Gc{ii} = Gm * Gdiag(smp(arg.mask),'mask',arg.mask);
	Gck{ii} = Gmk * Gdiag(smp(arg.mask),'mask',arg.mask);
	Gcp{ii} = Gmp * Gdiag(smp(arg.mask),'mask',arg.mask);
	Gcn{ii} = Gmn * Gdiag(smp(arg.mask).*tmp(arg.mask),'mask',arg.mask);
	
	Gmm = block_fatrix({Gcp{ii}, Gcn{ii}}, 'type', 'col');
	
% 	ticker reset, ticker print, ticker(4);
% 	xi = qpwls_pcg1(xinit(arg.mask), Gck{ii}, 1, vec(yk(:,:,:,ii)), 0, ...
% 		'niter', arg.niter);
% 	xxk(:,:,ii) = embed(xi(:,end), arg.mask);
% 	
% 	ticker reset, ticker print, ticker(4);
% 	ym = [vec(yp(:,:,:,ii));vec(yn(:,:,:,ii))];
% 	xi = qpwls_pcg1(xinit(arg.mask), Gmm, 1, ym, 0, 'niter', arg.niter);
% 	xxm(:,:,ii) = embed(xi(:,end), arg.mask);
% 	
% 	ticker reset, ticker print, ticker(4);
% 	xi = qpwls_pcg1(xinit(arg.mask), Gc{ii}, 1, vec(yc(:,:,:,ii)), 0, ...
% 		'niter', arg.niter);
% 	xxi(:,:,ii) = embed(xi(:,end), arg.mask);
end

Gb = block_fatrix(Gc, 'type', 'col');
Gbp = block_fatrix(Gcp, 'type', 'col');
Gbn = block_fatrix(Gcn, 'type', 'col');
Gbm = block_fatrix({Gbp, Gbn}, 'type', 'col');
Gbk = block_fatrix(Gck, 'type', 'col');

ticker reset, ticker print, ticker(4);
xi = qpwls_pcg1(xinit(arg.mask), Gbk, 1, yk(:), 0,'niter', arg.niter);
xf = embed(xi(:,end), arg.mask);

ticker reset, ticker print, ticker(4);
xi = qpwls_pcg1(xinit(arg.mask), Gbm, 1, [yp(:);yn(:)], 0,'niter', arg.niter);
xm = embed(xi(:,end), arg.mask);

ticker reset, ticker print, ticker(4);
xi = qpwls_pcg1(xinit(arg.mask), Gb, 1, yc(:), 0,'niter', arg.niter);
xd = embed(xi(:,end), arg.mask);

ticker reset, ticker print, ticker(4);
xi = qpwls_pcg1(xinit(arg.mask), Gb, 1, yi(:), 0,'niter', arg.niter);
xt = embed(xi(:,end), arg.mask);

out.xcimg = xxi;
out.xcksp = xxk;
out.xcmod = xxm;
out.ximg = xd;
out.xksp = xf;
out.xmod = xm;
out.xunc = xt;

figure(1)
im clf
im pl 2 2
im(1,abs(xt),'cbar','Uncorrected image')
im(2,abs(xf),'cbar','Correction in kspace')
im(3,abs(xd),'cbar','Correction in image domain')
im(4,abs(xm),'cbar','Model based correction')
figure(2)
im(abs([xf;xd;xm]),'cbar','Correction in kspace, image and model based')
figure(3)
im(abs([xxk;xxi;xxm]),'cbar','Correction for each coil');

end

function out = plot_results_shot(yi, yp, yn, yc, yk, ks, ksc, smap, pmap, arg)

arg.niter = 50;
nufft = {arg.N, [6 6], 2*arg.N, arg.N/2, 'table',2^10,'minmax:kb'};

Gmk = Gmri(ksc, arg.mask, 'fov', arg.fov, 'basis', {arg.basis},'nufft', nufft);

nufft = {arg.N, [1 1], arg.N, arg.N/2, 'linear'};

Gm = Gmri(ks, arg.mask, 'fov', arg.fov, 'basis', {arg.basis},'nufft', nufft);

[Gp, Gn] = create_system_shot_coil(yi, ks, smap, arg);

if ~isempty(arg.fmap)
	Gm = feval(Gm.arg.new_zmap, Gm, arg.ti(:), arg.fmap, arg.tseg);
	Gmk = feval(Gmk.arg.new_zmap, Gmk, arg.ti(:), arg.fmap, arg.tseg);
end

xxk = zeros(arg.nx,arg.ny,arg.nc);
xxm = zeros(arg.nx,arg.ny,arg.nc);
xxi = zeros(arg.nx,arg.ny,arg.nc);
xinit = zeros(arg.N);

Gc = cell(1,arg.nc);
Gck = cell(1,arg.nc);
Gcn = cell(arg.nc,1);
Gcp = cell(arg.nc,1);
for ii=1:arg.nc
	smp = smap(:,:,ii);
	
	Gsn = cell(1,arg.ns);
	Gsp = cell(1,arg.ns);
	for is=1:arg.ns
		if arg.indcoil
			tmpp = exp(-1i*pmap(:,:,is,ii));
			tmpn = exp(-1i*pmap(:,:,arg.ns+is,ii));
		else
			tmpp = exp(-1i*pmap(:,:,is));
			tmpn = exp(-1i*pmap(:,:,arg.ns+is));
		end
		
		Gsp{is} = Gp{is,ii}*Gdiag(tmpp(arg.mask),'mask',arg.mask);
		Gsn{is} = Gn{is,ii}*Gdiag(tmpn(arg.mask),'mask',arg.mask);
	end
	
	Gcp{ii} = block_fatrix(Gsp, 'type', 'col');
	Gcn{ii} = block_fatrix(Gsn, 'type', 'col');
	Gc{ii} = Gm * Gdiag(smp(arg.mask),'mask',arg.mask);
	Gck{ii} = Gmk * Gdiag(smp(arg.mask),'mask',arg.mask);
	
	Gmm = block_fatrix({Gcp{ii}, Gcn{ii}}, 'type', 'col');
	
% 	ticker reset, ticker print, ticker(4);
% 	xi = qpwls_pcg1(xinit(arg.mask), Gck{ii}, 1, vec(yk(:,:,:,ii)), 0, 'niter', arg.niter);
% 	xxk(:,:,ii) = embed(xi(:,end), arg.mask);
% 	
% 	ticker reset, ticker print, ticker(4);
% 	ym = [vec(yp(:,:,:,ii));vec(yn(:,:,:,ii))];
% 	xi = qpwls_pcg1(xinit(arg.mask), Gmm, 1, ym, 0, 'niter', arg.niter);
% 	xxm(:,:,ii) = embed(xi(:,end), arg.mask);
% 	
% 	ticker reset, ticker print, ticker(4);
% 	xi = qpwls_pcg1(xinit(arg.mask), Gc{ii}, 1, vec(yc(:,:,:,ii)), 0, 'niter', arg.niter);
% 	xxi(:,:,ii) = embed(xi(:,end), arg.mask);
end

Gb = block_fatrix(Gc, 'type', 'col'); % [G1; G2; ... ]
Gbp = block_fatrix(Gcp, 'type', 'col'); % [G1; G2; ... ]
Gbn = block_fatrix(Gcn, 'type', 'col'); % [G1; G2; ... ]
Gbm = block_fatrix({Gbp, Gbn}, 'type', 'col');
Gbk = block_fatrix(Gck, 'type', 'col'); % [G1; G2; ... ]

ticker reset, ticker print, ticker(4);
xi = qpwls_pcg1(xinit(arg.mask), Gbk, 1, yk(:), 0,'niter', arg.niter);
xf = embed(xi(:,end), arg.mask);

ticker reset, ticker print, ticker(4);
xi = qpwls_pcg1(xinit(arg.mask), Gbm, 1, [yp(:);yn(:)], 0,'niter', arg.niter);
xm = embed(xi(:,end), arg.mask);

ticker reset, ticker print, ticker(4);
xi = qpwls_pcg1(xinit(arg.mask), Gb, 1, yc(:), 0,'niter', arg.niter);
xd = embed(xi(:,end), arg.mask);

ticker reset, ticker print, ticker(4);
xi = qpwls_pcg1(xinit(arg.mask), Gb, 1, yi(:), 0,'niter', arg.niter);
xt = embed(xi(:,end), arg.mask);

out.xcimg = xxi;
out.xcksp = xxk;
out.xcmod = xxm;
out.ximg = xd;
out.xksp = xf;
out.xmod = xm;
out.xunc = xt;

figure(1)
im pl 2 2
im(1,abs(xt),'cbar','Uncorrected image')
im(2,abs(xf),'cbar','Correction in kspace')
im(3,abs(xd),'cbar','Correction in image domain')
im(4,abs(xm),'cbar','Model based correction')
figure(2)
im(abs([xf;xd;xm]),'cbar','Correction in kspace, image and model based')
figure(3)
im(abs([xxk;xxi;xxm]),'cbar','Correction for each coil');

end


