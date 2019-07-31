 function [xr, out] = image_restore(y, psf, varargin)
%function [xr, out] = image_restore(y, psf, varargin)
%|
%| Function to perform image restoration using Augmented Lagrangian methods and
%| "Split Bregman" approach (single split of regularizer term).
%|
%| in
%|	y	[nx ny]		Input noisy and/or blurred image
%|	psf	[nbx nby]	PSF of blur
%|
%| options
%|	est_type	string	Type of restoration. Valid choices
%|		'ncg'		Non-linear CG algorithm with rounding parameter
%|							and line search method
%|		'ista'		ISTA method with Chambolle type algorithm
%|		'fista'		FISTA method with inner Chambolle iterations
%|		'mfista'	MFISTA methode with inner Chambolle iterations
%|		'salsa'		SALSA algorithm with inner Chambolle iterations
%|		'alp1'		AL method with one split in regularization term
%|		'alp2'		AL method with two splits, one for regularization
%|						and one for data-fit term
%|		'cgq'		Conjugate gradient reconstruction using quadratic
%|						regularization
%|		'alq'		AL reconstruction using quadratic regularization
%|
%|	model_type	string	Type of reconstruction model
%|		'circ'			Unrealistic circulant model
%|		'circ-edge'		Circulant model with edge-taper
%|		'circ-zpad'		Circulant model with 0 padding of boundaries
%|		'circ-zpapedge'		Circulant model with 0 padding and edge taper
%|		'circ-rpad'		Circulant model with 0 order boundary extension
%|		'circ-rpadedge'		Circulant model with 0 order boundary
%|						 extension and edge taper
%|		'circ-ppad'		Circulant model with periodic padding of boundaries
%|		'circ-ppapedge'		Circulant model with periodic padding and edge taper
%|		'circ-spad'		Circulant model with symmetric padding of boundaries
%|		'circ-spapedge'		Circulant model with symmetric padding and edge taper
%|		'non-circ'		Accurate non-circulant model
%|
%|	reg_type	string	Type of regularizer
%|		'l1,wavelet'		l1 penalty with wavelet transform
%|		'tv,iso'			Isotropic TV
%|		'tv,aniso'			Anisotropic TV
%|
%|	psf_type	string		type of Gblur to be used
%|		'fft,same'			Circulant blur using fftn
%|		'conv,same'			Toeplitz blur using convn
%|
%|	wave_type	string		Type of wavelet for wavelet l1 norm regularization
%|							Valid parameters from wfilters
%|	lambda		scalar		regularization parameter
%|	mu			scalar		Lagrangian tuning parameter
%|	nu			[2 1]		Lagrangian tuning parameter
%|	xtrue		[nx ny]		True image for reference
%|	niter		scalar		# of outer iterations
%|	chiter		scaler		# of inner Chambolle iterations
%|	cgiter		scalar		# of ineer CG iterations
%|	lniter		scalar		# of inner line search iterations
%|	cgr_thres	scalar		# of iter to apply CG in SB-MIL method
%|	plot		0|1			plot results
%|	chat		0|1			show printouts
%|
%| Copyright 2011-10-15, Antonis Matakos, University of Michigan
%| Last revision 2012-01-04

if nargin == 1 && streq(y, 'test'), image_restore_test; return; end
if nargin < 2, help(mfilename), error(mfilename), end

%% General Setup

arg.est_type = 'alp2';
arg.model_type = 'non-circ';
arg.reg_type = 'tv,aniso';
arg.psf_type = 'fft,same';
arg.wave_type = 'haar';
arg.lambda = 2^-5;
arg.mu = 2^-1;
arg.nu = 2^-3/10;
arg.xtrue = [];
arg.niter = 500;
arg.chiter = 10;
arg.cgiter = 10;
arg.lniter = 10;
arg.cgr_thres = inf;
arg.cg_tol = 1e-6; % Tolerance of inner CG - hidden parameter (OK)
arg.epsilon = 1e-8; % Rounding parameter for NCG - hidden parameter
arg.precon = 1; % Preconditioner use for inner CG - hidden parameter (YES)
arg.L = 1; % Lipschitz constant - hidden parameter (Leave unchanged)
arg.extra = 0; % Extra update of eta params - hidden parameter (NO)
arg.skip = 50;
arg.plot = 0;
arg.chat = 0;

arg = vararg_pair(arg, varargin);

arg.clim = [0 max(arg.xtrue(:)) + 1];

[nx ny] = size(y);
[nbx nby] = size(psf);

if strcmpi(arg.model_type, 'non-circ')
	imask = true(nx+nbx-1, ny+nby-1);
	mask = false(nx+nbx-1, ny+nby-1);
	mask((nbx+1)/2:end-(nbx-1)/2, (nby+1)/2:end-(nby-1)/2) = true;
	dmask = true(nx,ny);
	omask = mask;
	arg.trim = 0;

elseif strcmpi(arg.model_type, 'circ-edge')
	mask = true(nx,ny);
	imask = mask;
	dmask = mask;
	omask = mask;
	y = edgetaper(y, psf);
	arg.trim = 0;

elseif strcmpi(arg.model_type, 'circ-rpad')
	mask = true(nx+nbx-1,ny+nby-1);
	imask = mask;
	dmask = false(nx+nbx-1,ny+nby-1);
	dmask((nbx+1)/2:end-(nbx-1)/2,(nby+1)/2:end-(nby-1)/2) = true;
	omask = dmask;
	y = padarray(y,[(nbx-1)/2 (nby-1)/2],'replicate');
	arg.trim = 1;

elseif strcmpi(arg.model_type, 'circ-rpadedge')
	mask = true(nx+nbx-1,ny+nby-1);
	imask = mask;
	dmask = false(nx+nbx-1,ny+nby-1);
	dmask((nbx+1)/2:end-(nbx-1)/2,(nby+1)/2:end-(nby-1)/2) = true;
	omask = dmask;
	y = padarray(y,[(nbx-1)/2 (nby-1)/2],'replicate');
	y = edgetaper(y, psf);
	arg.trim = 1;
%	mask = true(nx+2*nbx,nx+2*nby);
%	imask = mask;
%	dmask = false(nx+2*nbx,nx+2*nby);
%	dmask(nbx+1:end-nbx,nby+1:end-nby) = true;
%	omask = dmask;
%	y = padarray(y,[nbx nby],'replicate');
%	y = edgetaper(y, psf);

elseif strcmpi(arg.model_type, 'circ-zpad')
	mask = true(nx+nbx-1,ny+nby-1);
	imask = mask;
	dmask = false(nx+nbx-1,ny+nby-1);
	dmask((nbx+1)/2:end-(nbx-1)/2,(nby+1)/2:end-(nby-1)/2) = true;
	omask = dmask;
	y = padarray(y,[(nbx-1)/2 (nby-1)/2],0);
	arg.trim = 1;

elseif strcmpi(arg.model_type, 'circ-zpadedge')
	mask = true(nx+nbx-1,ny+nby-1);
	imask = mask;
	dmask = false(nx+nbx-1,ny+nby-1);
	dmask((nbx+1)/2:end-(nbx-1)/2,(nby+1)/2:end-(nby-1)/2) = true;
	omask = dmask;
	y = padarray(y,[(nbx-1)/2 (nby-1)/2],0);
	y = edgetaper(y, psf);
	arg.trim = 1;

elseif strcmpi(arg.model_type, 'circ-ppad')
	mask = true(nx+nbx-1,ny+nby-1);
	imask = mask;
	dmask = false(nx+nbx-1,ny+nby-1);
	dmask((nbx+1)/2:end-(nbx-1)/2,(nby+1)/2:end-(nby-1)/2) = true;
	omask = dmask;
	y = padarray(y,[(nbx-1)/2 (nby-1)/2],'circular');
	arg.trim = 1;

elseif strcmpi(arg.model_type, 'circ-ppadedge')
	mask = true(nx+nbx-1,ny+nby-1);
	imask = mask;
	dmask = false(nx+nbx-1,ny+nby-1);
	dmask((nbx+1)/2:end-(nbx-1)/2,(nby+1)/2:end-(nby-1)/2) = true;
	omask = dmask;
	y = padarray(y,[(nbx-1)/2 (nby-1)/2],'circular');
	y = edgetaper(y, psf);
	arg.trim = 1;

elseif strcmpi(arg.model_type, 'circ-spad')
	mask = true(nx+nbx-1,ny+nby-1);
	imask = mask;
	dmask = false(nx+nbx-1,ny+nby-1);
	dmask((nbx+1)/2:end-(nbx-1)/2,(nby+1)/2:end-(nby-1)/2) = true;
	omask = dmask;
	y = padarray(y,[(nbx-1)/2 (nby-1)/2],'symmetric');
	arg.trim = 1;

elseif strcmpi(arg.model_type, 'circ-spadedge')
	mask = true(nx+nbx-1,ny+nby-1);
	imask = mask;
	dmask = false(nx+nbx-1,ny+nby-1);
	dmask((nbx+1)/2:end-(nbx-1)/2,(nby+1)/2:end-(nby-1)/2) = true;
	omask = dmask;
	y = padarray(y,[(nbx-1)/2 (nby-1)/2],'symmetric');
	y = edgetaper(y, psf);
	arg.trim = 1;

else
	mask = true(nx,ny);
	imask = mask;
	dmask = mask;
	omask = mask;
	arg.trim = 0;
end

A = Gblur(imask, 'psf', psf, 'type', arg.psf_type);

if strcmpi(arg.reg_type,'l1,wave')
	C = Gwave2('mask', imask, 'nlevel', 2, 'wname', arg.wave_type);
	arg.creg = 1;
elseif strcmpi(arg.reg_type,'tv,iso') || strcmpi(arg.reg_type,'tv,aniso')
	C = Cdiffs(size(mask), 'offsets', [1 size(mask,1)],'type_diff','circshift');
	arg.creg = 8;
else
	error('Unknown regularization option');
end

cfft = double(fft2(reshape(C(:,1),C.odim)));

if strcmpi(arg.reg_type,'tv,iso')
	arg.is_iso = 1;
else
	arg.is_iso = 0;
end

yp = padarray(y,[(nbx-1)/2 (nby-1)/2],'replicate');
yp = edgetaper(yp, psf);

if strcmpi(arg.est_type, 'alp2')
	[xr isnr err time itrtime] = restore_alp2(y, yp, mask, cfft, A, C, arg);
elseif strcmpi(arg.est_type, 'alp1')
	[xr isnr err time itrtime] = restore_alp1(y, mask, dmask, omask, cfft, A, C, arg);
elseif strcmpi(arg.est_type, 'reeves')
	[xr isnr err time itrtime] = restore_reeves(y, yp, mask, dmask, omask, cfft, A, C, arg);
elseif strcmpi(arg.est_type, 'salsa')
	[xr isnr err time itrtime] = restore_salsa(y, mask, dmask, omask, cfft, A, C, arg);
elseif strcmpi(arg.est_type, 'mfista') && strcmpi(arg.model_type, 'refl-dct')
	[xr isnr err time itrtime] = restore_mfista_dct(y, mask, dmask, omask, cfft, A, C, arg);
elseif strcmpi(arg.est_type, 'mfista')
	[xr isnr err time itrtime] = restore_mfista(y, mask, dmask, omask, cfft, A, C, arg);
elseif strcmpi(arg.est_type, 'fista')
	[xr isnr err time itrtime] = restore_fista(y, mask, cfft, A, C, arg);
elseif strcmpi(arg.est_type, 'ista')
	[xr isnr err time itrtime] = restore_ista(y, mask, cfft, A, C, arg);
elseif strcmpi(arg.est_type, 'ncg')
	[xr isnr err time itrtime] = restore_ncg(y, mask, cfft, A, C, arg);
elseif strcmpi(arg.est_type, 'cgq')
	[xr isnr err time itrtime] = restore_cg_quad(y, mask, cfft, A, C, arg);
elseif strcmpi(arg.est_type, 'alq')
	[xr isnr err time itrtime] = restore_al_quad(y, mask, cfft, A, C, arg);
else
	error('Unknown estimation method');
end

xr = reale(xr,'warn','non-negligible imaginary part');
if arg.trim
	xr = xr((nbx+1)/2:end-(nbx-1)/2,(nby+1)/2:end-(nby-1)/2);
%	xr = xr(nbx+1:end-nbx,nby+1:end-nby);
end

out.isnr = isnr;
out.err = err;
out.time = time;
out.psf = psf;

out.cputime_iter = time(end,1)/arg.niter;
out.time_iter = time(end,2)/arg.niter;

if strcmpi(arg.est_type, 'salsa')
	out.cputime_cgiter = mean(itrtime(:,1));
	out.time_cgiter = mean(itrtime(:,2));
	out.cputime_chiter = mean(itrtime(:,3));
	out.time_chiter = mean(itrtime(:,4));
	out.cputime_lniter = 0;
	out.time_lniter = 0;
elseif strcmpi(arg.est_type, 'alp1') || strcmpi(arg.est_type, 'reeves')
	out.cputime_cgiter = mean(itrtime(:,1));
	out.time_cgiter = mean(itrtime(:,2));
	out.cputime_chiter = 0;
	out.time_chiter = 0;
	out.cputime_lniter = 0;
	out.time_lniter = 0;
elseif strcmpi(arg.est_type, 'fista') || strcmpi(arg.est_type, 'mfista')
	out.cputime_cgiter = 0;
	out.time_cgiter = 0;
	out.cputime_chiter = mean(itrtime(:,1));
	out.time_chiter = mean(itrtime(:,2));
	out.cputime_lniter = 0;
	out.time_lniter = 0;
elseif strcmpi(arg.est_type, 'ncg')
	out.cputime_cgiter = 0;
	out.time_cgiter = 0;
	out.cputime_chiter = 0;
	out.time_chiter = 0;
	out.cputime_lniter = mean(itrtime(:,1));
	out.time_lniter = mean(itrtime(:,2));
else
	out.cputime_cgiter = 0;
	out.time_cgiter = 0;
	out.cputime_chiter = 0;
	out.time_chiter = 0;
	out.cputime_lniter = 0;
	out.time_lniter = 0;
end

end

%% Proposed algorithms ADMM-P2
function [xr isnr err time itrtime] = restore_alp2(y, yp, mask, cfft, A, C, arg)

if strcmpi(arg.reg_type, 'l1,wave')
	arg.do_cfft = 1;
else
	arg.do_cfft = 0;
end

bfft = A.arg.psf_fft;

Hi = abs(bfft).^2 + arg.nu*sum(abs(cfft).^2,3);
Hd = 1./(double(mask)+arg.mu);

printm('Cond number of Hnu: %f', max(Hi(:))/min(Hi(:)))
printm('Cond number of Hmu: %f', (1+arg.mu)/arg.mu);

% Next 2 lines possibly redundant. Keep for safety
T = Gdiag2(ones(size(y)), 'mask_in', mask);
d = T'*y;

% xr = A'*(T'*y);
xr = yp;
Axr = A*xr;
Cxr = C*xr;

eta0 = zeros(C.idim);
eta1 = zeros(C.odim);

isnr = zeros(arg.niter+1,1);
err = zeros(arg.niter+1,1);
time = zeros(arg.niter+1,2);
itrtime = zeros(arg.niter+1,2);

isnr(1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xr(mask)-arg.xtrue(:)));
err(1) = 20*log10(nrms(xr(mask), arg.xtrue(:)));

tm1 = cputime;
tm2 = tic;

for ii=1:arg.niter

	if mod(ii,arg.skip) == 0 && arg.chat
		printm('iter: %d',ii)
	end

	u0 = Hd.*(d + arg.mu*(Axr + eta0));
	if arg.extra
		eta0 = eta0 - u0 + Axr;
	end

	if arg.is_iso
		u1 = shrink(Cxr + eta1, arg.lambda/arg.mu/arg.nu, 'vector');
	else
		u1 = shrink(Cxr + eta1, arg.lambda/arg.mu/arg.nu);
	end
	if arg.extra
		eta1 = eta1 - u1 + Cxr;
	end

	if arg.do_cfft
		xr = real(ifft2((conj(bfft).*fft2(u0-eta0) + ...
			arg.nu*sum(conj(cfft).*fft2(u1-eta1),3))./Hi));
	else
		xr = real(ifft2((conj(bfft).*fft2(u0-eta0) + arg.nu*fft2(C'*(u1-eta1)))./Hi));
	end

	Axr = A*xr;
	Cxr = C*xr;

	eta0 = eta0 - u0 + Axr;
	eta1 = eta1 - u1 + Cxr;

	time(ii+1,1) = cputime - tm1;
	time(ii+1,2) = toc(tm2);
	err(ii+1) = 20*log10(nrms(xr(mask), arg.xtrue(:)));
	isnr(ii+1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xr(mask)-arg.xtrue(:)));
end

if arg.plot
	figure(1), im(arg.xtrue,'Original image','cbar',arg.clim);
	figure(2), im(y,'Blured noisy image','cbar',arg.clim);
	figure(3), im(xr,'Restored image','cbar',arg.clim);
	figure(4), im(xr-arg.xtrue,'Difference image','cbar');
end

if arg.chat
	printm(['\nCPU time: %6.2f\tActual time: %6.2f\nNRMSE (db): %6.2f\t' ...
		'NRMSE: %6.4f%%\nISNR: %6.4f\n'], time(end,1), time(end,2), err(end),...
		100*nrms(xr(mask), arg.xtrue(:)), isnr(end));
end

end


%% Split-Bregman algorithm
function [xr isnr err time itrtime] = restore_alp1(y, mask, dmask, omask, cfft, A, C, arg)

if strcmpi(arg.model_type, 'non-circ')
	arg.do_cg = 1;
else
	arg.do_cg = 0;
end

if strcmpi(arg.reg_type, 'l1,wave')
	arg.do_cfft = 1;
else
	arg.do_cfft = 0;
end

T = Gdiag2(ones(size(y)), 'mask_in', mask);
bfft = A.arg.psf_fft;

dfft = conj(bfft) .* fft2(T' * y);
d = A' * (T' * y);
Hi = abs(bfft).^2 + arg.mu * sum(abs(cfft).^2, 3);

printm('Cond number of Hmu: %f', max(Hi(:))/min(Hi(:)))

xr = d;
Cxr = C*xr;
eta = zeros(C.odim);

isnr = zeros(arg.niter+1,1);
err = zeros(arg.niter+1,1);
time = zeros(arg.niter+1,2);
itrtime = zeros(arg.niter+1,2);

isnr(1) = 20*log10(norm(y(dmask)-arg.xtrue(:))/norm(xr(omask)-arg.xtrue(:)));
err(1) = 20*log10(nrms(xr(omask), arg.xtrue(:)));

tm1 = cputime;
tm2 = tic;

for ii=1:arg.niter

	if mod(ii,50) == 0 && arg.chat
		printm('iter: %d',ii)
	end

	if arg.is_iso
		u = shrink(Cxr + eta, arg.lambda/arg.mu, 'vector');
	else
		u = shrink(Cxr + eta, arg.lambda/arg.mu);
	end
	if arg.extra
		eta = eta - u + Cxr;
	end

	if arg.do_cg
		yd = d + arg.mu*(C'*(u-eta));
		[xr tm] = inner_cg(yd, xr, A, T, C, 1./Hi, cfft, arg);
		itrtime(ii+1,1) = tm(1);
		itrtime(ii+1,2) = tm(2);
	else
		if arg.do_cfft
			xr = real(ifft2((dfft + arg.mu*sum(conj(cfft).*fft2(u-eta),3))./Hi));
		else
			xr = real(ifft2((dfft + arg.mu*fft2(C'*(u-eta)))./Hi));
		end
	end

	Cxr = C*xr;
	eta = eta - u + Cxr;

	time(ii+1,1) = cputime - tm1;
	time(ii+1,2) = toc(tm2);
	err(ii+1) = 20*log10(nrms(xr(omask), arg.xtrue));
	isnr(ii+1) = 20*log10(norm(y(dmask)-arg.xtrue(:))/norm(xr(omask)-arg.xtrue(:)));

end

if arg.plot
	figure(1), im(arg.xtrue,'Original image','cbar',arg.clim);
	figure(2), im(y,'Blured noisy image','cbar',arg.clim);
	figure(3), im(xr,'Restored image','cbar',arg.clim);
	figure(4), im(xr-arg.xtrue,'Difference image','cbar');
end

if arg.chat
	printm(['\nCPU time: %6.2f\tActual time: %6.2f\nNRMSE (db): %6.2f\t' ...
		'NRMSE: %6.4f%%\nISNR: %6.4f\n'], time(end,1), time(end,2), err(end),...
		100*nrms(xr(omask), arg.xtrue(:)), isnr(end));
end

end


%% Bregman-MIL algorithm
function [xr isnr err time itrtime] = restore_reeves(y, yp, mask, dmask, omask, cfft, A, C, arg)

if strcmpi(arg.model_type, 'non-circ')
	arg.do_cg = 1;
else
	arg.do_cg = 0;
end

if strcmpi(arg.reg_type, 'l1,wave')
	arg.do_cfft = 1;
else
	arg.do_cfft = 0;
end

T = Gdiag2(ones(size(y)),'mask_in',mask);
bfft = A.arg.psf_fft;

dfft = conj(bfft).*fft2(T'*y);
d = A'*(T'*y);
Hi = abs(bfft).^2 + arg.mu*sum(abs(cfft).^2,3);
Ac = A(~mask(:),:);

printm('Cond number of Hmu: %f', max(Hi(:))/min(Hi(:)))

% xr = d;
xr = yp;
% wr = 0*yp(~mask);
wr = yp(~mask);
% wr = d(~mask);
Cxr = C*xr;
eta = zeros(C.odim);

isnr = zeros(arg.niter+1,1);
err = zeros(arg.niter+1,1);
time = zeros(arg.niter+1,2);
itrtime = zeros(arg.niter+1,2);

isnr(1) = 20*log10(norm(y(dmask)-arg.xtrue(:))/norm(xr(omask)-arg.xtrue(:)));
err(1) = 20*log10(nrms(xr(omask), arg.xtrue(:)));

tm1 = cputime;
tm2 = tic;

for ii=1:arg.niter

	if mod(ii,50) == 0 && arg.chat
		printm('iter: %d',ii)
	end

	if arg.is_iso
		u = shrink(Cxr + eta, arg.lambda/arg.mu, 'vector');
	else
		u = shrink(Cxr + eta, arg.lambda/arg.mu);
	end
	if arg.extra
		eta = eta - u + Cxr;
	end

	if arg.do_cg
		yc = arg.mu*fft2(C'*(u-eta));
		if ii <= arg.cgr_thres
			yd = Ac*vec(ifft2((dfft + yc)./Hi));
			[wr tm] = inner_cg_reeves(yd, wr, Ac, 1./Hi, size(yp), arg);
			itrtime(ii+1,1) = tm(1);
			itrtime(ii+1,2) = tm(2);
		end
		yr = (T'*y) + embed(wr, ~mask);
		xr = real(ifft2((conj(bfft).*fft2(yr) + yc)./Hi));

	else
		if arg.do_cfft
			xr = real(ifft2((dfft + arg.mu*sum(conj(cfft).*fft2(u-eta),3))./Hi));
		else
			xr = real(ifft2((dfft + arg.mu*fft2(C'*(u-eta)))./Hi));
		end
	end

	Cxr = C*xr;
	eta = eta - u + Cxr;

	time(ii+1,1) = cputime - tm1;
	time(ii+1,2) = toc(tm2);
	err(ii+1) = 20*log10(nrms(xr(omask), arg.xtrue));
	isnr(ii+1) = 20*log10(norm(y(dmask)-arg.xtrue(:))/norm(xr(omask)-arg.xtrue(:)));

end

if arg.plot
	figure(1), im(arg.xtrue,'Original image','cbar',arg.clim);
	figure(2), im(y,'Blured noisy image','cbar',arg.clim);
	figure(3), im(xr,'Restored image','cbar',arg.clim);
	figure(4), im(xr-arg.xtrue,'Difference image','cbar');
end

if arg.chat
	printm(['\nCPU time: %6.2f\tActual time: %6.2f\nNRMSE (db): %6.2f\t' ...
		'NRMSE: %6.4f%%\nISNR: %6.4f\n'], time(end,1), time(end,2), err(end),...
		100*nrms(xr(omask), arg.xtrue(:)), isnr(end));
end

end


%% SALSA
function [xr isnr err time itrtime] = restore_salsa(y, mask, dmask, omask, cfft, A, C, arg)

if strcmpi(arg.model_type, 'non-circ')
	arg.do_cg = 1;
else
	arg.do_cg = 0;
end

if strcmpi(arg.reg_type, 'l1,wave')
	arg.do_cfft = 1;
else
	arg.do_cfft = 0;
end

T = Gdiag2(ones(size(y)), 'mask_in', mask);
bfft = A.arg.psf_fft;

dfft = conj(bfft).*fft2(T'*y);
d = A'*(T'*y);
Hi = abs(bfft).^2 + arg.mu;

printm('Cond number of Hmu: %f', max(Hi(:))/min(Hi(:)))

xr = d;
u = d;
vk = zeros(C.odim);
eta = 0*d;

isnr = zeros(arg.niter+1,1);
err = zeros(arg.niter+1,1);
time = zeros(arg.niter+1,2);
itrtime = zeros(arg.niter+1,4);

isnr(1) = 20*log10(norm(y(dmask)-arg.xtrue(:))/norm(xr(omask)-arg.xtrue(:)));
err(1) = 20*log10(nrms(xr(omask), arg.xtrue(:)));

tm1 = cputime;
tm2 = tic;

for ii=1:arg.niter

	if mod(ii,50) == 0 && arg.chat
		printm('iter: %d',ii)
	end

	% Inner iterations for updating u
	dk = xr + eta;

	itm1 = cputime;
	itm2 = tic;
	for jj=1:arg.chiter
		Li = abs(C*u)*arg.mu/arg.lambda + arg.creg;
		if arg.do_cfft
			vk = real(arg.creg*vk - C*(ifft2(sum(conj(cfft).*fft2(vk),3)) - dk))./Li;
			u = dk - real(ifft2(sum(conj(cfft).*fft2(vk),3)));
		else
			vk = (arg.creg*vk - C*(C'*vk - dk))./Li;
			u = dk - C'*vk;
		end
	end
	itrtime(ii+1,3) = (cputime - itm1)/arg.chiter;
	itrtime(ii+1,4) = toc(itm2)/arg.chiter;
%	eta = eta - u + Cxr;

	if arg.do_cg
		yd = d + arg.mu*(u-eta);
		[xr tm] = inner_cg(yd, xr, A, T, 1, 1./Hi, cfft, arg); % use 1 for identity in hessian
		itrtime(ii+1,1) = tm(1);
		itrtime(ii+1,2) = tm(2);
	else
		xr = real(ifft2((dfft + arg.mu*fft2(u-eta))./Hi));
	end

	eta = eta - u + xr;

	time(ii+1,1) = cputime - tm1;
	time(ii+1,2) = toc(tm2);
	err(ii+1) = 20*log10(nrms(xr(omask), arg.xtrue(:)));
	isnr(ii+1) = 20*log10(norm(y(dmask)-arg.xtrue(:))/norm(xr(omask)-arg.xtrue(:)));

end

if arg.plot
	figure(1), im(arg.xtrue,'Original image','cbar',arg.clim);
	figure(2), im(y,'Blured noisy image','cbar',arg.clim);
	figure(3), im(xr,'Restored image','cbar',arg.clim);
	figure(4), im(xr-arg.xtrue,'Difference image','cbar');
end

if arg.chat
	printm(['\nCPU time: %6.2f\tActual time: %6.2f\nNRMSE (db): %6.2f\t' ...
		'NRMSE: %6.4f%%\nISNR: %6.4f\n'], time(end,1), time(end,2), err(end),...
		100*nrms(xr(omask), arg.xtrue(:)), isnr(end));
end

end


%% MFISTA
function [xk isnr err time itrtime] = restore_mfista(y, mask, dmask, omask, cfft, A, C, arg)

if strcmpi(arg.model_type, 'non-circ')
	arg.do_bfft = 0;
else
	arg.do_bfft = 1;
end

if strcmpi(arg.reg_type, 'l1,wave')
	arg.do_cfft = 1;
else
	arg.do_cfft = 0;
end

bfft = A.arg.psf_fft;

T = Gdiag2(ones(size(y)),'mask_in',mask);
d = A'*(T'*y);

xo = d;
uk = d;
zk = d;
vk = zeros(C.odim);
to = 1;
cost_old = inf;

isnr = zeros(arg.niter+1,1);
err = zeros(arg.niter+1,1);
time = zeros(arg.niter+1,2);
itrtime = zeros(arg.niter+1,2);

isnr(1) = 20*log10(norm(y(dmask)-arg.xtrue(:))/norm(xo(omask)-arg.xtrue(:)));
err(1) = 20*log10(nrms(xo(omask), arg.xtrue(:)));

tm1 = cputime;
tm2 = tic;

for ii=1:arg.niter

	if mod(ii,50) == 0 && arg.chat
		printm('iter: %d',ii)
	end

	% Inner iterations for updating zk
	if arg.do_bfft
		dk = uk - ifft2(abs(bfft).^2.*fft2(uk))/arg.L + d/arg.L;
	else
		dk = uk - A'*(T'*(T*(A*uk)))/arg.L + d/arg.L;
	end

	itm1 = cputime;
	itm2 = tic;
	for jj=1:arg.chiter
		Li = abs(C*zk)*arg.L/arg.lambda + arg.creg;
		if arg.do_cfft
			vk = real(arg.creg*vk - C*(ifft2(sum(conj(cfft).*fft2(vk),3)) - dk))./Li;
			zk = dk - real(ifft2(sum(conj(cfft).*fft2(vk),3)));
		else
			vk = (arg.creg*vk - C*(C'*vk - dk))./Li;
			zk = dk - C'*vk;
		end
	end
	itrtime(ii+1,1) = (cputime - itm1)/arg.chiter;
	itrtime(ii+1,2) = toc(itm2)/arg.chiter;

	tk = (1 + sqrt(1 + 4*to^2))/2;

	cost = sum(abs(vec(y-T*(A*zk))).^2)/2;
	if arg.is_iso
		cost = cost + arg.lambda*sum(vec(sqrt(sum(abs(C*zk).^2,3))));
	else
		cost = cost + arg.lambda*sum(abs(vec(C*zk)));
	end

	if cost < cost_old
		xk = zk;
	else
		xk = xo;
	end

	uk = xk + to/tk*(zk-xk) + (to-1)/tk*(xk-xo);

	to = tk;
	xo = xk;
	cost_old = cost;

	time(ii+1,1) = cputime - tm1;
	time(ii+1,2) = toc(tm2);
	err(ii+1) = 20*log10(nrms(xk(omask), arg.xtrue(:)));
	isnr(ii+1) = 20*log10(norm(y(dmask)-arg.xtrue(:))/norm(xk(omask)-arg.xtrue(:)));
end

if arg.plot
	figure(1), im(arg.xtrue,'Original image','cbar',arg.clim);
	figure(2), im(y,'Blured noisy image','cbar',arg.clim);
	figure(3), im(xk,'Restored image','cbar',arg.clim);
	figure(4), im(xk-arg.xtrue,'Difference image','cbar');
end

if arg.chat
	printm(['\nCPU time: %6.2f\tActual time: %6.2f\nNRMSE (db): %6.2f\t' ...
		'NRMSE: %6.4f%%\nISNR: %6.4f\n'], time(end,1), time(end,2), err(end),...
		100*nrms(xk(omask), arg.xtrue(:)), isnr(end));
end

end


%% DCT approach using MFISTA
function [xk isnr err time itrtime] = restore_mfista_dct(y, mask, dmask, omask, cfft, A, C, arg)

if strcmpi(arg.model_type, 'non-circ')
	arg.do_bfft = 0;
else
	arg.do_bfft = 0;
end

if strcmpi(arg.reg_type, 'l1,wave')
	arg.do_cfft = 1;
else
	arg.do_cfft = 0;
end

bfft = A.arg.psf_fft;

A = BlurMatrix(A.arg.psf,size(bfft),'refl');

T = Gdiag2(ones(size(y)),'mask_in',mask);
d = A'*(T'*y);

xo = d;
uk = d;
zk = d;
vk = zeros(C.odim);
to = 1;
cost_old = inf;

isnr = zeros(arg.niter+1,1);
err = zeros(arg.niter+1,1);
time = zeros(arg.niter+1,2);
itrtime = zeros(arg.niter+1,2);

isnr(1) = 20*log10(norm(y(dmask)-arg.xtrue(:))/norm(xo(omask)-arg.xtrue(:)));
err(1) = 20*log10(nrms(xo(omask), arg.xtrue(:)));

tm1 = cputime;
tm2 = tic;

for ii=1:arg.niter

	if mod(ii,50) == 0 && arg.chat
		printm('iter: %d',ii)
	end

	% Inner iterations for updating zk
	if arg.do_bfft
		dk = uk - ifft2(abs(bfft).^2.*fft2(uk))/arg.L + d/arg.L;
	else
		dk = uk - A'*(T'*(T*(A*uk)))/arg.L + d/arg.L;
	end

	itm1 = cputime;
	itm2 = tic;
	for jj=1:arg.chiter
		Li = abs(C*zk)*arg.L/arg.lambda + arg.creg;
		if arg.do_cfft
			vk = real(arg.creg*vk - C*(ifft2(sum(conj(cfft).*fft2(vk),3)) - dk))./Li;
			zk = dk - real(ifft2(sum(conj(cfft).*fft2(vk),3)));
		else
			vk = (arg.creg*vk - C*(C'*vk - dk))./Li;
			zk = dk - C'*vk;
		end
	end
	itrtime(ii+1,1) = (cputime - itm1)/arg.chiter;
	itrtime(ii+1,2) = toc(itm2)/arg.chiter;

	tk = (1 + sqrt(1 + 4*to^2))/2;

	cost = sum(abs(vec(y-T*(A*zk))).^2)/2;
	if arg.is_iso
		cost = cost + arg.lambda*sum(vec(sqrt(sum(abs(C*zk).^2,3))));
	else
		cost = cost + arg.lambda*sum(abs(vec(C*zk)));
	end

	if cost < cost_old
		xk = zk;
	else
		xk = xo;
	end

	uk = xk + to/tk*(zk-xk) + (to-1)/tk*(xk-xo);

	to = tk;
	xo = xk;
	cost_old = cost;

	time(ii+1,1) = cputime - tm1;
	time(ii+1,2) = toc(tm2);
	err(ii+1) = 20*log10(nrms(xk(omask), arg.xtrue(:)));
	isnr(ii+1) = 20*log10(norm(y(dmask)-arg.xtrue(:))/norm(xk(omask)-arg.xtrue(:)));
end

if arg.plot
	figure(1), im(arg.xtrue,'Original image','cbar',arg.clim);
	figure(2), im(y,'Blured noisy image','cbar',arg.clim);
	figure(3), im(xk,'Restored image','cbar',arg.clim);
	figure(4), im(xk-arg.xtrue,'Difference image','cbar');
end

if arg.chat
	printm(['\nCPU time: %6.2f\tActual time: %6.2f\nNRMSE (db): %6.2f\t' ...
		'NRMSE: %6.4f%%\nISNR: %6.4f\n'], time(end,1), time(end,2), err(end),...
		100*nrms(xk(omask), arg.xtrue(:)), isnr(end));
end

end


%% FISTA
function [xk isnr err time itrtime] = restore_fista(y, mask, cfft, A, C, arg)

if strcmpi(arg.model_type, 'non-circ')
	arg.do_bfft = 0;
else
	arg.do_bfft = 1;
end

if strcmpi(arg.reg_type, 'l1,wave')
	arg.do_cfft = 1;
else
	arg.do_cfft = 0;
end

bfft = A.arg.psf_fft;

T = Gdiag2(ones(size(y)),'mask_in',mask);
d = A'*(T'*y);

xo = d;
xk = d;
uk = d;
vk = zeros(C.odim);
to = 1;

isnr = zeros(arg.niter+1,1);
err = zeros(arg.niter+1,1);
time = zeros(arg.niter+1,2);
itrtime = zeros(arg.niter+1,2);

isnr(1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xk(mask)-arg.xtrue(:)));
err(1) = 20*log10(nrms(xk(mask), arg.xtrue(:)));

tm1 = cputime;
tm2 = tic;

for ii=1:arg.niter

	if mod(ii,50) == 0 && arg.chat
		printm('iter: %d',ii)
	end

	% Inner iterations for updating zk
	if arg.do_bfft
		dk = uk - ifft2(abs(bfft).^2.*fft2(uk))/arg.L + d/arg.L;
	else
		dk = uk - A'*(T'*(T*(A*uk)))/arg.L + d/arg.L;
	end

	itm1 = cputime;
	itm2 = tic;
	for jj=1:arg.chiter
		Li = abs(C*xk)*arg.L/arg.lambda + arg.creg;
		if arg.do_cfft
			vk = real(arg.creg*vk - C*(ifft2(sum(conj(cfft).*fft2(vk),3)) - dk))./Li;
			xk = dk - real(ifft2(sum(conj(cfft).*fft2(vk),3)));
		else
			vk = (arg.creg*vk - C*(C'*vk - dk))./Li;
			xk = dk - C'*vk;
		end
	end
	itrtime(ii+1,1) = (cputime - itm1)/arg.chiter;
	itrtime(ii+1,2) = toc(itm2)/arg.chiter;

	tk = (1 + sqrt(1 + 4*to^2))/2;
	uk = xk + (to-1)/tk*(xk-xo);

	to = tk;
	xo = xk;

	time(ii+1,1) = cputime - tm1;
	time(ii+1,2) = toc(tm2);
	err(ii+1) = 20*log10(nrms(xk(mask), arg.xtrue(:)));
	isnr(ii+1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xk(mask)-arg.xtrue(:)));
end

if arg.plot
	figure(1), im(arg.xtrue,'Original image','cbar',arg.clim);
	figure(2), im(y,'Blured noisy image','cbar',arg.clim);
	figure(3), im(xk,'Restored image','cbar',arg.clim);
	figure(4), im(xk-arg.xtrue,'Difference image','cbar');
end

if arg.chat
	printm(['\nCPU time: %6.2f\tActual time: %6.2f\nNRMSE (db): %6.2f\t' ...
		'NRMSE: %6.4f%%\nISNR: %6.4f\n'], time(end,1), time(end,2), err(end),...
		100*nrms(xk(mask), arg.xtrue(:)), isnr(end));
end


end

%% ISTA
function [xk isnr err time itrtime] = restore_ista(y, mask, cfft, A, C, arg)

if strcmpi(arg.model_type, 'non-circ')
	arg.do_bfft = 0;
else
	arg.do_bfft = 1;
end

if strcmpi(arg.reg_type, 'l1,wave')
	arg.do_cfft = 1;
else
	arg.do_cfft = 0;
end

bfft = A.arg.psf_fft;

T = Gdiag2(ones(size(y)),'mask_in',mask);
d = A'*(T'*y);

xk = d;
zk = zeros(C.odim);

isnr = zeros(arg.niter+1,1);
err = zeros(arg.niter+1,1);
time = zeros(arg.niter+1,2);
itrtime = zeros(arg.niter+1,2);

isnr(1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xk(mask)-arg.xtrue(:)));
err(1) = 20*log10(nrms(xk(mask), arg.xtrue(:)));

tm1 = cputime;
tm2 = tic;

for ii=1:arg.niter

	if mod(ii,50) == 0 && arg.chat
		printm('iter: %d',ii)
	end

	% Inner iterations for updating zk
	if arg.do_bfft
		dk = xk - real(ifft2(abs(bfft).^2.*fft2(xk)))/arg.L + d/arg.L;
	else
		dk = xk - A'*(T'*(T*(A*xk)))/arg.L + d/arg.L;
	end

	Li = abs(C*xk)*arg.L/arg.lambda + arg.creg;
	if arg.do_cfft
		zk = real(arg.creg*zk - C*(ifft2(sum(conj(cfft).*fft2(zk),3)) - dk))./Li;
		xk = dk - real(ifft2(sum(conj(cfft).*fft2(zk),3)));
	else
		zk = (arg.creg*zk - C*(C'*zk - dk))./Li;
		xk = dk - C'*zk;
	end

	time(ii+1,1) = cputime - tm1;
	time(ii+1,2) = toc(tm2);
	err(ii+1) = 20*log10(nrms(xk(mask), arg.xtrue(:)));
	isnr(ii+1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xk(mask)-arg.xtrue(:)));
end

if arg.plot
	figure(1), im(arg.xtrue,'Original image','cbar',arg.clim);
	figure(2), im(y,'Blured noisy image','cbar',arg.clim);
	figure(3), im(xk,'Restored image','cbar',arg.clim);
	figure(4), im(xk-arg.xtrue,'Difference image','cbar');
end

if arg.chat
	printm(['\nCPU time: %6.2f\tActual time: %6.2f\nNRMSE (db): %6.2f\t' ...
		'NRMSE: %6.4f%%\nISNR: %6.4f\n'], time(end,1), time(end,2), err(end),...
		100*nrms(xk(mask), arg.xtrue(:)), isnr(end));
end


end

%% NCG
function [xr isnr err time itrtime] = restore_ncg(y, mask, cfft, A, C, arg)

T = Gdiag2(ones(size(y)),'mask_in',mask);

d = A'*(T'*y);

if arg.precon
	bfft = real(A.arg.psf_fft);
	Hi = abs(bfft).^2 + arg.mu*sum(abs(cfft).^2,3);
	pfft = 1./Hi; % preconditioner
else
	pfft = 1;
end

xr = d;

Ax = T*(A*xr);
Cx = C*xr;

ddir = 0;
oldip = inf;

isnr = zeros(arg.niter+1,1);
err = zeros(arg.niter+1,1);
time = zeros(arg.niter+1,2);
itrtime = zeros(arg.niter+1,2);

isnr(1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xr(mask)-arg.xtrue(:)));
err(1) = 20*log10(nrms(xr(mask), arg.xtrue(:)));

tm1 = cputime;
tm2 = tic;

for ii=1:arg.niter

	if mod(ii,50) == 0 && arg.chat
		printm('iter: %d',ii)
	end

	wi = 1./sqrt(abs(Cx).^2 + arg.epsilon);
	Wi = Gdiag(wi);

	ngrad = A'*(T'*(y-Ax)) - arg.lambda*(C'*(Wi*Cx));

	if arg.precon
		pregrad = real(ifft2(pfft.*fft2(ngrad)));
	else
		pregrad = ngrad;
	end

	newip = ngrad(:)' * pregrad(:);

	gamma = newip / oldip;
	ddir = pregrad + gamma * ddir;

	oldip = newip;

	Adir = T*(A*ddir);
	Cdir = C*ddir;

	% Do line search
	dAWAd = real(Adir(:)'*Adir(:));
	dAWr = real(Adir(:)'*vec(y-Ax));
	step = 0;

	itm1 = cputime;
	itm2 = tic;
	for is=1:arg.lniter
		wi = 1./sqrt(abs(Cx+step*Cdir).^2 + arg.epsilon);
		Wi = Gdiag(wi);

		denom = dAWAd + arg.lambda*real(Cdir(:)'*(Wi*Cdir(:)));
		if denom == 0 || isinf(denom)
			warning('restore_ncg:denom','0 or inf denom?');
		end
		pdot = arg.lambda*real(Cdir(:)'*(Wi*vec(Cx + step*Cdir)));

		step = step - (-dAWr + step * dAWAd + pdot) / denom;
	end
	itrtime(ii+1,1) = (cputime - itm1)/arg.lniter;
	itrtime(ii+1,2) = toc(itm2)/arg.lniter;

	if step < 0
		warning('inner_cg:step','not descent direction, downhill?');
	end
	% End line search

	Ax = Ax + step * Adir;
	Cx = Cx + step * Cdir;
	xr = xr + step * ddir;

	time(ii+1,1) = cputime - tm1;
	time(ii+1,2) = toc(tm2);
	err(ii+1) = 20*log10(nrms(xr(mask), arg.xtrue(:)));
	isnr(ii+1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xr(mask)-arg.xtrue(:)));
end

if arg.plot
	figure(1), im(arg.xtrue,'Original image','cbar',arg.clim);
	figure(2), im(y,'Blured noisy image','cbar',arg.clim);
	figure(3), im(xr,'Restored image','cbar',arg.clim);
	figure(4), im(xr-arg.xtrue,'Difference image','cbar');
end

if arg.chat
	printm(['\nCPU time: %6.2f\tActual time: %6.2f\nNRMSE (db): %6.2f\t' ...
		'NRMSE: %6.4f%%\nISNR: %6.4f\n'], time(end,1), time(end,2), err(end),...
		100*nrms(xr(mask), arg.xtrue(:)), isnr(end));
end

end

%% Quadratic CG
function [xr isnr err time itrtime] = restore_cg_quad(y, mask, cfft, A, C, arg)

if strcmpi(arg.reg_type, 'l1,wave')
	arg.do_cfft = 1;
else
	arg.do_cfft = 0;
end

if strcmpi(arg.reg_type, 'l1,wave') && ~isscalar(C)
	arg.use_fft = 1;
else
	arg.use_fft = 0;
end

arg.converged = 0;
cfft = sum(abs(cfft).^2,3);

T = Gdiag2(ones(size(y)),'mask_in',mask);
bfft = A.arg.psf_fft;

yd = A'*(T'*y);
Hi = abs(bfft).^2 + arg.lambda*cfft;
pfft = 1./Hi;

printm('Cond number of Hmu: %f', max(Hi(:))/min(Hi(:)))

xr = yd;

Gx = A'*(T'*(T*(A*xr)));
if arg.use_fft
	Rx = arg.lambda*real(ifft2(cfft.*fft2(xr)));
else
	Rx = arg.lambda*(C'*(C*xr));
end

ddir = 0;
oldip = inf;

isnr = zeros(arg.niter+1,1);
err = zeros(arg.niter+1,1);
time = zeros(arg.niter+1,2);
itrtime = zeros(arg.niter+1,2);

isnr(1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xr(mask)-arg.xtrue(:)));
err(1) = 20*log10(nrms(xr(mask), arg.xtrue(:)));

tm1 = cputime;
tm2 = tic;

for ii=1:arg.niter

	if mod(ii,50) == 0 && arg.chat
		printm('iter: %d',ii)
	end

	ngrad = yd - Gx - Rx;

	if arg.precon
		pregrad = real(ifft2(pfft.*fft2(ngrad)));
	else
		pregrad = ngrad;
	end

	newip = ngrad(:)' * pregrad(:);

	gamma = newip / oldip;
	ddir = pregrad + gamma * ddir;

	oldip = newip;

	Gdir = A'*(T'*(T*(A*ddir)));
	if arg.use_fft
		Rdir = arg.lambda*real(ifft2(cfft.*fft2(ddir)));
	else
		Rdir = arg.lambda*(C'*(C*ddir));
	end

	denom = ddir(:)' * vec(Gdir + Rdir);

	if denom == 0
		warning('inner_cg:denom','found exact solution??? step=0 now!?');
		step = 0;
	else
		step = (ddir(:)' * ngrad(:)) / denom;
		step = real(step); % real step seems only natural
	end

	if step < 0
		warning('inner_cg:step','not descent direction, downhill?');
	end

	Gx = Gx + step * Gdir;
	Rx = Rx + step * Rdir;
	xr = xr + step * ddir;

	if norm(step*ddir(:)) < norm(xr(:))*arg.cg_tol && ~arg.converged
		printm('Convergence reached at iteration %d', ii);
		arg.converged = 1;
	end

	time(ii+1,1) = cputime - tm1;
	time(ii+1,2) = toc(tm2);
	err(ii+1) = 20*log10(nrms(xr(mask), arg.xtrue));
	isnr(ii+1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xr(mask)-arg.xtrue(:)));

end

if arg.plot
	figure(1), im(arg.xtrue,'Original image','cbar',arg.clim);
	figure(2), im(y,'Blured noisy image','cbar',arg.clim);
	figure(3), im(xr,'Restored image','cbar',arg.clim);
	figure(4), im(xr-arg.xtrue,'Difference image','cbar');
end

if arg.chat
	printm(['\nCPU time: %6.2f\tActual time: %6.2f\nNRMSE (db): %6.2f\t' ...
		'NRMSE: %6.4f%%\nISNR: %6.4f\n'], time(end,1), time(end,2), err(end),...
		100*nrms(xr(mask), arg.xtrue(:)), isnr(end));
end

end

%% Scalar and Vector Shrinkage function
function u = shrink(d, lam, type)

if nargin < 3
	type = 'scalar';
end

if strcmpi(type, 'scalar')
	u = d./abs(d).*max(abs(d) - lam,0);
elseif strcmpi(type, 'vector')
	nz = size(d,3);
	dn = repmat(sqrt(sum(abs(d).^2,3)),[1 1 nz]);
	u = d./dn.*max(dn-lam,0);
else
	error('Unknown shrink type');
end

if any(isnan(u(:)))
	u(isnan(u)) = 0;
end

end

%% AL algorithm with quadratic regularization
function [xr isnr err time itrtime] = restore_al_quad(y, mask, cfft, A, C, arg)

if strcmpi(arg.reg_type, 'l1,wave')
	arg.do_cfft = 1;
else
	arg.do_cfft = 0;
end

bfft = A.arg.psf_fft;

Hi = abs(bfft).^2 + arg.lambda/arg.mu*sum(abs(cfft).^2,3);
Hd = 1./(double(mask)+arg.mu);

printm('Cond number of Hnu: %f', max(Hi(:))/min(Hi(:)))
printm('Cond number of Hmu: %f', (1+arg.mu)/arg.mu);

% Next 2 lines possibly redundant. Keep for safety
T = Gdiag2(ones(size(y)),'mask_in',mask);
d = T'*y;

xr = A'*(T'*y);
Axr = A*xr;

eta0 = zeros(C.idim);

isnr = zeros(arg.niter+1,1);
err = zeros(arg.niter+1,1);
time = zeros(arg.niter+1,2);
itrtime = zeros(arg.niter+1,2);

isnr(1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xr(mask)-arg.xtrue(:)));
err(1) = 20*log10(nrms(xr(mask), arg.xtrue(:)));

tm1 = cputime;
tm2 = tic;

for ii=1:arg.niter

	if mod(ii,arg.skip) == 0 && arg.chat
		printm('iter: %d',ii)
	end

	u0 = Hd.*(d + arg.mu*(Axr + eta0));
	if arg.extra
		eta0 = eta0 - u0 + Axr;
	end

	xr = real(ifft2((conj(bfft).*fft2(u0-eta0))./Hi));
	Axr = A*xr;

	eta0 = eta0 - u0 + Axr;

	time(ii+1,1) = cputime - tm1;
	time(ii+1,2) = toc(tm2);
	err(ii+1) = 20*log10(nrms(xr(mask), arg.xtrue(:)));
	isnr(ii+1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xr(mask)-arg.xtrue(:)));
end

if arg.plot
	figure(1), im(arg.xtrue,'Original image','cbar',arg.clim);
	figure(2), im(y,'Blured noisy image','cbar',arg.clim);
	figure(3), im(xr,'Restored image','cbar',arg.clim);
	figure(4), im(xr-arg.xtrue,'Difference image','cbar');
end

if arg.chat
	printm(['\nCPU time: %6.2f\tActual time: %6.2f\nNRMSE (db): %6.2f\t' ...
		'NRMSE: %6.4f%%\nISNR: %6.4f\n'], time(end,1), time(end,2), err(end),...
		100*nrms(xr(mask), arg.xtrue(:)), isnr(end));
end

end

%% PCG algorithm for inner iterations
function [xr time] = inner_cg(y, xi, A, T, C, pfft, cfft, arg)

if strcmpi(arg.reg_type, 'l1,wave') && ~isscalar(C)
	arg.use_fft = 1;
else
	arg.use_fft = 0;
end

cfft = sum(abs(cfft).^2,3);

xr = xi;

Gx = A'*(T'*(T*(A*xr)));
if arg.use_fft
	Rx = arg.mu*real(ifft2(cfft.*fft2(xr)));
else
	Rx = arg.mu*(C'*(C*xr));
end

ddir = 0;
oldip = inf;

tm1 = cputime;
tm2 = tic;

for iter = 1:arg.cgiter

	ngrad = y - Gx - Rx;

	if arg.precon
		pregrad = real(ifft2(pfft.*fft2(ngrad)));
	else
		pregrad = ngrad;
	end

	newip = ngrad(:)' * pregrad(:);

	gamma = newip / oldip;
	ddir = pregrad + gamma * ddir;

	oldip = newip;

	Gdir = A'*(T'*(T*(A*ddir)));
	if arg.use_fft
		Rdir = arg.mu*real(ifft2(cfft.*fft2(ddir)));
	else
		Rdir = arg.mu*(C'*(C*ddir));
	end

	denom = ddir(:)' * vec(Gdir + Rdir);

	if denom == 0
		warning('inner_cg:denom','found exact solution??? step=0 now!?');
		step = 0;
	else
		step = (ddir(:)' * ngrad(:)) / denom;
		step = real(step); % real step seems only natural
	end

	if step < 0
		warning('inner_cg:step','not descent direction, downhill?');
	end

	Gx = Gx + step * Gdir;
	Rx = Rx + step * Rdir;
	xr = xr + step * ddir;

	if norm(step*ddir(:)) < norm(xr(:))*arg.cg_tol
%		printm('Convergence reached at iteration %d', iter);
		break
	end
end

time = [cputime-tm1 toc(tm2)]/iter;

end

%% PCG algorithm for residual boundaries in 'Reeves' method
function [xr time] = inner_cg_reeves(y, xi, Ac, hfft, imsiz, arg)

xr = xi;

Gx = xr - Ac*vec(ifft2(hfft.*fft2(reshape(Ac'*xr(:), imsiz))));

ddir = 0;
oldip = inf;

tm1 = cputime;
tm2 = tic;

for iter = 1:arg.cgiter

	ngrad = y - Gx;

	pregrad = ngrad;
	newip = ngrad(:)' * pregrad(:);

	gamma = newip / oldip;
	ddir = pregrad + gamma * ddir;

	oldip = newip;

	Gdir = ddir - Ac*vec(ifft2(hfft.*fft2(reshape(Ac'*ddir(:), imsiz))));

	denom = ddir(:)' * Gdir;

	if denom == 0
		warning('inner_cg:denom','found exact solution??? step=0 now!?');
		step = 0;
	else
		step = (ddir(:)' * ngrad(:)) / denom;
		step = real(step); % real step seems only natural
	end

	if step < 0
		warning('inner_cg:step','not descent direction, downhill?');
	end

	Gx = Gx + step * Gdir;
	xr = xr + step * ddir;

	if norm(step*ddir(:)) < norm(xr(:))*arg.cg_tol
%		printm('Convergence reached at iteration %d', iter);
		break
	end
end

time = [cputime-tm1 toc(tm2)]/iter;

end


%% Testing function (empty)
function image_restore_test

end
