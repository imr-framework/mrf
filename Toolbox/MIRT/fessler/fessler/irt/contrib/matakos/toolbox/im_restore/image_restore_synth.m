 function [xr, out] = image_restore_synth(y, psf, varargin)
%function [wk, out] = image_restore_synth(y, psf, varargin)
%|
%| Function to perform image restoration using Augmented Lagrangian methods and
%| FISTA for synthesis formulation. Regularization is L1 norm and transform is
%| orthonormal circulant wavelet.
%|
%| in
%|	y	[nx ny]		Input noisy and/or blurred image
%|	psf	[nbx nby]	PSF of blur
%|
%| options
%|	est_type	string	Type of restoration. Valid choices
%|		'ncg'				Non-linear CG algorithm with rounding parameter
%|							and line search method
%|		'ista'				ISTA method
%|		'fista'				FISTA method
%|		'mfista'			MFISTA method
%|		'salsa'				SALSA algorithm
%|		'alp2'				AL method with two splits, one for regularization
%|		'dualal'			Dual AL method
%|
%|	model_type	string	Type of reconstruction model
%|		'circ'				Unrealistic circulant model
%|		'non-circ'			Accurate non-circulant model
%|
%|	psf_type	string		type of Gblur to be used
%|		'fft,same'			Circulant blur using fftn
%|		'conv,same'			Toeplitz blur using convn
%|
%|	wave_type	string		Type of wavelet for wavelet l1 norm regularization
%|							Valid parameters from wfilters
%|	lambda		scalar		regularization parameter
%|	mu			scalar		Lagrangian tuning parameter
%|	nu			scalar		Lagrangian tuning parameter
%|	xtrue		[nx ny]		True image for reference
%|	niter		scalar		# of iterations
%|	aliter		scalar		# of inner iteration for dual AL
%|	cgiter		scalar		# of CG iterations for Newton method in dual AL
%|	lniter		scaalr		# of inner line search iterations
%|	plot		0|1			plot results
%|	chat		0|1			show printouts
%|
%| Copyright 2011-10-31, Antonis Matakos, University of Michigan

if nargin < 2
	image_restore_test;
	return;
end

arg.est_type = 'alp2';
arg.model_type = 'circ';
arg.psf_type = 'fft,same';
arg.wave_type = 'haar';
arg.lambda = 2^-12;
arg.mu = 2^-1;
arg.nu = 1;
arg.xtrue = 1;
arg.niter = 500;
arg.aliter = 5;
arg.cgiter = 4;
arg.lniter = 5;
arg.altol = 1e-21; % Dual AL tolerance parametner - hidden parameter
arg.cg_tol = 1e-6; % Tolerance of inner CG - hidden parameter (OK)
arg.epsilon = 1e-10; % Rounding parameter for NCG - hidden parameter
arg.precon = 1; % Preconditioner use for inner CG - hidden parameter (YES)
arg.L = 1; % Lipschitz constant - hidden parameter (Leave unchanged)
arg.extra = 0; % Extra update of eta params - hidden parameter (NO)
arg.reg_approx = 0; % Apply thresholding on approximation level - hidden parameter (NO)
arg.skip = 50;
arg.plot = 0;
arg.chat = 0;

arg = vararg_pair(arg, varargin);

[nx ny] = size(y);
[nbx nby] = size(psf);

if strcmpi(arg.model_type, 'non-circ')
	imask = true(nx+nbx-1,nx+nby-1);
	mask = false(nx+nbx-1,nx+nby-1);
	mask((nbx+1)/2:end-(nbx-1)/2,(nby+1)/2:end-(nby-1)/2) = true;
else
	imask = true(nx,ny);
	mask = true(nx,ny);
end

A = Gblur(imask, 'psf', psf, 'type', arg.psf_type);

C = Gwave2('mask', imask, 'nlevel', 2, 'wname', arg.wave_type, 'includeApprox', 1);

arg.clim = [0 max(arg.xtrue(:))+0.1];

if strcmpi(arg.est_type, 'alp2')
	[xr isnr err time] = restore_alp2(y, mask, A, C, arg);
elseif strcmpi(arg.est_type, 'salsa')
	[xr isnr err time itrtime] = restore_salsa(y, mask, A, C, arg);
elseif strcmpi(arg.est_type, 'mfista')
	[xr isnr err time] = restore_mfista(y, mask, A, C, arg);
elseif strcmpi(arg.est_type, 'fista')
	[xr isnr err time] = restore_fista(y, mask, A, C, arg);
elseif strcmpi(arg.est_type, 'ista')
	[xr isnr err time] = restore_ista(y, mask, A, C, arg);
elseif strcmpi(arg.est_type, 'ncg')
	[xr isnr err time itrtime] = restore_ncg(y, mask, A, C, arg);
else
	fail('Unknown estimation method %s', arg.est_type)
end

% xr = reale(C'*wk,'warn','non-negligible imaginary part');

out.isnr = isnr;
out.err = err;
out.time = time;
out.psf = psf;

out.cputime_iter = time(end,1)/arg.niter;
out.time_iter = time(end,2)/arg.niter;

if strcmpi(arg.est_type, 'salsa')
	out.cputime_cgiter = mean(itrtime(:,1));
	out.time_cgiter = mean(itrtime(:,2));
	out.cputime_lniter = 0;
	out.time_lniter = 0;
elseif strcmpi(arg.est_type, 'ncg')
	out.cputime_cgiter = 0;
	out.time_cgiter = 0;
	out.cputime_lniter = mean(itrtime(:,1));
	out.time_lniter = mean(itrtime(:,2));
else
	out.cputime_cgiter = 0;
	out.time_cgiter = 0;
	out.cputime_lniter = 0;
	out.time_lniter = 0;
end

end % image_restore_synth()



function [xr isnr err time] = restore_alp2(y, mask, A, C, arg)

bfft = A.arg.psf_fft;

if arg.reg_approx
	dg = 1;
else
	dg = ones(C.odim);
	dg(:,:,end) = 0;
end

% Ha = AC'CA' = A'A
Ha = abs(bfft).^2;
% Hb = A'(I - (AC'CA' + nI)^-1*AC'CA')
Hb = conj(bfft).*(1 - Ha./(Ha + arg.nu));
% Hc = A'(AC'CA' + nI)^-1*A
Hc = abs(bfft).^2./(Ha + arg.nu);
% Hd = (T'T + mI)^-1
Hd = 1./(double(mask)+arg.mu);

printm('Cond number of Hmu: %f', 1+1/arg.mu)
printm('Cond number of Hnu: %f', 1+1/arg.nu);

T = Gdiag2(ones(size(y)), 'mask_in', mask);
d = T'*y;

xr = A'*d;
wk = C*xr;
ACtwr = A*(C'*wk);

eta0 = zeros(C.idim);
eta1 = zeros(C.odim);

isnr = zeros(arg.niter+1,1);
err = zeros(arg.niter+1,1);
time = zeros(arg.niter+1,2);

isnr(1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xr(mask)-arg.xtrue(:)));
err(1) = 20*log10(nrms(xr(mask), arg.xtrue(:)));

tm1 = cputime;
tm2 = tic;

for ii=1:arg.niter

	if mod(ii,arg.skip) == 0 && arg.chat
		printm('iter: %d',ii)
	end

	% u0 = (T'T + mI)^-1(T'y + m(AC'w + eta0))
	u0 = Hd.*(d + arg.mu*(ACtwr + eta0));
%	eta0 = eta0 - u0 + ACtwr;

	u1 = shrink(wk + eta1, dg*arg.lambda/arg.mu/arg.nu);
%	eta1 = eta1 - u1 + wk;

	d0 = u0 - eta0;
	d1 = u1 - eta1;
	% wk = d1 + CA'((I - (AC'CA' + nI)^-1*AC'CA')d0 - (AC'CA' + nI)^-1*AC'd1)
	wk = d1 + C*(ifft2(Hb.*fft2(d0)/arg.nu - Hc.*fft2(C'*d1)));

	xr = C'*wk;
	wk = C*xr;

%	ACtwr = A*(C'*wk);
	ACtwr = A*xr;

	eta0 = eta0 - u0 + ACtwr;
	eta1 = eta1 - u1 + wk;

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

end % restore_alp2()



function [xr isnr err time itrtime] = restore_salsa(y, mask, A, C, arg)

if strcmpi(arg.model_type, 'non-circ')
	arg.do_cg = 1;
else
	arg.do_cg = 0;
end

if arg.reg_approx
	dg = 1;
else
	dg = ones(C.odim);
	dg(:,:,end) = 0;
end

bfft = A.arg.psf_fft;
T = Gdiag2(ones(size(y)),'mask_in',mask);

% Create circulant preconditioner
ej = zeros(size(y));
ej(1,1) = 1;
ej = fftshift(ej);
psf = fftshift(T*(A*(A'*(T'*ej))));
psf_fft = real(fft2(psf));
pfft = 1./(psf_fft + arg.mu);

% psf = conv2(A.arg.psf,A.arg.psf);
% bccb = create_bccb(psf,size(y),'psf');
% psf_fft = real(fft2(bccb));
% pfft = 1./(psf_fft + arg.mu);

% Ha = AC'CA' and C'C = I
Ha = abs(bfft).^2;
% Hb = A'(I - (AC'CA' + mI)^-1*AC'CA')
Hb = conj(bfft).*(1 - Ha./(Ha + arg.mu));
% Hc = A'(AC'CA' + mI)^-1*A
Hc = abs(bfft).^2./(Ha + arg.mu);

% d = CA'(I - (AC'CA' + mI)^-1*AC'CA')y/mu
d = C*ifft2(Hb.*fft2(T'*y))/arg.mu;
dd = C*(A'*(T'*y));

printm('Cond number of Hmu: %f', 1+1/arg.mu);

xr = A'*(T'*y);
wk = C*xr;
vk = T*(A*(C'*dd));
eta = zeros(C.odim);

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

	uk = shrink(wk + eta, dg*arg.lambda/arg.mu);
%	eta = eta - uk + wk;

	% 1/2||y-TAC'w|| + lam||u|| + mu/2||u-w-eta||
	% wkp1 = (CAT'TAC' + muI)^-1*(CA'T'y + mu(u-eta))

%	wk = d + uk - eta - CA'(AC'CA' + mI)^-1*AC'(uk-eta);
	if arg.do_cg
		yd = (dd + arg.mu*(uk-eta))/arg.mu;
		[vk tm] = inner_cg(T*(A*(C'*yd)), vk, A, T, pfft, arg);
		wk = yd - C*(A'*(T'*vk));
		itrtime(ii+1,1) = tm(1);
		itrtime(ii+1,2) = tm(2);
	else
		wk = d + uk - eta - C*ifft2(Hc.*fft2(C'*(uk-eta)));
	end

	xr = C'*wk;
	wk = C*xr;

	eta = eta - uk + wk;

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



function [xr isnr err time] = restore_mfista(y, mask, A, C, arg)

if arg.reg_approx
	dg = 1;
else
	dg = ones(C.odim);
	dg(:,:,end) = 0;
end

T = Gdiag2(ones(size(y)),'mask_in',mask);
Th = Gdiag(double(mask));

xo = A'*(T'*y);
d = C*xo;

wo = d;
uk = d;
to = 1;
cost_old = inf;

isnr = zeros(arg.niter+1,1);
err = zeros(arg.niter+1,1);
time = zeros(arg.niter+1,2);

isnr(1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xo(mask)-arg.xtrue(:)));
err(1) = 20*log10(nrms(xo(mask), arg.xtrue(:)));

tm1 = cputime;
tm2 = tic;

for ii=1:arg.niter

	if mod(ii,arg.skip) == 0 && arg.chat
		printm('iter: %d',ii)
	end

	dk = uk - C*(A'*(Th*(A*(C'*uk))))/arg.L + d/arg.L;
	zk = shrink(dk, dg*arg.lambda/arg.L);

	xr = C'*zk;
	zk = C*xr;

	tk = (1 + sqrt(1 + 4*to^2))/2;

	cost = sum(abs(vec(y-T*(A*xr))).^2)/2 + arg.lambda*sum(abs(zk(:)));

	if cost < cost_old
		wk = zk;
	else
		wk = wo;
		xr = xo;
	end

	uk = wk + to/tk*(zk-wk) + (to-1)/tk*(wk-wo);

	to = tk;
	wo = wk;
	xo = xr;
	cost_old = cost;

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



function [xr isnr err time] = restore_fista(y, mask, A, C, arg)

if arg.reg_approx
	dg = 1;
else
	dg = ones(C.odim);
	dg(:,:,end) = 0;
end

T = Gdiag2(ones(size(y)),'mask_in',mask);
Th = Gdiag(double(mask));

xr = A'*(T'*y);
d = C*xr;


wo = d;
uk = d;
to = 1;

isnr = zeros(arg.niter+1,1);
err = zeros(arg.niter+1,1);
time = zeros(arg.niter+1,2);

isnr(1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xr(mask)-arg.xtrue(:)));
err(1) = 20*log10(nrms(xr(mask), arg.xtrue(:)));

tm1 = cputime;
tm2 = tic;

for ii=1:arg.niter

	if mod(ii,arg.skip) == 0 && arg.chat
		printm('iter: %d',ii)
	end

	dk = uk - C*(A'*(Th*(A*(C'*uk))))/arg.L + d/arg.L;
	wk = shrink(dk, dg*arg.lambda/arg.L);

	xr = C'*wk;
	wk = C*xr;

	tk = (1 + sqrt(1 + 4*to^2))/2;
	uk = wk + (to-1)/tk*(wk-wo);

	to = tk;
	wo = wk;

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



function [xr isnr err time] = restore_ista(y, mask, A, C, arg)

if arg.reg_approx
	dg = 1;
else
	dg = ones(C.odim);
	dg(:,:,end) = 0;
end

T = Gdiag2(ones(size(y)),'mask_in',mask);
Th = Gdiag(double(mask));

xr = A'*(T'*y);
d = C*xr;

wk = d;

isnr = zeros(arg.niter+1,1);
err = zeros(arg.niter+1,1);
time = zeros(arg.niter+1,2);

isnr(1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xr(mask)-arg.xtrue(:)));
err(1) = 20*log10(nrms(xr(mask), arg.xtrue(:)));

tm1 = cputime;
tm2 = tic;

for ii=1:arg.niter

	if mod(ii,arg.skip) == 0 && arg.chat
		printm('iter: %d',ii)
	end

	dk = wk - C*(A'*(Th*(A*(C'*wk))))/arg.L + d/arg.L;
	wk = shrink(dk, dg*arg.lambda/arg.L);

	xr = C'*wk;
	wk = C*xr;

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



function [xr isnr err time itrtime] = restore_ncg(y, mask, A, C, arg)

if arg.reg_approx
	dg = 1;
else
	dg = ones(C.odim);
	dg(:,:,end) = 0;
end

T = Gdiag2(ones(size(y)),'mask_in',mask);


d = C*(A'*(T'*y));

% if arg.precon
%	bfft = real(A.arg.psf_fft);
%	Hi = abs(bfft).^2 + arg.mu*sum(abs(cfft).^2,3);
%	pfft = 1./Hi; % preconditioner
% else
%	pfft = 1;
% end

wk = d;
xr = C'*wk;

Gw = T*(A*xr);

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

	wi = 1./sqrt(abs(wk).^2 + arg.epsilon);
	Wi = Gdiag(wi);

	ngrad = C*(A'*(T'*(y-Gw))) - dg.*(arg.lambda*(Wi*wk));

	if arg.precon
%		pregrad = real(ifft2(pfft.*fft2(ngrad)));
		pregrad = ngrad;
	else
		pregrad = ngrad;
	end

	newip = ngrad(:)' * pregrad(:);

	gamma = newip / oldip;	% Fletcher-Reeves
	ddir = pregrad + gamma * ddir;

	oldip = newip;

	Gdir = T*(A*(C'*ddir));

	% Do line search
	dAWAd = real(Gdir(:)'*Gdir(:));
	dAWr = real(Gdir(:)'*vec(y-Gw));
	dd = dg.*ddir;
	step = 0;

	itm1 = cputime;
	itm2 = tic;
	for is=1:arg.lniter
		wi = 1./sqrt(abs(wk+step*ddir).^2 + arg.epsilon);
		Wi = Gdiag(wi);

		denom = dAWAd + arg.lambda*real(dd(:)'*(Wi*dd(:)));
		if denom == 0 || isinf(denom)
			warning('restore_ncg:denom','0 or inf denom?');
		end
		pdot = arg.lambda*real(dd(:)'*(Wi*vec(wk + step*dd)));

		step = step - (-dAWr + step * dAWAd + pdot) / denom;
	end
	itrtime(ii+1,1) = (cputime - itm1)/arg.lniter;
	itrtime(ii+1,2) = toc(itm2)/arg.lniter;

	if step < 0
		warning('inner_cg:step','not descent direction, downhill?');
	end
	% End line search
	xdir = C'*ddir;
	ddir = C*xdir;

	xr = xr + step * xdir;
	wk = wk + step * ddir;
	Gw = T*(A*xr);

%	xr = C'*wk;
%	wk = C*xr;
%	Gw = Gw + step * Gdir;
%	wk = wk + step * ddir;

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



function [wk time] = inner_cg(y, wi, A, T, pfft, arg)

bfft = A.arg.psf_fft;
bfft = abs(bfft).^2;

wk = wi;
% Gx = (T*(A*(A'*(T'*wk))));
Gx = real(T*ifft2(bfft.*fft2(T'*wk)));
Rx = arg.mu*wk;

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

	gamma = newip / oldip;	% Fletcher-Reeves
	ddir = pregrad + gamma * ddir;

	oldip = newip;

	Gdir = real(T*ifft2(bfft.*fft2(T'*ddir)));
%	Gdir = C*(A'*(T'*(T*(A*(C'*ddir)))));
	Rdir = arg.mu*ddir;

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
	wk = wk + step * ddir;

	if norm(step*ddir(:)) < norm(wk(:))*arg.cg_tol
%		printm('Convergence reached at iteration %d', iter);
		break
	end
end

time = [cputime-tm1 toc(tm2)]/iter;

end



function u = shrink(d, lam)

u = sign(d).*max(abs(d) - lam,0);

end



% function [wk err time] = restore_dualal(y, mask, cfft, A, C, arg)
%
% bfft = A.arg.psf_fft;
%
% T = Gdiag2(ones(size(y)),'mask_in',mask);
% % In = Gdiag(double(mask));
% d = C*(A'*(T'*y));
%
% ak = T*(A*(C'*d));
% wk = d;
% % pk = ak;
% muk = arg.mu;
%
% flag_converged = 0;
% err = zeros(arg.niter+1,1);
% time = zeros(arg.niter+1,1);
%
% err(1) = 20*log10(nrms(wk(:), arg.wtrue(:)));
% time(1) = 0;
%
% As = sparse_bccb(A);
% Cs = sparse_bccb(C);
%
% % Tf = double(full(T));
% % mkt = abs(Tf) > 1e-7;
% % Tf(~mkt) = 0;
% %
% % Ts = sparse(Tf);
% Ts = speye(prod(C.idim));
%
% B = Ts*(As*Cs');
%
% Iv = speye(prod(C.idim));
%
% H = muk*(B*B') + Iv;
%
% tm = cputime;
%
% % Ct = full(C');
% % At = T*A*Ct;
%
% for ii=1:arg.niter
%
%	if mod(ii,arg.skip) == 0 && arg.chat
%		printm('iter: %d',ii)
%	end
%
%	if ~flag_converged
%	qk = shrink(wk + muk*(C*(A'*(T'*ak))), arg.lambda*muk);
%	gk = ak - y + T*(A*(C'*qk));
%
%	nrg = sum(abs(vec(gk)).^2);
%	nrw = sum(abs(vec(qk-wk)).^2)/muk;
%
%	jj = 0;
%	while nrg > nrw && jj <= arg.aliter
%
%		mk = abs(qk) > 0;
%
% %		Th = Gdiag(double(mk(:)));
% %		Th = sparse(Th);
% %		Hk = eye(C.idim) + muk*T*A*C'*(Tt'*Tt)*C*A'*T';
% %		Hk = Iv + muk*At*(Th*At');
%		tm2 = cputime;
%		Bs = B(:,mk(:));
%		Hk = Iv + muk*(Bs*Bs');
%		disp(cputime-tm2)
% %		P = qpwls_precon('circ0', {In}, sqrt(muk)*C*A'*T', mask);
%
%		% Use CG to find pk as Hk * pk = gk;
% %		pk = zeros(C.idim);
% %		pk = qpwls_pcg2(pk(:), In, gk(:), sqrt(muk)*Tt*C*A'*T', ...
% %			'precon', P,'stop_threshold', 1e-3, 'niter', arg.cgiter);
%		tm2 = cputime;
%		pk = Hk\gk(:);
%		disp(cputime-tm2)
%		pk = reshape(pk, C.idim);
%
%
%		% find optimal step size
% %		gamma = 1;
%		if (nrg-nrw) < 1e-1
%			gamma = 1;
%		else
%			gamma = 1;
%			cost_old = sum(abs(vec(ak - y)).^2) + sum(abs(vec(qk)).^2);
%
%			cost = sum(abs(vec(ak - gamma*pk - y)).^2);
%			vk = shrink(wk + muk*(C*(A'*(T'*(ak - gamma*pk)))), arg.lambda*muk);
%			cost = cost + sum(abs(vec(vk)).^2);
%
%			while cost > cost_old
%				gamma = gamma/2;
%
%				cost = sum(abs(vec(ak - gamma*pk - y)).^2);
%				vk = shrink(wk + muk*(C*(A'*(T'*(ak - gamma*pk)))), arg.lambda*muk);
%				cost = cost + sum(abs(vec(vk)).^2);
%			end
%		end
%
%		ak = ak - gamma*pk;
%
%		qk = shrink(wk + muk*(C*(A'*(T'*ak))), arg.lambda*muk);
%		gk = ak - y + T*(A*(C'*qk));
%
%		nrg = sum(abs(vec(gk)).^2);
%		nrw = sum(abs(vec(qk-wk)).^2)/muk;
%
%		jj = jj+1;
%		if mod(jj,1) == 0 && arg.chat
%			printm('inner iter: %d',jj)
%			printm('diff: %e gamma: %2.4f', nrg-nrw, gamma)
%		end
%
%		if nrg-nrw > 0 && abs(nrg-nrw) < arg.altol
%			printm('Reached convergence at iter: %d, inner iter: %d\n', ii, jj);
%			flag_converged = 1;
%			break;
%		end
%	end
%
%	wk = qk;
%	muk = 2*muk;
%	end
%
%	time(ii+1) = cputime - tm;
%	err(ii+1) = 20*log10(nrms(wk(:), arg.wtrue(:)));
% end
%
% disp(cputime-tm);
%
% if arg.plot
%	figure(1), im(arg.wtrue,'Original coefficients','cbar',arg.clim);
%	figure(2), im(wk,'Restored coefficients','cbar',arg.clim);
%	figure(3), im(wk-arg.wtrue,'Difference coefficients','cbar');
% end
%
% if arg.chat && arg.plot
%	disp(nrms(wk(:), arg.wtrue(:)));
% end
%
% end



% function [wk err time] = restore_dualal_qs(y, mask, cfft, A, C, arg)
%
% bfft = A.arg.psf_fft;
%
% d = C*(A'*y);
%
% ak = A*(C'*d);
% wk = d;
% muk = arg.mu;
%
% Hk = sum(abs(cfft).^2,3).*abs(bfft).^2;
% Iv = ones(C.idim);
%
%
% flag_converged = 0;
% err = zeros(arg.niter+1,1);
% time = zeros(arg.niter+1,1);
%
% err(1) = 20*log10(nrms(wk(:), arg.wtrue(:)));
% time(1) = 0;
%
% tm = cputime;
%
% for ii=1:arg.niter
%
%	if mod(ii,arg.skip/10) == 0 && arg.chat
%		printm('iter: %d',ii)
%	end
%
%	if ~flag_converged
%	qk = shrink(wk + muk*(C*(A'*ak)), arg.lambda*muk);
%	gk = ak - y + A*(C'*qk);
%
%	nrg = sum(abs(vec(gk)).^2);
%	nrw = sum(abs(vec(qk-wk)).^2)/muk;
%
%	jj = 0;
% %	cgiter = arg.cgiter;
%	while nrg > nrw %&& jj <= arg.aliter
%
% %		qk = shrink(wk + muk*(C*(A'*ak)), arg.lambda*muk);
% %		gk = ak - y + A*(C'*qk);
%
% %		Pk = Iv + muk*Hk;
% %		pk = ifft2(fft2(gk)./Pk);
%		pk = gk/(1+muk);
%
%		ak = ak - pk;
%
%		qk = shrink(wk + muk*(C*(A'*ak)), arg.lambda*muk);
%		gk = ak - y + A*(C'*qk);
%
%		nrg = sum(abs(vec(gk)).^2);
%		nrw = sum(abs(vec(qk-wk)).^2)/muk;
%
%		jj = jj+1;
%		if (mod(jj,500) == 0 || nrg < nrw) && arg.chat
%			printm('inner iter: %d',jj)
%			printm('diff: %e', nrg-nrw)
%		end
%
%		if nrg-nrw > 0 && abs(nrg-nrw) < arg.altol
%			printm('Reached convergence at iter: %d, inner iter: %d\n', ii, jj);
%			flag_converged = 1;
%			break;
%		end
%	end
%
%	wk = qk;
% %	wk = shrink(wk + muk*(C*(A'*(T'*ak))), arg.lambda*muk);
%	muk = 2*muk;
%	end
%
%	time(ii+1) = cputime - tm;
%	err(ii+1) = 20*log10(nrms(wk(:), arg.wtrue(:)));
% end
%
% disp(cputime-tm);
%
% if arg.plot
%	figure(1), im(arg.wtrue,'Original coefficients','cbar',arg.clim);
%	figure(2), im(wk,'Restored coefficients','cbar',arg.clim);
%	figure(3), im(wk-arg.wtrue,'Difference coefficients','cbar');
% end
%
% if arg.chat && arg.plot
%	disp(nrms(wk(:), arg.wtrue(:)));
% end
%
% end



% function [xr isnr err time itrtime] = restore_salsa2(y, mask, A, C, arg)
%
% if strcmpi(arg.model_type, 'non-circ')
%	arg.do_cg = 1;
% else
%	arg.do_cg = 0;
% end
%
% if arg.reg_approx
%	dg = 1;
% else
%	dg = ones(C.odim);
%	dg(:,:,end) = 0;
% end
%
% bfft = A.arg.psf_fft;
% T = Gdiag2(ones(size(y)),'mask_in',mask);
%
% % Ha = AC'CA' and C'C = I
% Ha = abs(bfft).^2;
% % Hb = A'(I - (AC'CA' + mI)^-1*AC'CA')
% Hb = conj(bfft).*(1 - Ha./(Ha + arg.mu));
% % Hc = A'(AC'CA' + mI)^-1*A
% Hc = abs(bfft).^2./(Ha + arg.mu);
%
% % d = CA'(I - (AC'CA' + mI)^-1*AC'CA')y/mu
% d = C*ifft2(Hb.*fft2(T'*y))/arg.mu;
% dd = C*(A'*(T'*y));
%
% printm('Cond number of Hmu: %f', 1+1/arg.mu);
%
% xr = A'*(T'*y);
% wk = C*xr;
% eta = zeros(C.odim);
%
% isnr = zeros(arg.niter+1,1);
% err = zeros(arg.niter+1,1);
% time = zeros(arg.niter+1,2);
% itrtime = zeros(arg.niter+1,2);
%
% isnr(1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xr(mask)-arg.xtrue(:)));
% err(1) = 20*log10(nrms(xr(mask), arg.xtrue(:)));
%
% tm1 = cputime;
% tm2 = tic;
%
% for ii=1:arg.niter
%
%	if mod(ii,arg.skip) == 0 && arg.chat
%		printm('iter: %d',ii)
%	end
%
%	uk = shrink(wk + eta, dg*arg.lambda/arg.mu);
% %	eta = eta - uk + wk;
%
%	% 1/2||y-TAC'w|| + lam||u|| + mu/2||u-w-eta||
%	% wkp1 = (CAT'TAC' + muI)^-1*(CA'T'y + mu(u-eta))
%
% %	wk = d + uk - eta - CA'(AC'CA' + mI)^-1*AC'(uk-eta);
%	if arg.do_cg
%		yd = dd + arg.mu*(uk-eta);
%		[wk tm] = inner_cg(yd, wk, A, T, C, 1, arg);
%		itrtime(ii+1,1) = tm(1);
%		itrtime(ii+1,2) = tm(2);
%	else
%		wk = d + uk - eta - C*ifft2(Hc.*fft2(C'*(uk-eta)));
%	end
%
%	xr = C'*wk;
%	wk = C*xr;
%
%	eta = eta - uk + wk;
%
%	time(ii+1,1) = cputime - tm1;
%	time(ii+1,2) = toc(tm2);
%	err(ii+1) = 20*log10(nrms(xr(mask), arg.xtrue(:)));
%	isnr(ii+1) = 20*log10(norm(y(:)-arg.xtrue(:))/norm(xr(mask)-arg.xtrue(:)));
%
% end
%
% if arg.plot
%	figure(1), im(arg.xtrue,'Original image','cbar',arg.clim);
%	figure(2), im(y,'Blured noisy image','cbar',arg.clim);
%	figure(3), im(xr,'Restored image','cbar',arg.clim);
%	figure(4), im(xr-arg.xtrue,'Difference image','cbar');
% end
%
% if arg.chat
%	printm(['\nCPU time: %6.4f\tActual time: %6.4f\nNRMSE (db): %6.2f\t' ...
%		'NRMSE: %6.4f\nISNR: %6.4f\n'], time(end,1), time(end,2), err(end),...
%		nrms(xr(mask), arg.xtrue(:)), isnr(end));
% end
%
% end







% function [wk time] = inner_cg2(y, wi, A, T, C, pfft, arg)
%
% wk = wi;
% Gx = C*(A'*(T'*(T*(A*(C'*wk)))));
% Rx = arg.mu*wi;
%
% ddir = 0;
% oldip = inf;
%
% tm1 = cputime;
% tm2 = tic;
%
% for iter = 1:arg.cgiter
%
%	ngrad = y - Gx - Rx;
%
%	if arg.precon
%		pregrad = real(ifft2(pfft.*fft2(ngrad)));
%	else
%		pregrad = ngrad;
%	end
%
%	newip = ngrad(:)' * pregrad(:);
%
%	gamma = newip / oldip;	% Fletcher-Reeves
%	ddir = pregrad + gamma * ddir;
%
%	oldip = newip;
%
%	Gdir = C*(A'*(T'*(T*(A*(C'*ddir)))));
%	Rdir = arg.mu*ddir;
%
%	denom = ddir(:)' * vec(Gdir + Rdir);
%
%	if denom == 0
%		warning('inner_cg:denom','found exact solution??? step=0 now!?');
%		step = 0;
%	else
%		step = (ddir(:)' * ngrad(:)) / denom;
%		step = real(step); % real step seems only natural
%	end
%
%	if step < 0
%		warning('inner_cg:step','not descent direction, downhill?');
%	end
%
%	Gx = Gx + step * Gdir;
%	Rx = Rx + step * Rdir;
%	wk = wk + step * ddir;
%
%	if norm(step*ddir(:)) < norm(wk(:))*arg.cg_tol
% %		printm('Convergence reached at iteration %d', iter);
%		break
%	end
% end
%
% time = [cputime-tm1 toc(tm2)]/iter;
%
% end



% function u = proj(d, lam)
%
% u = sign(d).*min(abs(d),lam);
%
% end


function image_restore_test

end
