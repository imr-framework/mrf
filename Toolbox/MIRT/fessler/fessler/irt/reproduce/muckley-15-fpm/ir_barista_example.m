%| ir_barista_example.m
%| This is an example of how to run BARISTA. The script simulates a
%| multicoil MRI experiment with autocalibrated SENSE map calculations based
%| on the center of k-space acquired. Then, it estimates the image using
%| sparsity-promoting regularization via a BARISTA. Code optimizations are
%| minimal; this script is primarily meant to be a proof-of-concept.
%| Copyright 2014, Matthew Muckley, The University of Michigan

% generate synthetic data
if ~isvar('xtrue'), printm 'xtrue'
	load kmask;

	nx = 256;
	ny = 256;
	nc = 8;
	beta = 2^-13;

	f.dir = [path_find_dir('mri') '/../data/mri/'];
	f.xtrue = [f.dir 'brainweb_t1.jpg'];
	xtrue = imread(f.xtrue);
	xtrue = xtrue(1:256,1:256);
	xtrue = double(xtrue).';

	smap_true = mri_sensemap_sim('nx', nx, 'ny', ny, 'ncoil', nc);

	mask = abs(xtrue) > 0.1*max(abs(col(xtrue)));
	siglevel = xtrue .* smap_true(:,:,1);
	siglevel = mean(abs(siglevel(mask)));
	nsstd = siglevel/10^(30/20);

	mask = imdilate(mask, true(10,10));

	im plc 3 3
	im(1, xtrue)
	im(2, mask)
	im(3, 'row', 3, abs(smap_true) .* repmat(mask, [1 1 nc]))
prompt
end


% build (initial) system matrix
if ~isvar('A0'), printm 'A0'
	Qm = (1/sqrt(nx*ny)) * Gdft('ifftshift', 1, 'fftshift', 1, ...
		'samp', kmask);
	S = cell(nc,1);
	F = cell(nc,1);
	for i=1:nc
		S{i} = Gdiag(col(smap_true(:,:,i)));
		F{i} = Qm;
	end
	S = block_fatrix(S, 'type', 'col');
	F = block_fatrix(F, 'type', 'diag');
	A0 = F * S;
end


if ~isvar('phantdat'), printm 'k-space data'
	phantdat = A0 * xtrue(:);
	phantdat = phantdat ...
		+ nsstd * (randn(size(phantdat)) + 1i*randn(size(phantdat)));
end


if ~isvar('smap'), printm 'smap'
	lowim = complex(zeros(nx,ny,nc));
	tmp = repmat(kmask, [1 1 nc]);
	lowim(tmp(:)) = double(phantdat);
	centmask = false(nx,ny); centmask(113:144,113:144) = true;
	lowim(~centmask) = 0;
	Q1 = (1/sqrt(nx*ny))*Gdft('ifftshift', 1, 'fftshift', 1, ...
		'mask', true(nx,ny));
	lowim = Q1' * lowim;
	im(4, 'row', 3, abs(lowim))

	tmp = Q1 * xtrue + 1.3 * nsstd * (randn(nx,ny) + 1i*randn(nx,ny));
	tmp(~centmask) = 0; % noisy k-space data for body coil image
	bodyim = Q1' * tmp;

	% make "realistic" smap from noisy data
	smap = mri_sensemap_denoise(lowim, 'bodycoil', double(bodyim), ...
		'chol', 1, 'niter', 1, 'l2b', 5);
	im(6, 'row', 3, abs(smap) .* repmat(mask, [1 1 nc]))
	titlef 'estimate of sensitivity map'
	im(5, sqrt(sum(abs(smap).^2,3)) .* mask, 'SSoS estimate')
prompt
end


if ~isvar('A'), printm 'A'
%	Q = (1/sqrt(nx*ny))*Gdft('ifftshift', 1, 'fftshift', 1, ...
%		'samp', kmask);
	S = cell(nc,1);
	F = cell(nc,1);
	for i=1:nc
		S{i} = Gdiag(smap(:,:,i));
		F{i} = Qm;
	end
	S = block_fatrix(S, 'type', 'col');
	F = block_fatrix(F, 'type', 'diag');
	A = F * S;
end


if ~isvar('xinit'), printm 'xinit'
	xinit = A' * phantdat;
	xinit = reshape(xinit, [nx ny]);
	xinit = xinit ./ sum(abs(smap).^2,3);
	im(8, abs(xinit))
prompt
end


% run without an iteration counter, currently this ends with a normdiff
% cutoff of 1e-6
mask = true(nx,ny);
[xhat, info] = ir_mri_barista_undecwavelet(phantdat(:), mask, smap, kmask, ...
	'x0', xinit, 'beta', beta);
xhat = embed(xhat, mask);
im(9, xhat)
