 function [x, info] = ir_mri_barista_undecwavelet(data, mask, smap, kmask, varargin)
%function [x, info] = ir_mri_barista_undecwavelet(data, mask, smap, kmask, [options])
%|
%| Diagonal-majorized FISTA for MRI. This function is basically a wrapper
%| for ir_barista_analysis that sets up the wavelets, sensitivity maps,
%| and Fourier encoding. Use this for undecimated wavelets, see
%| ir_mri_barista_wavelet.m for orthogonal wavelets.
%|
%| tested for undecimated Haar wavelets
%|
%| cost(x) = (y-Ax)'W(y-Ax) / 2 + beta*R(W*x), R admits shrinkage
%|
%| in
%|	data	[nd 1]		raw data
%|	mask	[(N)]		logical support mask (only tested for 2D)
%|	smap	[(N) nc]	sensitivity maps
%|	kmask	[(N)]		k-space sampling mask (only Cartesian tested)
%|
%| options
%|	beta	double		regularization parameter (default: 2^-2)
%|	niter	int		# total iterations (default: 30, use stop_diff_tol)
%|	x0	[np 1]		initial estimate (default: zeros((N))
%|	nlevels int		number of wavelet levels (default: 2)
%|	wname	str		wavelet type (default: haar)
%|
%| out
%|	x	[(N)]		estimated result
%|	info	struct		random information (for debugging/testing)
%|
%| Copyright Dec 2013, Matthew Muckley, University of Michigan

if nargin < 4, help(mfilename), error(mfilename), end

arg.x0 = [];
arg.niter = 1;
arg.beta = 2^-2;
arg.xinf = [];
arg.wname = 'haar';
arg.nlevels = 2;

arg = vararg_pair(arg, varargin);

% fill out the rest of the parameters
if isempty(arg.x0)
	warn 'no initialization, using zeros';
	arg.x0 = zeros(size(mask));
end

arg.nk = numel(mask);

nx = size(mask, 1);
ny = size(mask, 2);

% estimate majorizer
coilind = numel(size(smap));
d = sum(abs(smap).^2,coilind);
if strcmp(arg.wname, 'haar')
	b = ir_wavmajundec(1./d, arg.nlevels, 2);
else
	error('Not coded for non-Haar wavelets');
end

% build matrices
nc = numel(size(smap));
nc = size(smap, nc);
Q = (1/sqrt(nx*ny))*Gdft('ifftshift', 1, 'fftshift', 1, ...
	'samp', kmask);
S = cell(nc,1);
F = cell(nc,1);
for i=1:nc
	S{i} = Gdiag(col(smap(:,:,i)));
	F{i} = Q;
end
S = block_fatrix(S, 'type', 'col');
F = block_fatrix(F, 'type', 'diag');
A = F*S;
W = Gwave2('mask', true(size(mask)), 'wname', arg.wname, ...
	'nlevel', arg.nlevels);

[x, info] = ir_barista_analysis(data, mask, A, W, d, b(:), 'niter', ...
	arg.niter, 'x0', arg.x0, 'xinf', arg.xinf, 'beta', arg.beta);
