 function [x, info] = ir_mri_barista_wavelet(data, mask, smap, kmask, varargin)
%function [x, info] = ir_mri_barista_wavelet(data, mask, smap, kmask, [options])
%|
%| Diagonal-majorized (variable shrinkage) FISTA for MRI. This function is
%| basically a wrapper for ir_barista_synthesis that sets up the wavelets,
%| sensitivity maps, and Fourier encoding. Use this for orthogonal
%| wavelets, see ir_mri_barista_undecwavelet.m for undecimated wavelets.
%|
%| tested for Haar and Daubechies D4 wavelets
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
%|	niter	int		# total iterations (default: 30)
%|	x0			[np 1]	initial estimate (default: zeros((N))
%|	nlevels int		number of wavelet levels (default: 3)
%|	wtype	str		wavelet type
%|				(default: haar, use db4 for Daubechies D4)
%|	shrink	@		shrinkage function (default: soft thresholding)
%|
%| out
%|	x	[(N)]	the estimated result
%|	info	struct	random information (typically for debugging/testing)
%|
%| Copyright Dec 2013, Matthew Muckley, University of Michigan

if nargin < 4, help(mfilename), error(mfilename), end

arg.x0 = [];
arg.niter = 1;
arg.nlevels = 3;
arg.beta = 2^-2;
arg.wtype = 'haar';
arg.niter = 30;
arg.kmask = [];
arg.shrink = @(t, a) (t-a .* sign(t)) .* (abs(t) > a);
arg.xinf = [];
arg.restart = 1;
arg.L = [];
arg.timeonly = 0;

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
coilind = length(size(smap));
d = sum(abs(smap).^2,coilind); d(~mask) = eps;
if strcmp(arg.wtype, 'haar')
	d = col(ir_wavmajdil(d, arg.nlevels, 2));
elseif strcmp(arg.wtype, 'db4')
	d = col(ir_wavmajdil(d, arg.nlevels, 4));
end
if ~isempty(arg.L)
	d = arg.L*ones(size(d));
end
% d = 2.98*ones(size(d));
% d = sum(abs(smap).^2,coilind);
% d = max(col(d(mask)));

% build matrices
nc = numel(size(smap));
nc = size(smap, nc);
Q = (1/sqrt(nx*ny))*Gdft('ifftshift', 1, 'fftshift', 1, ...
	'samp', kmask, 'mask', mask);
S = cell(nc,1);
F = cell(nc,1);
for i=1:nc
	S{i} = Gdiag(col(smap(:,:,i)));
	F{i} = Q;
end
S = block_fatrix(S, 'type', 'col');
F = block_fatrix(F, 'type', 'diag');
A = F*S;
W = Godwt1(mask, 'level', arg.nlevels, 'wname', arg.wtype);

% don't shrink the approximation coefficients; that would be obnoxious
arg.beta = arg.beta*ones(size(mask));
arg.beta(1:size(mask,1)/2^arg.nlevels,1:size(mask,2)/2^arg.nlevels) = 0;
arg.beta = col(arg.beta);

[x, info] = ir_barista_synthesis(data, mask, A, W, d, 'niter', ...
	arg.niter, 'x0', arg.x0, 'shrink', arg.shrink, 'xinf', arg.xinf, ...
	'beta', arg.beta, 'restart', arg.restart, 'timeonly', arg.timeonly);
