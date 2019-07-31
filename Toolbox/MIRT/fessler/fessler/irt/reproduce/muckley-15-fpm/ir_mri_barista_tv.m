 function [x, info] = ir_mri_barista_tv(data, mask, smap, kmask, varargin)
%function [x, info] = ir_mri_barista_tv(data, mask, smap, kmask, [options])
%|
%| Diagonal-majorized (variable shrinkage) FISTA for MRI. This function is
%| basically a wrapper for ir_barista_analysis that sets up the finite
%| differencing matrix and the De Pierro majorizer.
%|
%| cost(x) = (y-Ax)'W(y-Ax) / 2 + beta*R(W*x), R admits shrinkage
%|
%| in
%|	data	[nd 1]	the raw data
%|	mask	[(N)]	logical support mask (only tested for 2D)
%|	smap	[(N) nc] sensitivity maps
%|	kmask	[(N)]	k-space sampling mask (only tested for Cartesian)
%|
%| options
%|	beta	double	regularization parameter (default: 2^-2)
%|	niter	int	# total iterations (default: 30)
%|	x0		[np 1]	initial estimate (default: zeros((N))
%|	nlevels int	number of wavelet levels (default: 3)
%|	wtype	str	wavelet type (default: haar)
%|
%| out
%|	x	[(N)]	the estimated result
%|	info	struct	random information (typically for debugging/testing)
%|
%| Copyright Dec 2013, Matthew Muckley, University of Michigan

if nargin < 4, help(mfilename), error(mfilename), end

arg.x0 = [];
arg.niter = 1;
arg.beta = 2^-2;
arg.xinf = [];

arg = vararg_pair(arg, varargin);

% fill out the rest of the parameters
if isempty(arg.x0)
	warn 'no initialization, using zeros';
	arg.x0 = zeros(size(mask));
end

arg.nk = numel(mask);

nx = size(mask, 1);
ny = size(mask, 2);

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
C = Cdiffs(size(mask));
% C = [C.arg.Cc{1}.arg.matrix; C.arg.Cc{2}.arg.matrix; ...
%	C.arg.Cc{3}.arg.matrix; C.arg.Cc{4}.arg.matrix];

% estimate majorizer
coilind = numel(size(smap));
d = col(sum(abs(smap).^2,coilind));
% d = max(d)*ones(size(d));
Dinv = Gdiag(1./d);
b = abs(C)*(Dinv*(abs(C)'*ones(size(C,1),1)));


% use maxprec of -6 since Cdiffs is single precision
[x, info] = ir_barista_analysis(data, mask, A, C, d, b, 'niter', ...
	arg.niter, 'x0', arg.x0, 'xinf', arg.xinf, 'beta', arg.beta, ...
	'maxprec', -6);
