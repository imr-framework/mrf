function ob = zero_fatrix(dim)
%| function ob = zero_fatrix(dim)
%|
%| Creates a zero fatrix with the specified dimension.
%|
%| in:
%|	dim	[N 1]	The fatrix dimensions
%|
%| out:
%|	ob			The dim(1)x...xdim(N) zero fatrix
%|
%| Copyright 10-06-02, Antonis Matakos, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end

dim = dim(:)';

if ncol(dim) == 1
	dim = dim([1 1]);
end

arg.dim = dim;

ob = Fatrix(dim,arg,'caller','zero_fatrix','forw',@zero_fatrix_forw,...
	'back',@zero_fatrix_back,'gram',@zero_fatrix_gram);



function y = zero_fatrix_forw(arg, x)

if nrow(x) ~= arg.dim(2)
	error('x size=%d vs dim(2)=%d', nrow(x), arg.dim(2))
end

y = zeros(arg.dim(1),ncol(x));



function x = zero_fatrix_back(arg, y)

if nrow(y) ~= arg.dim(1)
	error('y size=%d vs dim(1)=%d', nrow(y), arg.dim(1))
end

x = zeros(arg.dim(2), ncol(y));


function T = zero_fatrix_gram(ob)

nc = ob.dim(2);

T = zero_fatrix([nc nc]);
