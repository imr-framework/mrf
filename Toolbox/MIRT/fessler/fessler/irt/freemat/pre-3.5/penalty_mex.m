 function out = penalty_mex(arg, varargin)
%function out = penalty_mex(arg, varargin)
% m-file interface to penalty_mex.mex*
% Usage for penalty_mex:
%	d = function('diff1,forw1', x, offsets, [lastdim]);
%		d = C * x
%		offsets is int32 array such as [1 nx nx-1 nx+1].
%		Usually lastdim is absent, and x is treated as single 'image',
%		and output d is of size [size(x) length(offsets)].
%		If lastdim is int32, then it specifies which dimension of x
%		is the last one associated with one 'image data'; other dims
%		that follow are further 'realizations', and output d has size
%		[size(x)(1:lastdim) length(offsets) size(s)(lastdim+1:end)].
%
%	d = function('diff1,forw2', ...);
%		d = C.^2 * x
%		like 'diff1,forw1' but with squared c_kj values.
%
%	x = function('diff1,back1', d, offsets, [lastdim]);
%		offsets is int32 array such as [-1 -nx -nx-1 -nx+1]
%		x = C' * d
%		output x is of size [numel(d)/length(offsets)]
%
%	x = function('diff1,back2', wk, offsets, [lastdim]);
%		like diff1,back1 but with squared c_kj values.
%		useful for diagonal of Hessian of penalty.
%		x = C' .^2 * wk
%
%	d = function('wk,tight,1', kappa, offsets, distance_power);
%		kappa is floats, often just a binary support mask in floats.
%		output is [size(kappa) length(offsets)]
%		contains weights {w_k} for each set of differences
%		where both neighbors are within the mask.
%			wk = kappa_j * kappa_k / distance^distance_power
%		Use distance_power=1 for the classical choice, leading
%		to 1/dis |a-b|^2, where dis=sqrt(2) for diagonal neighbors.
%		Use distance_power=2 for the possible improved choice,
%		leading to |(a-b)/dis|^2.
%
%	d = function('wk,leak,1', kappa, offsets, distance_power);
%		same but the 'leak' version that includes mask neighbors.
%
%	Use wk,*,2 for 2nd-order differences.
%
% Jeff Fessler and Samit Basu

if nargin < 1, help(mfilename), error(mfilename), end

if exist('penalty_mex') == 3
	warning 'penalty_mex.m called yet penalty_mex.mex* exists?'
end

if streq(arg, 'import'), penalty_mex_import, return, end

if streq(arg, 'diff1,forw', 10) | streq(arg, 'diff2,forw', 10)
	[out count] = sscanf(arg, 'diff%d,forw%d');
	if count ~= 2, error 'bad arg', end
	order = out(1); power = out(2);

	x = varargin{1};
	offsets = varargin{2};
	if length(varargin) == 3,
		repdim = varargin{3};
	else
		repdim = ndims(x);
	end
	xsize = size(x);
	outsize = [xsize(1:repdim),length(offsets),xsize((repdim+1):end)];
	out = zeros(outsize);
	nn = max(1,prod(xsize(1:repdim)));
	nrep = max(1,prod(xsize((repdim+1):end)));
	penalty_diff_forw(out, x, nn, offsets, length(offsets), nrep, order, power);

if streq(arg, 'diff1,back', 10) | streq(arg, 'diff2,back', 10)
	[out count] = sscanf(arg, 'diff%d,back%d');
	if count ~= 2, error 'bad arg', end
	order = out(1); power = out(2);

	d = varargin{1};
	offsets = varargin{2};
	if length(varargin) == 3,
		repdim = varargin{3};
	else
		repdim = ndims(d);
	end
	dsize = size(d);
	outsize = [dsize(1:repdim),dsize((repdim+2):end)];
	out = zeros(outsize);
	nn = max(1,prod(dsize(1:repdim)));
	nrep = max(1,prod(dsize((repdim+2):end)));
	penalty_diff_back(out, d, nn, offsets, length(offsets), nrep, order, power);

elseif streq(arg, 'wk,tight', 8)
	if streq(arg, 'wk,tight,1'), order = 1;
	elseif streq(arg, 'wk,tight,2'), order = 2;
	else, error 'bad arg', end

	kappa = varargin{1};
	offsets = varargin{2};
	distance_power = varargin{3};
	outsize = [size(kappa) length(offsets)];
	out = zeros(outsize);
	penalty_diff_wk_tight(out, kappa, size(kappa), ...
		ndims(kappa), offsets, length(offsets), distance_power, order);

elseif streq(arg, 'wk,leak', 7)
	if streq(arg, 'wk,leak,1'), order = 1;
	elseif streq(arg, 'wk,leak,2'), order = 2;
	else, error 'bad arg', end

	kappa = varargin{1};
	offsets = varargin{2};
	distance_power = varargin{3};
	outsize = [size(kappa) length(offsets)];
	out = zeros(outsize);
	penalty_diff_wk_leak(out, kappa, size(kappa), ...
		ndims(kappa), offsets, length(offsets), distance_power, order);


else
	error 'unknown argument'

end


%
% penalty_mex_import()
%
function penalty_mex_import

% locate directory containing the .a file
if streq(computer, 'MAC')
	libdir = [path_find_dir([filesep 'data']) filesep '../mex/lib,ppc']
elseif streq(computer, 'GLNX86') | streq(computer, 'UNIX')
	libdir = [path_find_dir([filesep 'data']) filesep '../mex/lib,i686']
else
	error 'unknown platform'
end
libname = [libdir filesep 'penalty,diff.so'];

type1f = ['float[nn * noffset * nrep] &po, float[nn * nrep] pi, int32 nn, ' ...
	'int32[noffset] offset, int32 noffset, int32 nrep, int32 order, int32 power'];
type1b = ['float[nn * nrep] &po, float[nn * noffset * nrep] pi, int32 nn, ' ...
	'int32[noffset] offset, int32 noffset, int32 nrep, int32 order, int32 power'];
type2 = ['float[prod(dim_i) * noffset] &po, float[prod(dim_i)] kappa, ' ...
	'int32[nod_i] dim_i, int32 nod_i, int32[noffset] offset, ' ...
	'int32 noffset, double distance_power, cint order'];

myimport(libname, 'penalty_diff_forw', 'void', type1f);
myimport(libname, 'penalty_diff_back', 'void', type1b);
myimport(libname, 'penalty_diff_wk_tight', 'int32', type2);
myimport(libname, 'penalty_diff_wk_leak', 'int32', type2);


function myimport(libname, name, type_return, type_call)
import(libname, name, name, type_return, type_call);
