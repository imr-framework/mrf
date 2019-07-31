  function ob = Gdiag2(diag, varargin)
%|function ob = Gdiag2(diag, [options])
%|
%| Construct diagonal "sparse" object, that does D * x = diag(d) * x via d .* x.
%|
%| in
%|	diag	[]	numeric array
%|
%| option
%|	'mask_in'	logical	with sum(mask(:)) = numel(diag) creates short-fat matrix
%|				wih smaller size output
%|	'mask_out'	logical	with sum(mask(:)) = numel(diag) creates tall-skinny
%|				matrix with larger size output
%|
%| out
%|	ob	fatrix2: [numel(mask) numel(diag)] for input or 
%|					[numel(diag) numel(mask)] for output 
%|
%| Copyright 2011-10-31, Antonis Matakos, University of Michigan
%| Based on Gdiag:
%| Copyright 2005-4-6, Jeff Fessler, University of Michigan

if nargin == 1 && streq(diag, 'test'), Gdiag_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

arg.mask_in = [];
arg.mask_out = [];

arg = vararg_pair(arg, varargin);

if ~isempty(arg.mask_in)
	arg.mask_out = true(size(diag));
elseif ~isempty(arg.mask_out);
	arg.mask_in = true(size(diag));
end

arg.diag_in = embed(diag(:), arg.mask_in);
arg.diag_out = embed(diag(:), arg.mask_out);

arg.idim = size(arg.diag_in);
if arg.idim(end) == 1, arg.idim = arg.idim(1:end-1); end % 1d trick
arg.odim = size(arg.diag_out);
if arg.odim(end) == 1, arg.odim = arg.odim(1:end-1); end % 1d trick

% forw = @(arg,x) embed(bsxfun(@times, arg.diag_in, x), arg.mask_in);
% back = @(arg,y) embed(bsxfun(@times, conj(arg.diag_out), y), arg.mask_in);

ob = fatrix2('arg', arg, 'does_many', 1, 'imask', true(size(arg.mask_in)), ...
	'omask', true(size(arg.mask_out)), 'forw', @Gdiag_forw, ...
	'back', @Gdiag_back, 'abs', @Gdiag_abs, 'power', @Gdiag_power, ...
	'gram', @Gdiag_gram, 'sparse', @Gdiag_sparse);


function y = Gdiag_forw(arg, x)

y = bsxfun(@times, arg.diag_in, x);
y = embed(y(arg.mask_in),arg.mask_out);


function x = Gdiag_back(arg, y)

x = bsxfun(@times, conj(arg.diag_out), y);
x = embed(x(arg.mask_out),arg.mask_in);



% Gdiag_abs(): |D|
function ob = Gdiag_abs(ob)
ob.arg.diag_in = abs(ob.arg.diag_in);
ob.arg.diag_out = abs(ob.arg.diag_out);


% Gdiag_power(): D .^ p
function ob = Gdiag_power(ob, p)
ob.arg.diag_in = ob.arg.diag_in .^ p;
ob.arg.diag_out = ob.arg.diag_out .^ p;


% Gdiag_sparse()
function sp = Gdiag_sparse(ob)
d = double(ob.arg.diag_in(ob.arg.imask));

np = size(ob);
idrow = 1:np(1);
idrow = idrow(arg.omask);
idcol = 1:np(2);
idcol = idcol(arg.imask);

sp = sparse(idrow, idcol, d, np(1), np(2), length(d));


% Gdiag_gram(): D'*W*D
% works only for fatrix2 or for Fatrix in case of all true mask
function [ob reuse] = Gdiag_gram(ob, W, reuse, varargin)
if isempty(W)
	ob = Gdiag(ob.arg.diag_in(ob.arg.imask),'mask',ob.arg.imask);
% 	ob.arg.diag_in = abs(ob.arg.diag_in).^2;
% 	ob.arg.diag_out = abs(ob.arg.diag_out).^2;
else
	fail('W not implemented yet')
end


%
% Gdiag_test()
%
function Gdiag_test
% 1D tests

% 1) short-fat
x = (1:10).';
d = ones(8,1);
msk = false(10,1);
msk(2:end-1) = true;

T = Gdiag2(d,'mask_in',msk);

y = T*x;

xh = T'*y;

figure(1), im(xh);
prompt

% figure(1), plot(x);
% figure(2), plot(y);
% figure(3), plot(xh);
% figure(4), im(full(T));
% figure(5), im(full(T'));
% figure(6), im(full(T'*T));
% figure(7), im(full(T*T'));

% 2) long-skinny
x = (1:8).';
d = ones(8,1);
msk = false(10,1);
msk(2:end-1) = true;

T = Gdiag2(d,'mask_out',msk);

y = T*x;

xh = T'*y;

figure(1), im(xh);
prompt

% 2D tests

% 1) short-fat
x = [1 2 4 6 4 2 1]'*[1 2 4 6 4 2 1];
d = ones(4);
msk = false(7);
msk(2:end-2,3:end-1) = true;

T = Gdiag2(d,'mask_in',msk);

y = T*x;

xh = T'*y;

figure(1), im(xh);
prompt

% figure(1), plot(x);
% figure(2), plot(y);
% figure(3), plot(xh);
% figure(4), im(full(T));
% figure(5), im(full(T'));
% figure(6), im(full(T'*T));
% figure(7), im(full(T*T'));

% 2) long-skinny
x = [1 2 4 6 4 2 1]'*[1 2 4 6 4 2 1];
d = ones(7);
msk = false(10);
msk(2:end-2,3:end-1) = true;

T = Gdiag2(d,'mask_out',msk);

y = T*x;

xh = T'*y;

figure(1), im(xh);
prompt



