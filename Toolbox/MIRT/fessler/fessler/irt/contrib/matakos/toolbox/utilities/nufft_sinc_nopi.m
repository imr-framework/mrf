  function y = nufft_sinc_nopi(x)
%|function y = nufft_sinc_nopi(x)
%|
%| no pi version of "sinc" function, because matlab's sinc() is in a toolbox
%|
%| Copyright 2001-12-8, Jeff Fessler, University of Michigan
%| Modified by M Allison to not have a pi.

if nargin < 1, help(mfilename), error(mfilename), end
if streq(x, 'test'), nufft_sinc_test, return, end

iz = find(x == 0); % indices of zero arguments
x(iz) = 1;
y = sin(x) ./ (x);
y(iz) = 1;


% test
function nufft_sinc_test

x = linspace(-4, 4, 2^26+1)';

nufft_sinc(0); % warm up
cpu etic
y1 = nufft_sinc(x);
cpu etoc 'nufft_sinc time'

if 2 == exist('sinc')
	sinc(0); % warm up
	cpu etic
	y2 = sinc(x);
	cpu etoc 'matlab sinc time'
	jf_equal(y1, y2)
end

if im, plot(x, y1, '-'), end
