 function data = poisson(xm, seed)
%function data = poisson(xm, seed)
%
%	Generate Poisson random vector with mean xm.
%	For small, use poisson1.m
%	For large, use poisson2.m
%	see num. rec. C, P. 222
%
%	Copyright 1997-4-29	Jeff Fessler	The University of Michigan


%
%	if no arguments, run a timing race against matlab's poissrnd
%
if nargin == 0
	n = 100;
	t = reshape(linspace(1, 1000, n^2), n, n);
	tic
	poisson(t);
	Sprintf('fessler poisson time %g', toc)
	tic 
	poissrnd(t);
	Sprintf('matlab poissrnd time %g', toc)
return
end

if ~isvar('seed')
	seed = [];
end

data	= xm;
xm	= xm(:);

data(  xm < 12 ) = poisson1(xm(  xm < 12 ), seed);
data(~(xm < 12)) = poisson2(xm(~(xm < 12)), seed);
