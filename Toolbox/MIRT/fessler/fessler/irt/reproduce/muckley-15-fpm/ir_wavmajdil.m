 function d = ir_wavmajdil(midd, nlevels, patchsize)
%function d = ir_wavmajdil(midd, nlevels, patchsize)
%|
%| calculates a majorizer for WDW', where D is a diagonal matrix and W is a
%| nlevels orthogonal Haar wavelet transform.
%|
%| Copyright Dec 2013 Matthew Muckley, University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

[nx, ny] = size(midd);
d = zeros(nx,ny);
curnx = nx; curny = ny;
curim = midd;
for i=1:nlevels
	curim = wextend('2D', 'ppd', curim, [patchsize patchsize]);
	curim = imdilate(curim, ones([patchsize patchsize]), 'same');
%	curim = imdilate(curim, ones([patchsize patchsize]), 'same');
	curD = curim(patchsize+2:2:end-patchsize,patchsize+2:2:end-patchsize);
%	curD = imdilate(curD, ones([3 3]), 'same');
	curim = reshape(curD, [curnx/2 curny/2]);

	curD = reshape(curD, [curnx/2 curny/2]);
	d(curnx/2+1:curnx,1:curny/2) = curD;
	d(1:curnx/2, curny/2+1:curny) = curD;
	d(curnx/2+1:curnx, curny/2+1:curny) = curD;
	if (i == nlevels)
		d(1:curnx/2,1:curny/2) = curD;
	end
	curnx = curnx/2; curny = curny/2;
end
