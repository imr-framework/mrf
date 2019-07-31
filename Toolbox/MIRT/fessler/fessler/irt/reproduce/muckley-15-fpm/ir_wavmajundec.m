 function d = ir_wavmajundec(midd, nlevels, patchsize)
%function d = ir_wavmajundec(midd, nlevels, patchsize)
%|
%| Calculates a majorizer for WDW', where D is a diagonal matrix and W is a
%| nlevels stationary wavelet transform. Ignores approximation
%| coefficients
%|
%| Copyright Dec 2013 Matthew Muckley, University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

[nx, ny] = size(midd);
d = zeros(nx,ny,3*nlevels);
bigim = wextend('2D', 'ppd', midd, [nx/2 ny/2]);
for i=1:nlevels
	curim = imdilate(bigim, ...
		ones([patchsize*2^(i-1) patchsize*2^(i-1)]), 'same');

	if i == 1
		curim = curim(nx/2+2:end-nx/2+1,:);
		curim = curim(:,ny/2+2:end-ny/2+1);
	else
		curim = curim(nx/2+3:end-nx/2+2,:);
		curim = curim(:,ny/2+3:end-ny/2+2);
	end
%	curim = wkeep2(curim, [nx ny]);

	d(:,:,i) = curim;
	d(:,:,i+nlevels) = curim;
	d(:,:,i+nlevels*2) = curim;
end
