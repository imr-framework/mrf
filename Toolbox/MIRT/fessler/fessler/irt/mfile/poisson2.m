 function data = poisson2(xm, seed)
%function data = poisson2(xm, seed)
% Generate Poisson random column vector with mean xm.
% Uses rejection method - good for large values of xm.
% See "Numerical Recipes in C", P. 222.

if isvar('seed') & ~isempty(seed)
	rand('state', seed)
end

data = zeros(size(xm));
if any(xm < 0), error 'negative poisson means?', end
data(xm > 0) = poisson2_positive(xm(xm > 0));


%
% poisson2_positive()
%
function data = poisson2_positive(xm)
sx = sqrt(2.0 * xm);
lx = log(xm);
gx = xm .* lx - gammaln(1+xm);

data = zeros(size(xm));
id = [1:length(xm)]'; % indicates which data left to do

while any(id)
	Tss = sx(id);
	Tll = lx(id);
	Tgg = gx(id);
	Txx = xm(id);

	tmp = pi * rand(size(id));
	yy = tan(tmp);
	em = Tss .* yy + Txx;
	ib = find(em < 0);

	while ib
		tmp = pi * rand(nrow(ib),ncol(ib));
		yy(ib) = tan(tmp);
		em(ib) = Tss(ib) .* yy(ib) + Txx(ib);
		ib = find(em < 0);
	end

	em = floor(em);
	tt = 0.9 * (1+yy.*yy) .* exp(em .* Tll - gammaln(em+1) - Tgg);
	if (any(tt > 1)), error('ERROR: 0.9 factor is too large!!!!'), end

	ig = rand(size(id)) < tt;
	data(id(ig(:))) = em(ig);
	id = id(~ig(:));
end
