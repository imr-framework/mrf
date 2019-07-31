 function [x, info] = ir_barista_synthesis(data, mask, A, R, d, varargin)
%function [x, info] = ir_barista_synthesis(data, mask, A, R, d, [options])
%|
%| Diagonal-majorized (variable shrinkage) FISTA.
%| cost(x) = (y-Ax)'W(y-Ax) / 2 + beta*phi(x), phi admits shrinkage
%|
%| in
%|	data	[nd 1]		raw data
%|	mask	[(N)]		logical support mask
%|	A	[nd (N)]	system matrix
%|	R	[m (N)]		sparsity promoting regularizer, R'R = I
%|	d	[(N) 1]		elements of diagonal matrix that upper bounds RA'AR'
%|
%| options
%|	beta	double		regularization parameter
%|	niter	int		# total iterations (default: 30)
%|	x0	[np 1]		initial estimate in synthesis domain
%|				(default: zeros((N))
%|	alpha	double		restart angle parameter
%|				(default -1*cos(80*pi/180)
%|				recommend not changing)
%|	shrink handle	shrinkage function (can change for non-l1)
%|
%| out
%|	x	[(N)]	the estimated result
%|	info	struct	random information (typically for debugging/testing)
%|
%| Copyright Dec 2013, Matthew Muckley, University of Michigan

if nargin < 5, help(mfilename), error(mfilename), end

arg.x0 = [];
arg.niter = 1;
arg.beta = 2^-2;
arg.shrink = @(t, a) (t-a .* sign(t)) .* (abs(t) > a);
arg.xinf = [];
arg.restart = 1;
arg.timeonly = 0;
arg.obarista = 1;
arg.alpha = -1*cos(80*pi/180);
arg.norm_diff_tol = 1e-6;

arg = vararg_pair(arg, varargin);

% check for initialization
if isempty(arg.x0)
	warn 'no initialization, using zeros';
	arg.x0 = R*zeros(size(mask));
end

arg.x0 = R*arg.x0;

% build diagonal majorizer
Dinv = Gdiag(1./d);

% fista parameters
tau = 1;

% vectorize everything for simplicity
data = col(data); arg.x0 = col(arg.x0);

% initializations
% arg.x0(1./d > 1e6*median(1./d)) = 0;
x = arg.x0;
z = x;

% check how fast we go
info.time(1) = 0;
if ~isempty(arg.xinf)
	info.dist(1) = norm(col(R'*x) - col(arg.xinf(mask)));
end

tic;

loopflag = true;
i = 0;

% start the iterations
while loopflag
	i = i + 1;
	ticker(mfilename, i, arg.niter)

	ngrad = R*(A'*(data - A*(R'*z)));

	xold = x;
	oldtau = tau;
	tau = (1/2)*(1+sqrt(1+4*tau^2));

	x = arg.shrink(z + Dinv*ngrad, arg.beta ./ d);

	xdiff = x - xold;

	%restart check
	if arg.obarista
		zdiff = z - x;

		if real(zdiff'*xdiff) / ...
			sqrt(real(zdiff'*zdiff)*real(xdiff'*xdiff)) > arg.alpha
			z = x;
			tau = 1;
		else
			z = x + ((oldtau-1)/tau)*xdiff + (oldtau/tau)*(x - z);
		end
		clear zdiff;
	elseif arg.restart
		zdiff = z - x;

		if real(zdiff'*xdiff) / ...
			sqrt(real(zdiff'*zdiff)*real(xdiff'*xdiff)) > arg.alpha
			z = x;
			tau = 1;
		else
			z = x + ((oldtau-1)/tau)*xdiff;
		end
		clear zdiff;
	else
		z = x + ((oldtau-1)/tau)*xdiff;
	end

	normdifftest = norm(xdiff)/norm(xold);
	clear xdiff;

	info.time(i+1) = toc;
	if ~isempty(arg.xinf)
		if ~arg.timeonly
			info.dist(i+1) = norm(col(R'*x) - col(arg.xinf(mask)));
		end
	end
	if arg.niter ~= 1
		if i == arg.niter
			loopflag = false;
		end
	else
		if normdifftest < arg.norm_diff_tol
			loopflag = false;
		end
	end
end

x = embed(R'*x, mask);
