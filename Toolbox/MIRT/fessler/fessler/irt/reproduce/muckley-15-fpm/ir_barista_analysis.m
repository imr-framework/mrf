 function [x, info] = ir_barista_analysis(data, mask, A, R, d, d_reg, varargin)
%function [x, info] = ir_barista_analysis(data, mask, A, R, d, d_reg, [options])
%|
%| Diagonal-majorized (variable shrinkage) FISTA.
%| cost(x) = (y-Ax)'W(y-Ax) / 2 + beta*phi(x), phi admits shrinkage
%|
%| in
%|	data	[nd 1]	the raw data
%|	mask	[(N)]		logical support mask
%|	A	[nd (N)]	system matrix
%|	R	[m (N)]	real-valued, sparsity promoting regularizer
%|	d	[(N) 1]	elements of diagonal matrix that upper bounds A'A
%|	d_reg	[m 1]	elements of diagonal matrix that upper bounds
%|			Rdiag(1./d)R'
%|
%| options
%|	beta	double	regularization parameter
%|	niter	int	# total iterations (default: inf, so use stop rule)
%|	x0	[np 1]	initial estimate in synthesis domain
%|				(default: zeros((N))
%|	minprec	double	epsilon^(0), initial precision level for prox
%|			subproblem (default 10^-1, recommend not changing)
%|	diffprec double	epsilon_diff, code solves supbroblem to
%|			epsilon_diff*(outer problem step size)
%|			(default 10^-1, recommend not changing)
%|	maxprec	double	epsilon_max, maximum precision level for prox
%|			subproblem (default 10^-12, recommend not changing)
%|	alpha	double	restart angle parameter (default -1*cos(80*pi/180),
%|			recommend not changing)
%|	maxinniter int	maximum number of iterations for prox subproblem
%|				(default 30)
%|
%| out
%|	x	[(N)]	the estimated result
%|	info	struct	random information (typically for debugging/testing)
%|
%| Copyright Dec 2013, Matthew Muckley, University of Michigan

if nargin < 5, help(mfilename), error(mfilename), end

arg.x0 = [];
arg.niter = 1;
arg.beta = 2^-3;
arg.xinf = [];
arg.restart = 1;
arg.obarista = 1;
arg.stopslope = 1/80;
arg.minprec = -1;
arg.maxprec = -12;
arg.maxinniter = 30;
arg.diffprec = -1;
arg.alpha = -1*cos(80*pi/180);
arg.norm_diff_tol = 1e-6;

arg = vararg_pair(arg, varargin);

% check for initialization
if isempty(arg.x0)
	warn 'no initialization, using zeros';
	arg.x0 = zeros(size(mask));
end

% build diagonal majorizer
Dinv = Gdiag(1./d);
Binv = Gdiag(1./d_reg);

% fista parameters
tau = 1;

% vectorize everything for simplicity
data = col(data); arg.x0 = col(arg.x0);

% initializations
x = col(arg.x0);
z = x;
betainv = 1./arg.beta;
betainv(isinf(betainv)) = 0;
q = betainv.*(Binv*(R*x));
threshinds = abs(q) > 1;
q(threshinds) = sign(q(threshinds));

% check how fast we go
info.time(1) = 0;
if ~isempty(arg.xinf)
	info.dist(1) = norm(col(x) - col(arg.xinf));
end

stopval = 10^arg.minprec;
tic;

loopflag = true;
i = 0;

% start the iterations
while loopflag
	i = i + 1;
	ticker(mfilename, i, arg.niter)
	% printm('stopval is %d', stopval);

	ngrad = (A'*(data - A*z));

	xold = x;
	oldtau = tau;
	tau = (1/2)*(1+sqrt(1+4*tau^2));

	b = z + Dinv*ngrad;

	normdiffinner = 5;
	j = 0;
	innertau = 1;
	w = q;
	x = b - Dinv*(R'*(arg.beta.*w));
	x(~mask) = 0;

	while (normdiffinner > stopval) && (j < arg.maxinniter)
		qold = q;
		oldinnertau = innertau;
		innertau = (1/2)*(1 + sqrt(1+4*innertau^2));

		ngrad = R*x;

		% update main variable
		q = w + betainv.*(Binv*ngrad);

		% project on to convex set
		threshinds = abs(q) > 1;
		q(threshinds) = sign(q(threshinds));

		qdiff = q - qold;

		% update momentum variable
		if arg.obarista
			wdiff = w - q;
			if real(wdiff'*qdiff) / ...
				sqrt(real(wdiff'*wdiff)*real(qdiff'*qdiff)) > arg.alpha
				w = q;
				innertau = 1;
			else
				w = q + ((oldinnertau - 1)/innertau)*qdiff + ...
					(oldinnertau/innertau)*(q - w);
			end
			clear wdiff;
		elseif arg.restart
			wdiff = w - q;
			if real(wdiff'*qdiff) / ...
				sqrt(real(wdiff'*wdiff)*real(qdiff'*qdiff)) > arg.alpha
				w = q;
				innertau = 1;
			else
				w = q + ((oldinnertau - 1)/innertau)*qdiff;
			end
			clear wdiff;
		elseif arg.obarista

		else
			w = q + ((oldinnertau - 1)/innertau)*qdiff;
		end

		clear qdiff;

		% calc x for gradient
		innerxold = x;
		x = b - Dinv*(R'*(arg.beta.*w));
		x(~mask) = 0;

		normdiffinner = norm(x - innerxold)/norm(innerxold);
		j = j + 1;
	end

	x = b - Dinv*(R'*(arg.beta.*q));
	x(~mask) = 0;

	xdiff = x - xold;

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

	stopval = max(min(normdifftest*10^arg.diffprec, ...
		10^arg.minprec), 10^arg.maxprec);
	arg.minprec = max(arg.maxprec, log10(stopval));

	clear xdiff;

	info.time(i+1) = toc;
	if ~isempty(arg.xinf)
		info.dist(i+1) = norm(col(x) - col(arg.xinf));
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
