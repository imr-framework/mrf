 function x = pcg_bls(A, C, yi, mi, xi, beta, niter, varargin)
%function x = pcg_bls(A, C, yi, mi, xi, beta, niter, varargin)
%
% Minimize the cost function with regularizer 1 by using preconditioned conjugate
% gradient with backtracking line search method.
% cost(xi) = |yi-Ami.*exp(xi)|_W^2 + beta * |C xi|^2
%
% in
%	A	[Nd Np]		system fatrix
%	C	[Nc*Np Np]	finite differencing fatrix
%	yi	[Nd 1]		k-space data
%	mi	[Np 1]		magnitude image
%	xi	[Np 1]		phase image
%	beta	[1]		regularization parameter
%	niter	[1]		# of iterations for updating x
%
% option
%	Pre	[1] or [Np Np]	preconditioner
%	thr	[1]		stop criterion for line search
%
% out
%	x	[1]		updated phase image
%
% Copyright 2012-06-15, Feng Zhao, The University of Michigan
% Written based on pl_pcg_qs_ls.m in Image Reconstruction Toolbox by
% Jeffrey Fessler, The University of Michigan
% 2013-04-11 modified to enable octave by JF

arg.thr = 0.01;
arg.Pre = 1;
arg.btr = 1; % do backtracking or not
arg.mask = true(size(xi));
arg.chat = 0;
arg = vararg_pair(arg, varargin);

oldinprod = 0;

for iter = 1:niter
	if arg.chat
		ticker(mfilename, iter, niter)
	end

	% update xi by roughness regularization
	EX_xi = exp(1i*xi);
%	Am0 = A * diag_sp(mi);
	Am0 = A * Gdiag(mi, 'mask', arg.mask); % jf
	Am0_xn = Am0 * EX_xi;
	GLn = 2*real( 1i*conj(EX_xi) .* (Am0'*(yi - Am0_xn)) );
	GRn = 2*beta*C*xi; % gradient for beta1*R(x(n))
	Gn = GLn + C'*GRn; % gradient for Psi(x(n))
	pre_G = arg.Pre * Gn; % preconditioned gradient
	newinprod = Gn' * pre_G;

	% update search direction
	if oldinprod == 0
		dn = -pre_G;
		gamma = 0;
	else
		gamma = newinprod / oldinprod; % Fletcher-Reeves
		dn = -pre_G + gamma * dn;
	end
	oldinprod = newinprod;

	% check if descent direction
	if dn' * Gn > 0
		warning 'wrong direction so resetting'
		printm('<ddir,grad>=%g, |ddir|=%g, |grad|=%g', ...
				dn' * Gn, norm(dn), norm(Gn))
		% reset
		dn = -pre_G;
		oldinprod = 0;
	end

	pn = C * dn;

	% optimize step size ak
	ak = 0;
	ak_old = 1e5;
	k = 1;

	% jf added ak==0:
	while ((ak==0) || (abs(ak-ak_old)/abs(ak) > arg.thr)) && (k <=4)
		ak_old = ak;
		xi_dn = xi + ak*dn;
		EX_xi_dn = exp(1i*(xi_dn));
		Am0_xn_dn = Am0 * EX_xi_dn;
		tmp = 1i*conj(EX_xi_dn) .* (Am0'*(yi - Am0_xn_dn));
		if ak ~= 0
			GLn = 2*real( tmp );
			C_xi_dn = C * xi_dn;
			GRn = 2*beta*(C_xi_dn);
			% gradient of beta1*R(x(n)+ak*d(n))
		else
			C_xi_dn = C*xi;
		end
		Gf = dn'*GLn + pn'*GRn; % gradient of f(a)

		tmp1 = Am0 * (EX_xi_dn .* dn); % jf
		CLn = 2 * real(tmp1'*tmp1 + dn'*(-1i*tmp.*dn)); % curvature of L(x(n)+ak*d(n))

		Cf = CLn + 2*beta*(pn'*pn); % curvature of f(a)
		Cf = abs(Cf); % take its absolute value
		stp = Gf/Cf;

		if arg.btr == 1 % backtracking
			R1x = norm(C_xi_dn)^2;
			fa0 = norm(yi-Am0_xn_dn)^2 + beta*R1x;

			stp = 2*stp;
			while 1
				stp = stp/2;
				ak1 = ak - stp;
				xi_dn1 = xi + ak1*dn;
				EX_xi_dn1 = exp(1i*xi_dn1);
				R1x = norm(C*xi_dn1)^2;
				fa1 = norm(yi-Am0*EX_xi_dn1)^2 + beta*R1x;
				if (fa1-fa0) < 1e-8
					break
				end
				fa0 = fa1;
				printm 'backtracking'
			end
		end

		if Cf == 0
			warning 'found exact solution???  step=0 now!?'
			ak = 0;
		else
			ak = ak - stp;
		end

		if ak < 0
			warning 'downhill step?'
		end

		k = k + 1;
	end

	xi = xi + ak*dn; % update x
end
x = xi;
end
