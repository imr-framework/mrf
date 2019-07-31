function anew = get_MaxEigenValue_1by1(params, strMat)
% 2011-11-30, Sathish Ramani, University of Michigan
% Find the maximum eigenvalue of the system sum_i Si TT Si using power method
% Refer Matthieu's ISBI 2009 paper "WAVELET-REGULARIYED RECONSTRUCTION FOR RAPID MRI"

%% Load Parameters
N = params.N;
tol = params.eigtol;

maxitr = params.eigpowermaxitr;
dispitrnum = params.eigdispitr;

%% Initialization
randn('state', 0);
v = randn(N); v = v/sqrt(sum(v(:).^2));
v = v + ones(N);

anew = 1;
aold = 0;

atol = abs(anew-aold)/abs(anew);

itr = 1;
while ((itr <= maxitr) && (atol > tol))
	% Do A' * W * A * v
	Mv = compute_GradDataTerm_CT2D(v, params);

	%% Compute Eigenvalue update
	aold = anew;
%	anew = real(sum(v(:).*Mv(:)))/sum(abs(v(:)).^2); % wrong if v complex!
	anew = real(sum(conj(v(:)).*Mv(:)))/sum(abs(v(:)).^2); % 2016-03-10 change!

	%% Compute Eigenvector update
	v = Mv./sqrt(sum(abs(Mv(:)).^2));

	%% Book keeping
	itr = itr + 1;

	atol = abs(anew-aold)/abs(anew);

	if(~mod(itr, dispitrnum))
		disp(['Itr=' int2str(itr) '; Eval_' strMat ' = ', num2str(anew) '; adiff = ' num2str(atol) '; Tol = ' num2str(tol)]);
	end
end
disp(['Itr=' int2str(itr) '; Eval_' strMat ' = ', num2str(anew) '; adiff = ' num2str(atol) '; Tol = ' num2str(tol)]);
