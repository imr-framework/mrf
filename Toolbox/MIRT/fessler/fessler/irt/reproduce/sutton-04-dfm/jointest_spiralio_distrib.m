function [fm] = jointest_spiralio(data,kx1,ky1,kx2,ky2,dts,tau)

% code used in 2004 MRM paper sutton:04:dfm http://dx.doi.org/10.1002/mrm.20079

%--------------------------------------------
% jointestIFR.m
% created November 2001
% Last modified June 28, 2004
% Brad Sutton, bsutton@uiuc.edu
% Work done at University of Michigan
% Added comments on 4/24/2013 BPS
%
% Code is to jointly estimate image and field map from spiral-in/ spiral-out
% data. Needs the Fessler toolbox.
%
% Paper is: Sutton BP, Noll DC, Fessler JA.
%		Dynamic field map estimation using a spiral-in/spiral-out
%		acquisition. Magn Reson Med. 2004 Jun;51(6):1194-204.
%
%
%
% Relevant publications:
% Sutton, B. P. and S. J. Peltier and J. A. Fessler and D. C. Noll.
%	"Simultaneous estimation of I0, R2*, and field map using a
%	multi-echo spiral acquisition." 10th Int. Society of Mag. Res.
%	Med., 1323, 2002.
%
% J A Fessler, B P Sutton. Nonuniform fast Fourier transforms using
%	min-max interpolation. IEEE Tr. Sig. Proc., 51(2):560-574, 2003.
%
% B. P. Sutton, D. C. Noll, J. A. Fessler. Fast, iterative,
%	field-corrected image reconstruction for MRI. IEEE Tr. Med. Im.,
%	22(2):178-188, 2003.
%

% Parameters for time segmented approach.
N = 64;
numloop = 50; % number of outer iteration loops
L = 9; % number of time segments
FOV = 24; % in cm
tau1 = 3e-3;

%Other variables used
datseg1 = length(kx1);
datseg2 = length(kx2);
npm = N*N;
niter = 10;
tt1 = [1:datseg1]'*dts+tau1;
tt2 = [1:datseg2]'*dts+max(tt1(:))+tau;
tt_ext = [tt1;tt2];

%we_b holds the initial estimates for the field
we_b = col(zeros(N)) + col(10*2*pi*randn(64));

%Saving intermediate results
costresult = zeros(1,numloop);
imgresult = zeros(N*N,numloop);
weresult = zeros(N*N,numloop);

for loop = 1:numloop
	if ((loop<4)|(~mod(loop,10))) % Only recompute interpolators for
				% for iterations with large changes
				% in parameters.
		A1 = fast_mr(kx1,ky1,FOV,N,2*N,5,tt1,we_b(:),1,L,1,[],0);
		A2 = fast_mr(kx2,ky2,FOV,N,2*N,5,tt2,we_b(:),1,L,1,[],0);
		A_full = spiral_inout(A1,A2);
	else
		% create mri object
		A1 = fast_mr(A1,'we',we_b(:));
		A2 = fast_mr(A2,'we',we_b(:));
		A_full = spiral_inout(A1,A2);
	end

	mask = zeros(N,N);
	C = C2sparse('leak', mask, 8, 0);
	C = sqrt(2^(18)) * C; % IMPORTANT PARAMETER TO BE TUNED

	xinit = zeros(npm,1);
	im_rec = qpwls_pcg(xinit(:),A_full,1,data,0,C,1,niter+1,[],0);
	%keyboard
	mskwe = logical(ones(N,N));

	mobj = im_rec(:,end);
	B = mobj;
	nup = 8;

	powbeta = 0; % IMPORTANT PARAMETER TO BE TUNED

	betwe = (2^powbeta)/((2*pi)^2);
	wetemp = zeros(npm,nup+1);
	wetemp(:,1) = we_b;
	Rwe = Rbuild('tight',mskwe,8,betwe,'quad',0,0);
	for jj = 1:nup
		powalp = 0;
		wenow = wetemp(:,jj);
		if 1
			A1 = fast_mr(A1,'we',wenow(:));
			A2 = fast_mr(A2,'we',wenow(:));
			A_full = spiral_inout(A1,A2);
		end
		if (jj == 1)
			resid = data-A_full*B;
			gg = conj(B).*(A_full'*((tt_ext).*resid));
			gradwe = imag(gg)+Rwe.cgrad(Rwe,wetemp(:,jj));
			roughwe = Rwe.penal(Rwe,wetemp(:,jj));
			cost_old = 1/2*(real(resid'*resid))+roughwe
			cost_new = cost_old;
		else
			cost_old = cost_new;
			gg = conj(B).*(A_full'*((tt_ext).*resid));
			gradwe = imag(gg)+Rwe.cgrad(Rwe,wetemp(:,jj));
		end

		while ((cost_new>=cost_old)|(isnan(cost_new)))
			if (powalp <= -100)
				sprintf('cannot decrease cost function')
				break % keyboard
			end
			powalp = powalp -1;
			alp = 2^powalp;
			wetemp(:,jj+1) = wetemp(:,jj)-alp*(gradwe);
			if 1
				A1 = fast_mr(A1,'we',wetemp(:,jj+1));
				A2 = fast_mr(A2,'we',wetemp(:,jj+1));
				A_full = spiral_inout(A1,A2);
			end
			resid = data-A_full*B;
			roughwe = Rwe.penal(Rwe,wetemp(:,jj+1));
			cost_new = 1/2*(real(resid'*resid))+roughwe
		end
	end % end for jj=1:nup
	we_b = wetemp(:,jj+1);

	imgresult(:,loop) = im_rec(:,end);
	weresult(:,loop) = we_b;
	costresult(:,loop) = cost_new;
end

if 0
	save imgresult imgresult
	save weresult weresult
	save costresult costresult
end

fm = we_b;
