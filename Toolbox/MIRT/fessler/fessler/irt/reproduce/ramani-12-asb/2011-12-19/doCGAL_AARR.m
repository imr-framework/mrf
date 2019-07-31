function [xnew rho rhsn itrCG] = doCGAL_AARR(rhs, xnew, params)
% 2011-11-30, Sathish Ramani, University of Michigan
% Load parameters
% fstr = params.formatstringC;
tol = params.CG.restol;
rhsn = sqrt( sum( abs( rhs(:) ) .^ 2 ) );
if(rhsn < 1)
    tol = tol * rhsn;
end

%%  Do Ax to prepare for the CG iteration
DT = ( params.A )' * params.A * xnew;
RT = params.AL.nu1 * params.ig.mask .* real( ifft2( params.RR .* fft2( params.ig.mask .* xnew ) ) ); % Accounts for the fact that R has a mask
grd = DT + RT;

%% Prepare for the CG iterations
grd = rhs - grd; % Gradient of quadratic surrogate / residue of the linear system
if params.CG.precon % SIA preconditioner
    u = real( ifft2( fft2( grd ) ./ ( params.CAA + params.AL.nu1 * params.RR) ) ); 
else
    u = grd; % Direction vector
end
tau_old = sum( u(:) .* grd(:) );
rho = sqrt( sum( abs( grd(:) ) .^ 2 ) );
pgrd = u;
itrCG = 1; % CG iteration index

%% CG Iterations
while((itrCG <= params.CG.nCG) && (rho >= tol))
    % Do Ax on pgrd
    DT = ( params.A )' * params.A * pgrd;
    RT = params.AL.nu1 * params.ig.mask .* real( ifft2( params.RR .* fft2( params.ig.mask .* pgrd ) ) ); % Accounts for the fact that R has a mask

    Mp = DT + RT;
    
    % Calculate update scaling factor
    temp = sum( pgrd(:) .* Mp(:) );
    alphaCG = tau_old / temp; %% Adapted to the complex setting <x,y> = sum(x.*y);
    
    % Update estimate and gradient
    xnew = xnew + alphaCG * pgrd;
    grd = grd - alphaCG * Mp;
    rho = sqrt(sum(abs( grd(:) ) .^ 2 ) );
    
    % New residue
    if params.CG.precon % SIA preconditioner
        u = real( ifft2( fft2( grd ) ./ ( params.CAA + params.AL.nu1 * params.RR) ) ); 
    else
        u = grd; % Direction vector
    end
    
    % Norm of new residue
    tau_new = sum( u(:) .* grd(:) );
    betaCG = tau_new / tau_old;
    
    % Update Search Direction
    pgrd = u + betaCG * pgrd;
    
    % Book-keeping
    tau_old = tau_new;

    itrCG = itrCG + 1;
end
itrCG = itrCG - 1;

xnew = xnew .* params.ig.mask; % Account for mask