function [z rho rhsn] = doCGMIL(Dinv, bi, z, params)
% 2011-11-30, Sathish Ramani, University of Michigan
% Do CG for solving (Dinv + R * Gamma^-1 * R') z = R * Gamma^-1 * bi
% No preconditioning (Jacobi preconditioner for (diagonal + RGR) may work, but not implemented!)

% Compute rhs = R * Gamma^-1 * bi
Cbi = bi / params.mEAWA; % Update solution
rhs =  params.R * Cbi;

rhsn = sqrt(sum(abs(rhs(:)).^2));

%% Compute Az
Dz = Dinv .* z;
RAdjz = ( params.R )' * z;

GRz = RAdjz / params.mEAWA; % Update solution

RGRz =  params.R * GRz;
Az = Dz + RGRz;

g = rhs - Az; % Gradient of quadratic surrogate / residue of the linear system
rho_old = sum(abs(g(:)).^2);
p = g;
itrCG = 1; % CG iteration index

%% CG Iterations
while(itrCG <= params.MFISTA.nCG)
    % Do Ax on p
    Dz = Dinv .* p;
    RAdjz = ( params.R )' * p;
    
    GRz = RAdjz / params.mEAWA; % Update solution

    RGRz =  params.R * GRz;
    Az = Dz + RGRz;

    % Calculate update scaling factor
    temp = sum(conj(p(:)).*Az(:));
    alphaCG = rho_old/real(temp); %% Adapted to the complex setting <x,y> = sum(x.*y);
    
    % Update estimate and gradient
    z = z + alphaCG*p;
    g = g - alphaCG*Az;
    rho_new = sum(abs(g(:)).^2);
    
    % New residue
    betaCG = rho_new/rho_old;
    
    % Update Search Direction
    p = g + betaCG*p;
    
    % Book-keeping
    rho_old = rho_new;
    itrCG = itrCG + 1;
    
end
rho = sqrt(rho_old);