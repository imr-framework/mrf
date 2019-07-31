function [xnew C T E RMSE strStatus] = runMFISTA(data, rhs, xnew, params) 
% (M)FSITA from 
% Beck et al., "Fast gradient-based algorithms for 
% constrained total variation image denoising and deblurring problems,"
% IEEE Trans. Image Processing, vol. 18, no. 11, 2009.
%
% Uses Chambolle-like algo 
% from Figueiredo et al., "Signal Restoration with overcomplete 
% wavelet transforms - comparison of analysis and synthesis priors,"
%
% for solving an intermediate denoising problem
%
% 2011-11-30, Sathish Ramani, University of Michigan
%
% Load Parameters
dispitr = params.dispitr;
dispitrnum = params.dispitrnum;
dispfig = params.dispfig;
Nmask = params.Nmask;

clim = params.clim;

scale = params.scale;

figno = params.figno;

maxitr = params.maxitr;

dxtol = params.dxtol;
dcosttol = params.dcosttol;

xini = xnew;

zr = params.zoomr;
zc = params.zoomc;

fstrC = params.formatstringC;
fstrT = params.formatstringT;
fstrO = params.formatstringO;

%% Compute Error in estimate
err = params.ig.mask .* ( xnew - params.img );
err = sqrt( sum( abs( err(:) ).^2 ) / Nmask );

time_elapse = 0;
itr = 1;

%% Prepare cost and thresholds
cost_new = compute_Cost_CT2D(data, xnew, params); % Reinitialize cost_old for current sigma
cost_old = 2*cost_new;
diffcost = cost_new - cost_old;

normx = sqrt( sum( abs( xnew(:) ) .^ 2 ) );
diffx = normx;

C(itr) = cost_new;
T(itr) = time_elapse;
E(itr) = sqrt( sum( params.ig.mask(:) .* abs( xnew(:) - params.xinf(:) ) .^ 2)) / params.xinfnorm;
RMSE(itr) = err;

%% Show Images
if(dispfig)
    figure(figno); clf;
    subplot(2,2,1); imagesc(scale * params.img(zr,zc), clim); colormap(gray); colorbar; axis off; axis equal; title('Noisefree phantom, x');
    subplot(2,2,2); imagesc(scale * err(zr,zc)); colormap(gray); colorbar; axis off; axis equal; title(cat(2,'x - x', int2str(itr), '; ', num2str(RMSE(itr), fstrC)));
    subplot(2,2,3); imagesc(scale * xini(zr,zc), clim); colormap(gray); colorbar; axis off; axis equal; title('Initial estimate x0');
    subplot(2,2,4); imagesc(scale * xnew(zr,zc), clim); colormap(gray); colorbar; axis off; axis equal; title(cat(2, 'Current estimate x', int2str(itr)));
end

%% Display cost and other parameters
strStatus{itr} = ['MFISTA: ' int2str(itr-1) '/' int2str(maxitr) ...
                  '; C=' num2str(cost_new, fstrC) '; DC=' num2str(diffcost, fstrC) '; Dx/||x||=' num2str(diffx/normx, fstrO) ...
                  '; l2D=' num2str(E(itr), fstrC) '; RMSE=' num2str(err, fstrC) '; ' num2str(time_elapse, fstrT) 's'];
if(dispitr)
    disp('--------------------------------------------------------------------');
    disp(strStatus{itr});
end

%% MFISTA Initialization
tk = 1;
yk = xnew;

Rx = params.R * xnew;
z = zeros(size(Rx));

%% Inner while loop for minimizing cost corresponding to a given sigma
while((itr <= maxitr) && (diffx/normx >= dxtol) && (diffcost <= dcosttol))
    tic;
    
    %% Book Keeping
    cost_old = cost_new;
    xold = xnew;

    %% Compute the gradient of the data term
    AWAy = compute_GradDataTerm_CT2D(yk, params);
    gradyn = AWAy - rhs;
    
    Cyk = params.mEAWA * yk;
    bi =  Cyk - gradyn;
    
    %% Get D (diagonal weights in the regularization) for the current x
    Rx = params.R * xnew;
    Dinv = get_Dinv(Rx, params);
        
    %% Do CG for solving (Dinv + R * Gamma^-1 * R') z = R * Gamma^-1 * bi
    z = doCGMIL(Dinv, bi, z, params);
    
    bRz = bi - ( params.R )' * z;
    
    zk = bRz / params.mEAWA; % Update solution
    
    %% Update tk
    if( params.doMFISTA )
        cost = compute_Cost_CT2D(data, zk, params); % Compute Cost at current estimate
        dcost = cost - cost_old;
        if(dcost < 0)
            xk = zk;
        else
            xk = xold;
        end
        tkp1 = (1 + sqrt(1 + 4*tk^2))/2;
        yk = xk + tk/tkp1*(zk - xk) + (tk-1)/tkp1*(xk - xold); % Two step update
        tk = tkp1;
        xnew = xk;
    else
        xnew = zk;
        yk = xnew;
    end
    itr = itr + 1; %% Book Keeping

    time_elapse = toc;

    %% Compute Cost at current estimate
    cost_new = compute_Cost_CT2D(data, xnew, params); % Reinitialize cost_old for current sigma
    diffcost = cost_new - cost_old;
    
    %% Compute the norm difference between successive iterates
    diffx = sqrt(sum(abs(xnew(:)-xold(:)).^2));

    %% Compute RMSE at current estimate
    err = params.ig.mask.*( xnew - params.img );
    err = sqrt( sum( abs( err(:) ).^2 ) / Nmask );
    
    %% Show updates of cost, estimate, etc
    C(itr) = cost_new;
    T(itr) = time_elapse;
    E(itr) = sqrt( sum( params.ig.mask(:) .* abs( xnew(:) - params.xinf(:) ) .^ 2)) / params.xinfnorm;
    RMSE(itr) = err;

    if(dispfig)
        if(~mod(itr-1, dispitrnum))
            figure(figno); clf;
            subplot(2,2,1); imagesc(scale * params.img(zr,zc), clim); colormap(gray); colorbar; axis off; axis equal; title('Noisefree phantom, x');
            subplot(2,2,2); imagesc(scale * err(zr,zc), clim); colormap(gray); colorbar; axis off; axis equal; title(cat(2,'x - x', int2str(itr), '; ', num2str(RMSE(itr), fstrC)));
            subplot(2,2,3); imagesc(scale * xini(zr,zc), clim); colormap(gray); colorbar; axis off; axis equal; title('Initial estimate x0');
            subplot(2,2,4); imagesc(scale * xnew(zr,zc), clim); colormap(gray); colorbar; axis equal; axis off; title(cat(2, 'Current estimate x', int2str(itr)));
            pause(0.1);
        end
    end
    
    %% Display cost and other parameters
    strStatus{itr} = ['MFISTA: ' int2str(itr-1) '/' int2str(maxitr) ...
                     '; C=' num2str(cost_new, fstrC) ';DC=' num2str(diffcost, fstrC) ';Dx/||x||=', num2str(diffx/normx, fstrO) ...
                     '; l2D=' num2str(E(itr), fstrC) '; RMSE=' num2str(err, fstrC) '; ' num2str(time_elapse, fstrT) 's'];
    
     if(dispitr && ~mod(itr-1, dispitrnum))
         disp('--------------------------------------------------------------------');
         disp(strStatus{itr});
     end
end

if(dispitr && mod(itr-1, dispitrnum))
    disp(strStatus{itr});
    disp('--------------------------------------------------------------------');
end