function [xnew CAL TAL EAL RMSE strStatus] = runADMM(data, xnew, params)
% ADMM algorithm in 
% Ramani et al., "A splitting based iterative algorithm 
% for statistical X-ray CT Reconstruction,"
% IEEE TMI, in press
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

mu = params.AL.mu;
nu1 = params.AL.nu1;

zr = params.zoomr;
zc = params.zoomc;

fstrC = params.formatstringC;
fstrT = params.formatstringT;
fstrO = params.formatstringO;

%% Compute Error in estimate
err = params.ig.mask .* ( xnew - params.img );
errnorm = sqrt( sum( abs( err(:) ).^2 ) / Nmask );

time_elapse = 0;
itr = 1;

%% Initialization
Ax = params.A * xnew;
Rx = params.R * xnew;

%% Prepare cost and thresholds
cost_new = compute_Cost_CT2D(data, xnew, params, Ax, Rx); % Reinitialize cost_old for current sigma
cost_old = 2*cost_new;
diffcost = cost_new - cost_old;

normx = sqrt( sum( abs( xnew(:) ) .^ 2 ) );
diffx = normx;

CAL(itr) = cost_new;
TAL(itr) = time_elapse;
EAL(itr) = sqrt( sum( params.ig.mask(:) .* abs( xnew(:) - params.xinf(:) ) .^ 2)) / params.xinfnorm;
RMSE(itr) = errnorm;

%% Show Images
if(dispfig)
    figure(figno); clf;
    subplot(2,2,1); imagesc(scale * params.img(zr,zc), clim); colormap(gray); colorbar; axis off; axis equal; title('Noisefree phantom, x');
    subplot(2,2,2); imagesc(scale * err(zr,zc)); colormap(gray); colorbar; axis off; axis equal; title(cat(2,'x - x', int2str(itr), '; ', num2str(RMSE(itr), fstrC)));
    subplot(2,2,3); imagesc(scale * xini(zr,zc), clim); colormap(gray); colorbar; axis off; axis equal; title('Initial estimate x0');
    subplot(2,2,4); imagesc(scale * xnew(zr,zc), clim); colormap(gray); colorbar; axis off; axis equal; title(cat(2, 'Current estimate x', int2str(itr)));
end

%% Display cost and other parameters
strStatus{itr} = ['ADMM: ' int2str(itr-1) '/' int2str(maxitr) ...
                  '; C=' num2str(cost_new, fstrC) '; DC=' num2str(diffcost, fstrC) '; Dx/||x||=' num2str(diffx/normx, fstrO) ...
                  '; l2D=' num2str(EAL(itr), fstrC) '; RMSE=' num2str(errnorm, fstrC) '; ' num2str(time_elapse, fstrT) 's'];
if(dispitr)
    disp('--------------------------------------------------------------------');
    disp(strStatus{itr});
end

%% ADMM
d = cell(1,2); d{1} = zeros(size(Ax)); d{2} = zeros(size(Rx));
z = d;
z{1} = params.Wy + mu * (Ax + d{1});
z{1} = params.AL.iWmu .* z{1};
z{2} = Rx + d{2};
z{2} = thresholdRx(z{2}, params, nu1*mu); % Perform thresholding

while((itr <= maxitr) && (diffx/normx >= dxtol) && (diffcost <= dcosttol)) % 0)) %
    tic;
    
    cost_old = cost_new;
    xold = xnew;
    
    % Compute x = (A'A + nu1*R'R)^(-1) (A'(u0-d0) + nu1*R'(u1-d1))
    Ru1d1 = ( params.R )' * (z{2} - d{2}); % R'(u1-d1)
    rhsx = ( params.A )' * (z{1} - d{1}) + nu1 * Ru1d1;
    xnew = doCGAL_AARR(rhsx, xnew, params);
    
    % Book Keeping
    Ax = params.A * xnew;
    Rx = params.R * xnew;
    
    % Compute u0 = (W + mu*I)^-1 (W*y + mu*(Ax + d0))
    rhsu0 = params.Wy + mu * (Ax + d{1});
    z{1} = params.AL.iWmu .* rhsu0;
    
    % Compute u1 = arg min_{u1} lambda/mu* Psi(u1) + 0.5*||u1 - Rx - d1||^2
    rhsu1 = Rx + d{2};
    z{2} = thresholdRx(rhsu1, params, nu1*mu); % Perform thresholding

    % Update d
    d{1} = d{1} - (z{1} - Ax);
    d{2} = d{2} - (z{2} - Rx);
    itr = itr + 1;
    time_elapse = toc;

    %% Compute Cost at current estimate
    cost_new = compute_Cost_CT2D(data, xnew, params, Ax, Rx); % Reinitialize cost_old for current sigma
    diffcost = cost_new - cost_old;

    %% Compute error
    err = params.ig.mask .* ( xnew - params.img );
    errnorm = sqrt( sum( abs( err(:) ).^2 ) / Nmask );
    
    CAL(itr) = cost_new;
    TAL(itr) = time_elapse;
    EAL(itr) = sqrt(sum(params.ig.mask(:).*abs(xnew(:) - params.xinf(:)).^2)) / params.xinfnorm;
    RMSE(itr) = errnorm;

    %% Compute the norm difference between successive iterates
    diffx = sqrt(sum(abs(xnew(:)-xold(:)).^2));

    %% Show updates of cost, estimate, etc
    if(dispfig)
        if(~mod(itr-1, dispitrnum))
            figure(figno); clf;
    subplot(2,2,1); imagesc(scale * params.img(zr,zc), clim); colormap(gray); colorbar; axis off; axis equal; title('Noisefree phantom, x');
            subplot(2,2,2); imagesc(scale * err(zr,zc)); colormap(gray); colorbar; axis off; axis equal; title(cat(2,'x - x', int2str(itr), '; ', num2str(RMSE(itr), fstrC)));
            subplot(2,2,3); imagesc(scale * xini(zr,zc), clim); colormap(gray); colorbar; axis off; axis equal; title('Initial estimate x0');
            subplot(2,2,4); imagesc(scale * xnew(zr,zc), clim); colormap(gray); colorbar; axis equal; axis off; title(cat(2, 'Current estimate x', int2str(itr)));
            pause(0.1);
        end
    end
    
    %% Display cost and other parameters
    strStatus{itr} = ['ADMM: ' int2str(itr-1) '/' int2str(maxitr) ...
                      '; C=' num2str(cost_new, fstrC) ';DC=' num2str(diffcost, fstrC) ';Dx/||x||=' num2str(diffx/normx, fstrO) ...
                      '; l2D=' num2str(EAL(itr), fstrC) '; RMSE=' num2str(errnorm, fstrC) '; ' num2str(time_elapse, fstrT) 's' ...
                      ];
    if(dispitr && ~mod(itr-1, dispitrnum))
        disp('--------------------------------------------------------------------');
        disp(strStatus{itr});
    end
end

if(dispitr && mod(itr-1, dispitrnum))
    disp('--------------------------------------------------------------------');
    disp(strStatus{itr});
end