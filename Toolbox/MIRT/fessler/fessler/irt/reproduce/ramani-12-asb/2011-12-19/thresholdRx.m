function v = thresholdRx(Rx, params, tval)
% 2011-11-30, Sathish Ramani, University of Michigan
%% Prepare coefficients
switch(params.Operator)
    case{'FD'} % gradient-norm of finite difference
        Dtemp = repmat( sqrt( sum( abs( Rx ) .^ 2, 3 ) ), [1 1 size(Rx, 3) ] );
        
    case{'AFD', 'W'} % Anisotropic TV or shift-invariant wavelets
        Dtemp = abs( Rx );
        
    otherwise
        error('Only available options are FD, AFD, W; others will be included in developments');
end

%% Perform coordinate-wise thresholding
lammu = params.lambda / tval * repmat( params.rw, [1 1 size(Rx, 3)] ); % Space varying weights;

switch params.PriorType
    case{'l1'} % l1 norm of the sparsifying transform coefficients
        D = Dtemp - lammu; D = ( D + abs( D ) ) / 2; % Max(x, 0) operation
        v = Rx ./ Dtemp .* D;
        
    case{'FP'} % Fair Potential
        kappa =  Dtemp - lammu * params.Prior.alpha - params.Prior.alpha;
        ckappa = 0.5 * ( kappa + sqrt( kappa .^ 2 + 4 * params.Prior.alpha * Dtemp ) );
        calpha = 1 + lammu * params.Prior.alpha ./ ( params.Prior.alpha + ckappa );
        v = Rx ./ calpha;
        
    otherwise
        error('Available options are l1 and FP; Other potential functions will be included in developments');
end

%% Ensure v does not contain NaN obtained from 0/0 operations
v( isnan( v ) ) = 0;
