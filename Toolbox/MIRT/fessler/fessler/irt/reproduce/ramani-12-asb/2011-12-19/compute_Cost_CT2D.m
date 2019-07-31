function [cost costDT costRT] = compute_Cost_CT2D(data, x, params, varargin)
% 2011-11-30, Sathish Ramani, University of Michigan
% Compute the cost function sum_i ||yi - TSi x||^2 + lambda * sum(rho(|grad(x)|))
% varargin{1} has Ax
% varargin{1} has Rx
lvarargin = length(varargin);

if( lvarargin == 0 )
    doDT = 1;
    doRT = 1;
elseif( lvarargin == 1 )
    doDT = 0;
    doRT = 1;
else
    doDT = 0;
    doRT = 0;
end

x = x .* params.ig.mask;

%% Compute the Data Term
if( doDT )
    costDT = compute_CostDataTerm_CT2D(data, x, params);
else
    Ax = data - varargin{1};
    costDT = sum( params.W(:) .* abs( Ax(:) ) .^ 2 ) / 2;
end

%% Compute Regularization term
if( doRT )
    Rx = params.R * x;
else
    Rx = varargin{2};
end

switch params.Operator
    case{'FD'} % gradient-norm of finite difference
        D = repmat( sqrt( sum( abs( Rx ) ) .^ 2, 3 ), [ 1 1 size( Rx, 3 ) + 1 ] );
        
    case{'AFD', 'W'} % Bilateral finite-differences or shift-invariant wavelets
        D = abs( Rx );
        
    otherwise
        error('Only available options are FD, AFD, W; others will be included in developments');
end

switch ( params.PriorType )
    case{'l1'}
        D = D .* repmat( params.rw, [ 1 1 size( Rx, 3 ) ] );
        
    case{'FP'}
        D = ( D / params.Prior.alpha - log( 1 + D / params.Prior.alpha ) )  .* repmat( params.rw, [1 1 size( Rx, 3 ) ] ) * params.Prior.alpha ^ 2;
        
    otherwise
        error('Available options are l1 and FP; Other potential functions will be included in developments');
end
costRT = params.lambda * sum( D(:) );
cost = costDT + costRT;