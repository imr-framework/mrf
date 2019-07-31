function [Dinv Dtemp] = get_Dinv(Rx, params)
% 2011-11-30, Sathish Ramani, University of Michigan
%% Computer the Diagonal Weighting Term
switch params.Operator
    case{'FD'} % gradient-norm of finite difference
        Dinv = repmat( sqrt( sum( abs( Rx ) .^ 2, 3 ) ), [1 1 size( Rx, 3 ) ] );
        
    case{'AFD', 'W'} % Bilateral finite-differences or shift-invariant wavelets
        Dinv = abs( Rx );
        
    otherwise
        error('Only available options are FD, AFD, W; others will be included in developments');
end

%% Compute diagonal weighting matrix
switch params.PriorType
    case{'l1'} % l1 norm of the sparsifying transform coefficients
        Dinv = Dinv ./ ( params.lambda * repmat( params.rw, [1 1 size( Rx, 3 ) ] ) ); % Originally kappas (and the rw's) are designed for bilateral kinda reg-operatos only
        
    case{'FP'} % Fair potential
        Dinv = ( Dinv + params.Prior.alpha ) ./ ( params.lambda * repmat( params.rw, [1 1 size( Rx, 3 ) ] ) * params.Prior.alpha );
        
    otherwise
        error('Available options are l1 and FP; Other potential functions will be included in developments');
end
