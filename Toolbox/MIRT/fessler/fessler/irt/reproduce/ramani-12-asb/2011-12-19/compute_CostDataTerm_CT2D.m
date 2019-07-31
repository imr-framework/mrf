function DT = compute_CostDataTerm_CT2D(data, x, params)
% 2011-11-30, Sathish Ramani, University of Michigan
% Compute the data part of the cost function
% 1/2 * (y-Ax)^H * W * (y-Ax);

Ax = params.A * x;
Ax = data - Ax;

DT = sum( params.W(:) .* abs( Ax(:) ) .^ 2 ) / 2;