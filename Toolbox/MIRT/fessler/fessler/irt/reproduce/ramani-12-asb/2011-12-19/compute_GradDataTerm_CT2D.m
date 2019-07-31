function [DT Ax] = compute_GradDataTerm_CT2D(x, params)
% 2011-11-30, Sathish Ramani, University of Michigan
% Do the matrix product A' * W * A

%% Data Term
Ax = params.A * x;
DT = ( params.A )' * ( params.W .* Ax );