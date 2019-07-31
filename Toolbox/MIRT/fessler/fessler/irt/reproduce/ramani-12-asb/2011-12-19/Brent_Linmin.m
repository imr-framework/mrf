function [minmeas_temp lamopt_temp] = Brent_Linmin(ptemp)
% 2011-11-30, Sathish Ramani, University of Michigan
CAA = ptemp.CAA;
np = real( max( CAA(:) ) / min( CAA(:) ) );

% Bracket the minimum
disp('Bracketing the optimal Step-size');
abcf = Bracket_Minimum(ptemp);
disp(['             Bracket = ' num2str(abcf(1)) '; ' num2str(abcf(2)) '; ' num2str(abcf(3))]);

% Brent 1D line search
disp('             Finding the optimal Step-size using Brent 1D search within the bracket');
[minmeas_temp lamopt_temp] = Brent1D(abcf, ptemp);

% Check if that is the real minimum
if(minmeas_temp>=np)
    lamopt_temp = 0;
    minmeas_temp = np;
end