function measval = computemeasure(lamcurr, ptemp)
% 2011-11-30, Sathish Ramani, University of Michigan
CAARR = ptemp.CAA + lamcurr * ptemp.RR; % Eig.vals of the circ.matrix CAA + nu * R'R

mxE = max(CAARR(:));
mnE = min(CAARR(:));

measval = mxE/mnE; % Condition number