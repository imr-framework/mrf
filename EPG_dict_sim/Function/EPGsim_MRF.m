function [om_store,echoes,seq] = EPGsim_MRF(phis,alphas,TRs,rlx)
% [om_store, echoes, seq] = EPGsim_MRF(phis,alphas,TRs,rlx)
% Simulates variable flip angle & variable timing EPG
% Author: Gehua Tong 
% Date: Feb 25, 2019

% Note: 
%    (1) phis and alphas must have the same size & be in degrees
%    (2) TRs are in milliseconds; TRs must be the same size as phis &
%    alphas
%    (3) rlx = [T1 T2] in milliseconds(set both to zero for no relaxation) 
% Important: all gradients = 0, so we are always at coherence order 0
%            (i.e. only RF flipping and relaxation happen)

N = length(TRs);
seq.name = 'Variable FA & TR for MRF';
seq.rf = [phis(:)';alphas(:)'];
diffTime = [zeros(1,N);repmat(TRs,2,1)];
seq.time = repelem(cumsum([0 TRs(1:end-1)]),1,3) + diffTime(:)';
seq.events = repmat({'rf','grad','relax'},1,length(TRs));
seq.grad = zeros(1,length(TRs));
seq.T1 = rlx(1); seq.T2 = rlx(2); 

[om_store,echoes] = EPG_custom(seq);


end