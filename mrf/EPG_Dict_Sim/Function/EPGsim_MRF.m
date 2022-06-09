function [Echo_Final,echoes,seq] = EPGsim_MRF(phis,alphas,TRs,TEs,T1,T2,TI)
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
seq.grad = ones(1,length(TRs));  % add gradient
seq.T1 = T1; seq.T2 = T2; 

%% Check if 180 pulse needs to be added
if TI ~= 0  % if TI is not 0, add 180 pulse to Flip angles
    seq.rf = ([0,180;seq.rf'])';
    seq.events = [{"rf","grad","relax"}, seq.events];
    seq.time = [0,TI,TI,seq.time+TI];
    seq.grad = [0,seq.grad];
end


[om_store,echoes] = EPG_custom(seq);
% Scale with T2 relaxation
F_plus = cellfun(@(v)v(1), om_store);
F_plus = F_plus(1:3:end);
Echo_Final = F_plus(1,2:1:end).*exp(-TEs./T2);

end