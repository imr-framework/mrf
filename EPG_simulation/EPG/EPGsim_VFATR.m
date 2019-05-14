function [om_store,echoes,seq] = EPGsim_VFATR(phis,alphas,TRs,dt,rlx)
% [om_store, echoes, seq] = EPGsim_VFATR(phis,alphas,TRs,dt,rlx)
% Simulates variable flip angle & variable timing EPG
% Gehua Tong; Feb 12 2019

% Note: 
%    (1) phis and alphas must have the same size & be in degrees
%    (2) TRs are in milliseconds; TRs must be the same size as phis &
%    alphas
%    (3) dt - smallest time interval that TRs are discretized into (ms)
%    (4) rlx = [T1 T2] in milliseconds(set both to zero for no relaxation) 
% Display of the EPG is not recommended for small discretizations.
% Instead, just plot the echoes (col 1: timing; col 2: intensity)


T1 = rlx(1); T2 = rlx(2);

seq.name = 'Variable FA & TR';
seq.rf = [phis(:)';alphas(:)'];

% Discretize TR
TRs_n = round(TRs/dt);
seq.time = [];
seq.events = {};
seq.grad = ones(1,sum(TRs_n));
for q = 1:length(TRs(:))
    if q == 1
        seq.time = [0];
    else
        seq.time = [seq.time, seq.time(end)];
    end
    seq.events{end+1} = 'rf';
    for k = 1:TRs_n(q)
        seq.time = [seq.time, seq.time(end)+dt, seq.time(end)+dt];
        seq.events{end+1} = 'grad';
        seq.events{end+1} = 'relax';
    end
end
seq.T1 = T1; seq.T2 = T2; 

[om_store,echoes] = EPG_custom(seq);

end