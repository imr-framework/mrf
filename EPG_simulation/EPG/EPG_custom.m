function [om_store,echoes] = EPG_custom(seq)
%[om_store,echoes] = EPG_custom(seq)
% Performs EPG simulation of a general pulse sequence
% INPUTS 
%     seq: sequence struct with the required fields
%         seq.rf - 2 x P matrix with each column = [phi,alpha]'
%         seq.grad - vector of k-shifts (must be integers)
%         seq.events - cell of strings: 'rf','grad', or 'relax'
%         seq.time - vector of timing for each event in seq.events
%         seq.T1, seq.T2 : T1 and T2 values for relaxation 
%                         (set both to 0 for no relaxation)
% OUTPUTS
%
%      om_store - 1 x Q cell array of matrices with size (3 x K).
%      Each matrix is a record of configuration states (F+(k),F-(k),Z(k)
%      at all k's generated up to that point. Each matrix corresponds to
%      the elements in seq.events and seq.timing at the same index.
%
%      echoes - output of echo information using findEchoes 

%%  Inputs
rf = seq.rf;
grad = seq.grad;
timing = seq.time; % in ms
uniqtimes = unique(timing);
events = seq.events; % 3 types of events: 'rf','grad', and 'relax'
N = length(events);
T1 = seq.T1; % must be populated - set to 0 for no relaxation
T2 = seq.T2; % must be populated - set to 0 for no relaxation
%% Initialize
omega = [0 0 1].'; % initial magnetization is always at equilibrium (+z)
%delk=1; %Unit dephasing
rf_index = 1;
om_index = 1;
grad_index = 1;
%% Describe pulse sequence in steps
om_store = cell(1,N);
for n = 1:N
   switch events{n}
       case 'rf'
           omega = rf_rotation(rf(1,rf_index),rf(2,rf_index))*omega;
           rf_index = rf_index + 1;
       case 'grad'
           omega = shift_grad(grad(grad_index),omega);     
           grad_index = grad_index + 1;
       case 'relax'
           q = find(uniqtimes==timing(n));
           tau = uniqtimes(q) - uniqtimes(q-1);
           omega = relax(tau,T1,T2,omega);
   end
   om_store{om_index} = omega;
   om_index = om_index + 1;
end
% Find and store echos
echoes = findEchoes(seq,om_store);

end

