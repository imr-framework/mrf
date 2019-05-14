function [om_store,echoes,seq] = EPGsim_HE(rf,esp,rlx)
%[om_store,echoes,seq] = EPGsim_HE(RFseq,esp,rlx)
% Generates EPG for the hyperecho pulse sequence, a sensitive test for
% the accuracy of the EPG framework. With an initial 90 pulse and 
% a midway 180 pulse, the hyperecho sequence mirror the phase, flip
% angle, and timing of the first half in the second half to generate a
% perfect echo (if relaxation is neglected) at the end of the sequence.
%
% 
% INPUTS
%       rf - 2 x N matrix representing N arbitrary RF pulses
%            which are reversed for the second half of the sequence
%       esp - spacing between RF pulses
%       rlx - relaxation ([T1,T2]; set to [0,0] for no relxation)
% alphas : length m array of flip angles.
% TR: spacing between RF pulses
% rlx: mode of relaxation. Default : no relaxation
%
% Gehua Tong, Nov 19 2018


if nargin < 3
    rlx = 'none';
end
if ischar(rlx)
    switch rlx
        case 'none'
             seq.T1 = 0; seq.T2 = 0;% zero T1,T2 means no relaxation in relax.m 
    
        case 'gm' % gray matter
             seq.T1 = 1300; seq.T2 = 110;
       
        case 'wm' % white matter
             seq.T1 = 960; seq.T2 = 80; 
     
        case 'csf'
             seq.T1 = 3600; seq.T2 = 1800;
        
        case 'default'
             seq.T1 = 1000; seq.T2 = 100; % same as in EPG paper for reference
    end 
elseif isequal(size(rlx),[1,2])
       seq.T1 = rlx(1); seq.T2 = rlx(2);
else
    seq.T1 = 0; seq.T2 = 0; % zero T1,T2 means no relaxation in relax.m 
end

seq.name = 'Hyperecho test';

dt = esp/2;
m = size(rf,2);
N = 4*m + 2; % number of time intervals
seq.rf = [[90 90]',rf,[0 180]',fliplr(-rf)];
seq.time = [0 dt];
seq.events = {'rf','grad'};

for n = 1:2*m
    % Order of operators : T(rf)->S(grad)->E(relax) "TSE",easy to remember!
    seq.events{end+1} = 'rf';
    seq.events{end+1} = 'grad';
    seq.events{end+1} = 'grad';
    ta = seq.time(end);
    seq.time = [seq.time ta ta+dt ta+2*dt];
end
seq.time = [seq.time seq.time(end) seq.time(end)+dt];
seq.events{end+1} = 'rf';
seq.events{end+1} = 'grad';
seq.grad = ones(1,N);

[om_store,echoes] = EPG_custom(seq);


end
