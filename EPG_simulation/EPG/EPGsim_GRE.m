function [om_store,echoes,seq] = EPGsim_GRE(alpha,N,TR,rlx)
%[om_store,echoes,seq] = EPGsim_GRE(alpha,N,TR,rlx)
% EPG simulation of simple GRE 
% 
% The same RF phase (=0) & flip angle is repeated at each TR
%     starting from t = 0
% alpha : flip angle
% N : number of RF pulses
% TR : interval between neighboring RF pulses
% rlx: relaxation mode. Default is no relaxation.
%
% Gehua Tong, Nov 19 2018

if nargin < 4
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


seq.name = 'Simple gradient-recalled echo';
seq.rf(:,1:N) = repmat([0,alpha]',1,N); 
seq.time = [];
rftimes = 0:TR:TR*(N-1);
seq.events = {};
seq.grad = ones(1,N+1);

for n = 1:N
    % Order of operators : T(rf)->S(grad)->E(relax) "TSE",easy to remember!
    seq.events{end+1} = 'rf';
    seq.events{end+1} = 'grad';
    seq.events{end+1} = 'relax';
    seq.time = [seq.time,rftimes(n),rftimes(n)+TR,rftimes(n)+TR];
        
end


[om_store,echoes] = EPG_custom(seq);

end

