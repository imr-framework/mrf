function [om_store,echoes,seq] = EPGsim_TSE(alpha,N,esp,use_y90,rlx)
%[om_store,echoes,seq] = EPGsim_TSE(alphas,N,esp,use_y90,rlx)
% Simulates EPG of Turbo Spin Echo sequence with an initial 
%                 (90,90) pulse followed by N (0,alpha) pulses
%  alpha : flip angle of 2nd to (N+1)th pulses  
%  N : number of pulses after the initial (90,90) pulse
%  esp: TR between pulses (also )
%  use_y90: boolean indicating whether we use a (90,90) pulse at the
%           beginning (if false, all pulses are (0,alpha))
%  rlx: mode of relaxation

if nargin < 4
    use_y90 = 1;
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

dt = esp/2; % time evolves in 0.5*esp steps (dt = 0.5*esp -> dk = 1) 

seq.name = 'Turbo Spin Echo';
if use_y90 == 1
    seq.rf(:,1) = [90,90]';
    seq.rf(:,2:N) = repmat([0,alpha]',1,N-1);
else 
    seq.rf(:,1:N) = repmat([0,alpha]',1,N); 
end

seq.time = [0 dt dt];
seq.events = {'rf','grad','relax'};

for n = 1:N-1
    % Order of operators : T(rf)->S(grad)->E(relax) "TSE",easy to remember!
    seq.events{end+1} = 'rf';
    seq.events{end+1} = 'grad';
    seq.events{end+1} = 'relax';
    seq.events{end+1} = 'grad';
    seq.events{end+1} = 'relax';
    seq.time = [seq.time (2*n-1)*dt 2*n*dt 2*n*dt (2*n+1)*dt (2*n+1)*dt];
end
seq.grad = ones(1,2*N-1);
[om_store,echoes] = EPG_custom(seq);

end