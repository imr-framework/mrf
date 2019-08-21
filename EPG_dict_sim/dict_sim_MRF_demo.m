clc
clear 

%This script is used as a demo for dict_sim_MRF.m. All TR, FA, and TE are
%loaded based on the sequence generated in Sequence Design part. 
% 7/2019 Enlin Qian
% # Copyright of the Board of Trustees of Columbia University in the City of New York
addpath(genpath('.'));
load('TR_TE_FA.mat')
%% Dictionary Simulation
phis = zeros(1,1001); % the sequence starts with a 180 pulse, thus 1001 points
alphas = [180, FA_all']; % all flip angles
TRs = [18, TR_all'.*1000]; % all TRs, in ms
T1_range=0.02:0.005:0.065;
T1_range=[T1_range, 0.09:0.03:0.36];
T1_range=[T1_range, 0.5:0.1:4];
T2_range=0.005:0.002:0.011;
T2_range=[T2_range, 0.015:0.005:0.045];
T2_range=[T2_range, 0.065:0.03:0.275];
T2_range=[T2_range, 0.4:0.1:2];
n3 = 1;
Echo_Final=zeros(length(T1_range)*length(T2_range), 1000);
for n1 = 1:length(T1_range)
    for n2 = 1:length(T2_range)
        disp (n3)
        [om_store,echoes,seq] = EPGsim_MRF(phis,alphas,TRs,[T1_range(n1), T2_range(n2)]);
        Echo_Final(n3, :) = Var_TE(om_store, TE_all, T2_range(n2));
        n3=n3+1;
    end
end
