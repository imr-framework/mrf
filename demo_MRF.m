%% This script functions as a wrapper for all components of MRF package
clc
clear

%% Sequence Design
addpath(genpath('.'));
[kshot, dcf, ind, TR_all, FA_all, TE_all] = gen_MRF_sequence_pulseq();

%% Image Reconstruction
image_data_final_Complex = MRF_recon();

%% Dictionary Simulation
phis = zeros(1,1001); % the sequence starts with a 180 pulse, thus 1001 points
alphas = [180, FA_all']; % all flip angles
TRs = [18, TR_all'.*1000]; % all TRs, in ms
T1_range = 1:1:4000;
T2_range = 1:1:2000;
[om_store_WM,echoes,seq] = EPGsim_MRF(phis,alphas,TRs,rlxWM);
Echo_Final_WM = Var_TE(om_store_WM, TE_all, rlxWM(2));
%% Dictionary Matching

%% ROI Analysis Tool