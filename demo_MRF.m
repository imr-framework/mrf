%% This script functions as a wrapper for all components of MRF package
clc
clear

%% Design a sequence
addpath(genpath('.'));
[kshot, dcf, ind, TR_all, FA_all, TE_all] = gen_MRF_sequence_pulseq();

%% Reconstruct image
image_data_final_Complex = MRF_recon();

%% Dictionary simulation
phis = zeros(1,1001); % the sequence starts with a 180 pulse, thus 1001 points
alphas = [180, FA_all']; % all flip angles
TRs = [18, TR_all'.*1000]; % all TRs, in ms
T1_range = [];
T2_range = [];
%% Dictionary matching

%% ROI analysis tool