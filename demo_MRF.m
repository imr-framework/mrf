%% This script functions as a wrapper for all components of MRF package
clc
clear

%% Design a sequence
addpath(genpath('.'));
[kshot, dcf, ind, TR_all, FA_all, TE_all] = gen_MRF_sequence_pulseq();

%% Reconstruct image
image_data_final_Complex = MRF_recon();

%% Dictionary simulation

%% Dictionary matching

%% ROI analysis tool