%% This script functions as a wrapper for all components of MRF package
clc
clear

%% Sequence Design
addpath(genpath('.'));
[kshot, dcf, ind, TR_all, FA_all, TE_all] = gen_MRF_sequence_pulseq();

%% Image Reconstruction
image_data_final_Complex = MRF_recon();

%% Dictionary Simulation
T1=0.02:0.005:0.065;
T1=[T1, 0.09:0.03:0.36];
T1=[T1, 0.5:0.1:4];
T2=0.005:0.002:0.011;
T2=[T2, 0.015:0.005:0.045];
T2=[T2, 0.065:0.03:0.245];
T2=[T2, 0.35:0.05:0.5];
T2=[T2, 0.7:0.1:2];
TI=18;
kshots=48;
MRF_dict = gen_MRF_dict(TR_all, TE_all, FA_all, T1, T2, TI, kshots);

%% Dictionary Matching
load('dict_example.mat')
Par_map = MRF_Dict_Matching(MRF_dict, image_data_final_Complex);

%% ROI Analysis Tool