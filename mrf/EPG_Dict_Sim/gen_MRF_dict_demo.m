clc
clear 

%This script is used as a demo for dict_sim_MRF.m. All TR, FA, and TE are
%loaded based on the sequence generated in Sequence Design part. 
% 7/2019 Enlin Qian
% # Copyright of the Board of Trustees of Columbia University in the City of New York
addpath(genpath('.'));
load('TR_TE_FA.mat')
%% Dictionary Simulation
%   Modified based on ISMRM/NIST phantom specs - T1 array
T1 = 0.02:0.005:0.1;
T1 = [T1 0.1:0.01:0.5];
T1 = [T1 0.5:0.2:2];
% Modified based on ISMRM/NIST phantom specs - T2 array
T2 = 0.005:0.002:0.02;
T2 = [T2 0.02:0.01:0.1];
T2 = [T2 0.1:0.025:0.3];
T2 = [T2 0.3:0.1:1];
TI = 0.018;
kshots=48;
MRF_dict = gen_MRF_dict(TR_all, TE_all, FA_all, T1, T2, TI, kshots);

