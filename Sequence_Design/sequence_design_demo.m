clc
clear
%This script is used as a demo for gen_MRF_sequence_pulse.m. The function
%uses TR, TE and FA from Spiral_Design/Data/method_orig.mat to construct a
%sequence. 
% 8/2019 Enlin Qian
% # Copyright of the Board of Trustees of Columbia University in the City of New York

%% Sequence Design
addpath(genpath('.'));
[kshot, dcf, ind, TR_all, FA_all, TE_all] = gen_MRF_sequence_pulseq();