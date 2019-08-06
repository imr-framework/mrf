%% This script functions as a wrapper for all components of MRF package
clc
clear

%% Design a sequence
addpath(genpath('.'));
% addpath(genpath('.\Spiral_Design\')) %add Spiral_Design and subfolders
% addpath(genpath('.\')) %add Spiral_Design path and sublfolders
gen_MRF_sequence_pulseq %generate MRF sequence

%% Reconstruct image
