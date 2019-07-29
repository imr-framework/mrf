clc
clear
%This script is used as a demo for NIST_ROI_ana.m. For now, four sets of 
%data are used for testing.
%       T1_map.mat  T1_map data
% T1_map_rot45.mat  T1_map data rotated counterclockwise 45 degrees
%       T2_map.mat  T2_map data
% T2_map_rot30.mat  T2_map data rotated counterclockwise 30 degree
% 7/2019 Enlin Qian
% # Copyright of the Board of Trustees of Columbia University in the City of New York

%% Load data
% Determine where your m-file's folder is.
folder = fileparts(which('T2_map_rot30.mat')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
load('T2_map_rot30.mat')
sphere_par = NIST_ROI_ana(T2_map_rot, 'T2', 128, 210);