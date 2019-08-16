clc
clear
%This script is used as a demo for MRF_recon.m. The function takes in a
%.dat raw data file and .mat trajectory file and using sliding windown
%reconstruction and complex coil combination to reconstruct the image.
%These two files can both be found under Image_Reconstruction/Data
% 8/2019 Enlin Qian
% # Copyright of the Board of Trustees of Columbia University in the City of New York

addpath(genpath('.'));
image_data_final_Complex = MRF_recon();