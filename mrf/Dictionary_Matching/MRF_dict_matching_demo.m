% This script is used as a demo for MRF_Dict_Matching.m. The function takes in a
% .mat dictionary file and .mat recontructed image file and using vector dot
% product to find the most similar dictionary entry. 
% The .mat dictionary file can be located in Dictionary_Matching/Data, while
% the reconstructed image data file is too large to be included in a Github
% repository. They can be found in Google drive. Shareable link are below.
% T1 Array slice: https://drive.google.com/open?id=1G4hsc00CS5ycbnQ6KQa8YDlkis_Yc-5f
% G2 Array slice: https://drive.google.com/open?id=199m40D-eiPWnAMtcNZj6gIM0aDZ-VZoB
%
% 8/2019 Enlin Qian
% # Copyright of the Board of Trustees of Columbia University in the City of New York

%% Image Reconstruction
addpath(genpath('.'));
load('dict_example.mat')
Par_map = MRF_Dict_Matching(MRF_dict, image_data_final_Complex);