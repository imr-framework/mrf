% This script is used as a demo for MRF_Dict_Matching.m. The function takes in a
% .mat dictionary file and .mat recontructed image file and using vector dot
% product to find the most similar dictionary entry. 
% The .mat dictionary file can be located in Dictionary_Matching/Data, while
% the reconstructed image data file is too large to be included in a Github
% repository. They can be found in Google drive. Shareable link are below.
% T1 Array slice: https://drive.google.com/open?id=1G4hsc00CS5ycbnQ6KQa8YDlkis_Yc-5f
% T2 Array slice: https://drive.google.com/open?id=199m40D-eiPWnAMtcNZj6gIM0aDZ-VZoB
%
% 8/2019 Enlin Qian
% # Copyright of the Board of Trustees of Columbia University in the City of New York

%% Image Reconstruction
addpath(genpath('.'));
load('MRF_dict.mat')
Par_map = MRF_Dict_Matching(MRF_dict, image_data_final_Complex);

%% Some testing (sphere 1)
plot(mat2gray(abs(squeeze(image_data_final_Complex(124,86,:)))))
hold on
plot(mat2gray(abs(MRF_dict.dict_SW_norm(3096,:))))
hold off
%% Some testing (sphere 2)
plot(mat2gray(abs(squeeze(image_data_final_Complex(158,96,:)))))
hold on
plot(mat2gray(abs(MRF_dict.dict_SW_norm(2840,:))))
hold off