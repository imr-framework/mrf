clc
clear 

%This script is used as a demo for dict_sim_MRF.m. All TR, FA, and TE are
%loaded based on the sequence generated in Sequence Design part. 
% 7/2019 Enlin Qian
% # Copyright of the Board of Trustees of Columbia University in the City of New York
addpath(genpath('.'));
load('TR_TE_FA.mat')
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

%% Plot on three spheres 
figure(1)
plot(abs(squeeze(MRF_dict.dict_SW_norm(1716,:))),'r','LineWidth',2)
hold on
plot(abs(squeeze(MRF_dict.dict_SW_norm(1029, :))),'b','LineWidth',2)
plot(abs(squeeze(MRF_dict.dict_SW_norm(627, :))),'k','LineWidth',2)
hold off
legend('Sphere 1', 'Sphere 5', 'Sphere 10');
