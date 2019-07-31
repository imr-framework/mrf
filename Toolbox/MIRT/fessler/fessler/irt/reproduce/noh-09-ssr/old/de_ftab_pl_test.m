%Test for PL for DE-CT
clear all;load de_ct_setup.mat; %saved data upto here
[shat cost_eval] = de_ftab_s_pl(ftab.fit,ftab.xray,ftab.mac.mac,fhat.raw,ymi,...
                                'niter',300,'beta',2^(-7),'curvtype',1); %by PL

shat = reshape(shat,[n.b n.a 2]);
imagesc(shat(:,:,1)');colorbar;colormap gray;
