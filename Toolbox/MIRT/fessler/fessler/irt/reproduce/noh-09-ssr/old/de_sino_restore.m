%This is the code for DE-CT sinogram restoration
%Modified from spie02.m

clear all;
im off;
%%%%%%%%%%%%
%Load data %
%%%%%%%%%%%%
load xtrue.mat; %128 by 104
%load xtrue1024.mat; %1024 by 1024
[n.x n.y] = size(xtrue(:,:,1));
s_factor = 1;
n.b = 140*s_factor; n.a = 128*s_factor;f.dx = 0.16 * 2;
f.stype = 'ps1'; %setting dual spectra
ftab = de_ftab_build({f.stype});
%mask
mask = sum(xtrue, 3) > 0;
mask = imdilate(mask, strel('disk', 5));
mask = logical(mask);
%true density map display
im clf
c.soft = [0 1.2]; %max color range
c.bone = [0 2.2];
c.dens = [0 2.2];
im(221, xtrue(:,:,1), 'SoftTissue Density', c.soft), cbar([0 1])
im(222, xtrue(:,:,2), 'Bone Density', c.bone), cbar([0 2])
im(223, sum(xtrue,3), 'Density Map', c.dens), cbar([0 2])
im(224, (sum(xtrue,3)>0.5) + 20*double(mask), 'Reconstruction Support')

%%%%%%%%%%%%%%%%
%System matrix %
%%%%%%%%%%%%%%%%
ig = image_geom('nx', n.x, 'ny', n.y, 'dx', f.dx, 'mask', mask);
sg = sino_geom('par', 'nb', n.b, 'na', n.a, 'dr', f.dx);
G = Gtomo2_strip(sg, ig); %small size CT
%G = Gtomo2_dscmex(sg, ig); %full size CT
gi = reshape(sum(G'), [n.b n.a]); %simple projection
im clf, im(121, gi, 'gi (simple projection)'), im(122, gi>0, 'gi>0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ideal and Noisy Measurement % y_{mi}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Noiseless sinogram and f
%
tmp = reshape(xtrue, n.x*n.y, 2);
tic %noiseless forward projection
strue = reshape(G * tmp(mask,:), [n.b n.a 2]);
printm('forward projection time %0.3g', toc)

%sinogram display
im clf, im(421, strue(:,:,1), 'noiseless s1'), cbar, im(422, strue(:,:,2), 'noiseless s2'), cbar
s1max = max(col(strue(:,:,1)));
s2max = max(col(strue(:,:,2)));

ftrue = ftab.fm_fun(ftab, {strue(:,:,1), strue(:,:,2)}); %noiseless f eval (fn call)
im(423, ftrue(:,:,1), 'noiseless f1'), cbar,im(424, ftrue(:,:,2), 'noiseless f2'), cbar

%
%Noisy Measurement
%
ybi = zeros(n.b, n.a, 2); %bar{y}_{mi}
ybi(:,:,1) = ftab.xray.I(1) * exp(-ftrue(:,:,1)); %true mean
ybi(:,:,2) = ftab.xray.I(2) * exp(-ftrue(:,:,2));
im(425, ybi(:,:,1), 'ybar1'), cbar
im(426, ybi(:,:,2), 'ybar2'), cbar
%f.scale = inf; %for noiseless profile
%f.scale = 10^6 / ybi(1,1,2); %high dose
f.scale = 5*10^4 / ybi(1,1,2); %low dose
if isinf(f.scale)
	ymi = ybi;	% noiseless
else
	ymi = poisson(f.scale*ybi, 0) / f.scale; %corrupted by poisson
end %measurements look like sinogram

f.title1 = sprintf('%2.0f kVp', ftab.xray.kvp(1));
f.title2 = sprintf('%2.0f kVp', ftab.xray.kvp(2));
if isinf(f.scale)
	f.title1 = [f.title1 ' (Noiseless)'];
	f.title2 = [f.title2 ' (Noiseless)'];
else
	printm('total counts %0.3g', sum(ymi(:)) * f.scale)
end % total counts 5.96034e+009

im clf, 
t.y1 = ymi(:,:,1) * f.scale;	% show counts!
t.y2 = ymi(:,:,2) * f.scale;
c.ymi1 = [0 floor(max(t.y1(:))/1e3)*1e3];
c.ymi2 = [0 floor(max(t.y2(:))/1e4)*1e4];
im(211, t.y1, f.title1, c.ymi1), cbar(c.ymi1)
im(212, t.y2, f.title2, c.ymi2), cbar(c.ymi2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pre-filtering before FBP reconstruction %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isinf(f.scale) | 1
	f.kernel = [1]';
else
	f.kernel = [1 8 1]';
end
f.kernel = f.kernel / sum(f.kernel);	% radial smoothing
ymi_filt = convn(ymi, f.kernel, 'same');

% borrow neighbor(s) for any log(0) values
off = 1;
while any(ymi_filt(:) == 0)
	ii = find(ymi_filt(:) == 0);
	printm('fixing %d zeros in ymi_filt with %d', length(ii), off)
	ymi_filt(ii) = max(ymi_filt(ii+off), ymi_filt(ii-off));
	off = off + 1;
end, clear off
%
% Estimation of f ; fhat = -log(ymi / Imi) in eq(19)
%
fhat.raw(:,:,1) = -log(ymi_filt(:,:,1) / ftab.xray.I(1)); 
fhat.raw(:,:,2) = -log(ymi_filt(:,:,2) / ftab.xray.I(2));
fhat.raw(isinf(fhat.raw)) = 0;

im clf;c.fhat = [0 ceil(max(fhat.raw(:)))];
im(221, fhat.raw(:,:,1), f.title1, c.fhat), cbar
im(222, fhat.raw(:,:,2), f.title2, c.fhat), cbar

%Error plot between fhat and ftrue
plot(fhat.raw(:,:,1)-ftrue(:,:,1), fhat.raw(:,:,2)-ftrue(:,:,2), 'y.')
axis([-1 1 -0.5 0.5]);title('Error plot between fhat and ftrue');
ytick, xlabel 'f_1 error', ylabel 'f_2 error'
fhat.err.raw(1) = max_percent_diff(fhat.raw(:,:,1), ftrue(:,:,1));
fhat.err.raw(2) = max_percent_diff(fhat.raw(:,:,2), ftrue(:,:,2));
printm('fhat.raw err %0.3g%% %0.3g%%', fhat.err.raw(1), fhat.err.raw(2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DECT sinogram restoration %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all;load de_ct_setup.mat; %high dose
%clear all;load de_ct_low4_setup.mat; %low dose : f.scale = 5*10^4 / ybi(1,1,2); 
%clear all;load de_ct_setup_noiseless.mat; %noiseless setup
%clear all;load de_ct_low1024_setup.mat;
%regbeta = [2^(-5) 2^(-5)]; reg_sel = 0;%old regularization
regbeta = [2^(-5) 2^(-5)]; reg_sel = 1;%modified regularization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conventional decomposition %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[shat_decomp cost_decomp] = de_ftab_s_iter(ftab.fit,fhat.raw,'data',ymi,'niter',10);  %by (W)LS
shat_decomp = reshape(shat_decomp,[n.b n.a 2]);
%Post smoothing
f.fbp_kernel = [1 2 1]';
f.fbp_kernel = f.fbp_kernel / sum(f.fbp_kernel);
shat_decomp(:,:,1) = conv2(shat_decomp(:,:,1), f.fbp_kernel, 'same'); 
shat_decomp(:,:,2) = conv2(shat_decomp(:,:,2), f.fbp_kernel, 'same');
shat_decomp = reshape(shat_decomp, [], 2); %to initialize PWLS and PL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two statistical restoration %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[shat_pwls cost_eval_pwls nrms_pwls] = de_ftab_s_pwls(ftab.fit,fhat.raw,ymi,strue,'init',shat_decomp,...
                                                      'niter',300,'beta',regbeta,'regsel',reg_sel); 
[shat_pl cost_eval_pl nrms_pl] = de_ftab_s_pl(ftab.fit,ftab.xray,ftab.mac.mac,fhat.raw,ymi,strue,'init',shat_decomp,...
                                              'niter',300,'beta',regbeta,'curvtype','pc','regsel',reg_sel);

shat_decomp = reshape(shat_decomp,[n.b n.a 2]);
shat_pwls = reshape(shat_pwls,[n.b n.a 2]); 
shat_pl = reshape(shat_pl,[n.b n.a 2]); 

%%%%%%%%%%%%%%%%%%
% Post-smoothing %
%%%%%%%%%%%%%%%%%%
post_kernel = [0.8 2 0.8]';
post_kernel = post_kernel / sum(post_kernel);
if reg_sel == 1 %modified penalty with post-smoothing
    shat_pwls(:,:,1) = conv2(shat_pwls(:,:,1), post_kernel, 'same');
    shat_pwls(:,:,2) = conv2(shat_pwls(:,:,2), post_kernel, 'same');
    shat_pl(:,:,1) = conv2(shat_pl(:,:,1), post_kernel, 'same');
    shat_pl(:,:,2) = conv2(shat_pl(:,:,2), post_kernel, 'same');
end

%Cost functions    ?? a mismatch in cost function ??
figure;
subplot(1,2,1);plot(cost_eval_pwls(2:end)-max(cost_eval_pwls),'r-o','linewidth',1);grid on;
title('Cost function of PWLS');xlabel('Iteration #');ylabel('\Phi(x^n)-\Phi(x^0)');
subplot(1,2,2);plot(cost_eval_pl(2:end)-max(cost_eval_pl),'b-o','linewidth',1);grid on;
title('Cost function of PL');xlabel('Iteration #');ylabel('\Psi(x^n)-\Psi(x^0)');

%NRMS errors for each iteration
% ?? strange behavior for modified regularizer (ask Fessler)
nrms_pwls_soft = nrms_pwls(:,1);nrms_pwls_bone = nrms_pwls(:,2);
nrms_pl_soft = nrms_pl(:,1);nrms_pl_bone = nrms_pl(:,2);
leng_tmp = size(nrms_pwls,1);
figure;plot((1:leng_tmp),nrms_pwls_soft(1:leng_tmp),'b-.','linewidth',2);
hold on;plot((1:leng_tmp),nrms_pl_soft(1:leng_tmp),'r-','linewidth',2);
legend('PWLS','PL');title('Restored sinogram (Soft tissue)');
xlabel('Number of iterations');ylabel('NRMS error');
figure;plot((1:leng_tmp),nrms_pwls_bone(1:leng_tmp),'b-.','linewidth',2);
hold on;plot((1:leng_tmp),nrms_pl_bone(1:leng_tmp),'r-','linewidth',2);
legend('PWLS','PL');title('Restored sinogram (Bone)');
xlabel('Number of iterations');ylabel('NRMS error');
clear nrms_pwls_soft nrms_pwls_bone nrms_pl_soft nrms_pl_bone leng_tmp

%Restored sinogram plotting
c.s1 = [0 30];c.s2 = [0 20];
figure;im(241,strue(:,:,1),'strue1',c.s1),cbar(c.s1),ylabel('Soft Tissue');
im(245,strue(:,:,2),'strue2',c.s2),cbar(c.s2),ylabel('Bone');
im(242,shat_decomp(:,:,1),'shat1 by decomp',c.s1),cbar(c.s1)
im(246,shat_decomp(:,:,2),'shat2 by decomp',c.s2),cbar(c.s2)
im(243,shat_pwls(:,:,1),'shat1 by pwls',c.s1),cbar(c.s1)
im(247,shat_pwls(:,:,2),'shat2 by pwls',c.s2),cbar(c.s2)
im(244,shat_pl(:,:,1),'shat1 by pl',c.s1),cbar(c.s1)
im(248,shat_pl(:,:,2),'shat2 by pl',c.s2),cbar(c.s2)

%Nrms Error for Sinogram
printm('Sinogram restoration by decomp %0.3g %0.3g', ...
    Nrms(shat_decomp(:,:,1), strue(:,:,1)), Nrms(shat_decomp(:,:,2), strue(:,:,2)))
printm('Sinogram restoration by pwls %0.3g %0.3g', ...
    Nrms(shat_pwls(:,:,1), strue(:,:,1)), Nrms(shat_pwls(:,:,2), strue(:,:,2)))
printm('Sinogram restoration by pl %0.3g %0.3g', ...
    Nrms(shat_pl(:,:,1), strue(:,:,1)), Nrms(shat_pl(:,:,2), strue(:,:,2)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FBP reconstruction from pre-processed data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = fbp2(sg, ig); %for setup
%%%%%%%%%
% DE-CT %
%%%%%%%%%
%The conventional decomposition
xfbp_decomp(:,:,1) = fbp2(shat_decomp(:,:,1), tmp);
xfbp_decomp(:,:,2) = fbp2(shat_decomp(:,:,2), tmp);
xfbp_decomp = max(xfbp_decomp, 0); %negative out
%PWLS
xfbp_pwls(:,:,1) = fbp2(shat_pwls(:,:,1), tmp);
xfbp_pwls(:,:,2) = fbp2(shat_pwls(:,:,2), tmp);
xfbp_pwls = max(xfbp_pwls, 0); 
%PL
xfbp_pl(:,:,1) = fbp2(shat_pl(:,:,1), tmp);
xfbp_pl(:,:,2) = fbp2(shat_pl(:,:,2), tmp);
xfbp_pl = max(xfbp_pl, 0); 

%Reconed image plotting  
figure;im(241,xtrue(:,:,1),'Truth',c.soft),ylabel('Soft Tissue'),cbar([0 1],'h')
im(242,xfbp_decomp(:,:,1),'xhat1 by decomp',c.soft),cbar([0 1],'h')
xlab = sprintf('NRMS = %0.3g',Nrms(xfbp_decomp(:,:,1), xtrue(:,:,1)));xlabel(xlab);
im(243,xfbp_pwls(:,:,1),'xhat1 by pwls',c.soft),cbar([0 1],'h')
xlab = sprintf('NRMS = %0.3g',Nrms(xfbp_pwls(:,:,1), xtrue(:,:,1)));xlabel(xlab);
im(244,xfbp_pl(:,:,1),'xhat1 by pl',c.soft),cbar([0 1],'h')
xlab = sprintf('NRMS = %0.3g',Nrms(xfbp_pl(:,:,1), xtrue(:,:,1)));xlabel(xlab);
im(245,xtrue(:,:,2),'',c.bone),ylabel('Bone'),cbar([0 2],'h')
im(246,xfbp_decomp(:,:,2),'xhat2 by decomp',c.bone),cbar([0 2],'h')
xlab = sprintf('NRMS = %0.3g',Nrms(xfbp_decomp(:,:,2), xtrue(:,:,2)));xlabel(xlab);
im(247,xfbp_pwls(:,:,2),'xhat2 by pwls',c.bone),cbar([0 2],'h')
xlab = sprintf('NRMS = %0.3g',Nrms(xfbp_pwls(:,:,2), xtrue(:,:,2)));xlabel(xlab);
im(248,xfbp_pl(:,:,2),'xhat2 by pl',c.bone),cbar([0 2],'h')
xlab = sprintf('NRMS = %0.3g',Nrms(xfbp_pl(:,:,2), xtrue(:,:,2)));xlabel(xlab);

%Nrms Error for Recon Image
printm('FBP by decomp %0.3g %0.3g', Nrms(xfbp_decomp(:,:,1), xtrue(:,:,1)), Nrms(xfbp_decomp(:,:,2), xtrue(:,:,2)))
printm('FBP by pwls %0.3g %0.3g', Nrms(xfbp_pwls(:,:,1), xtrue(:,:,1)), Nrms(xfbp_pwls(:,:,2), xtrue(:,:,2)))
printm('FBP by pl %0.3g %0.3g', Nrms(xfbp_pl(:,:,1), xtrue(:,:,1)), Nrms(xfbp_pl(:,:,2), xtrue(:,:,2)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAC Computation in image domain %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mac_soft511 = xray_read_atten('soft',511);
mac_bone511 = xray_read_atten('bone',511);
mux_decomp_pet = mac_soft511*xfbp_decomp(:,:,1) + mac_bone511*xfbp_decomp(:,:,2);
mux_pwls_pet = mac_soft511*xfbp_pwls(:,:,1) + mac_bone511*xfbp_pwls(:,:,2);
mux_pl_pet = mac_soft511*xfbp_pl(:,:,1) + mac_bone511*xfbp_pl(:,:,2);
muxtrue_pet = mac_soft511*xtrue(:,:,1) + mac_bone511*xtrue(:,:,2);

im(221,muxtrue_pet,'True LAC(PET)');colorbar horiz;
im(222,mux_decomp_pet,'LAC by de-decomp');colorbar horiz;
xlab = sprintf('NRMS = %0.3g',Nrms(mux_decomp_pet,muxtrue_pet));xlabel(xlab);
im(223,mux_pwls_pet,'LAC by de-pwls');colorbar horiz;
xlab = sprintf('NRMS = %0.3g',Nrms(mux_pwls_pet,muxtrue_pet));xlabel(xlab);
im(224,mux_pl_pet,'LAC by de-pl');colorbar horiz;
xlab = sprintf('NRMS = %0.3g',Nrms(mux_pl_pet,muxtrue_pet));xlabel(xlab);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attenuation Correction Factor I % Only direct methods w/o bilinear scaling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%acfs
acf_true = mac_soft511*strue(:,:,1) + mac_bone511*strue(:,:,2);
acf_shat_decomp = mac_soft511*shat_decomp(:,:,1) + mac_bone511*shat_decomp(:,:,2);
acf_shat_pwls = mac_soft511*shat_pwls(:,:,1) + mac_bone511*shat_pwls(:,:,2);
acf_shat_pl = mac_soft511*shat_pl(:,:,1) + mac_bone511*shat_pl(:,:,2);

figure;
im(221,acf_true,'True ACF');colorbar;
im(222,acf_shat_decomp,'ACF by de-decomp');colorbar;
xlab = sprintf('NRMS = %0.3g',Nrms(exp(acf_shat_decomp),exp(acf_true)));xlabel(xlab);
im(223,acf_shat_pwls,'ACF by de-pwls');colorbar;
xlab = sprintf('NRMS = %0.3g',Nrms(exp(acf_shat_pwls),exp(acf_true)));xlabel(xlab);
im(224,acf_shat_pl,'ACF by de-pl');colorbar;
xlab = sprintf('NRMS = %0.3g',Nrms(exp(acf_shat_pl),exp(acf_true)));xlabel(xlab);

%Nrms Error for ACF
printm('NRMS of ACF by decomp : %0.3g', Nrms(exp(acf_shat_decomp),exp(acf_true)))
printm('NRMS of ACF by pwls : %0.3g', Nrms(exp(acf_shat_pwls),exp(acf_true)))
printm('NRMS of ACF by pl : %0.3g', Nrms(exp(acf_shat_pl),exp(acf_true)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combining with PET model %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load xtrue_pet.mat; %xtrue meas for pet
xtrue_pet = xtrue(:,:,1) + 0.1*xtrue(:,:,2);
ker = [1/2 1 1/2;1 1 1;1/2 1 1/2];
%ker = ones(5,5);
ker = ker / sum(ker(:));
xtrue_pet = conv2(xtrue_pet,ker,'same');
temp = reshape(xtrue_pet, n.x*n.y, 1);
strue_pet = reshape(G * temp(mask,:), [n.b n.a]);
%DE-CT
lam_fbp_decomp = max(fbp2(exp(acf_shat_decomp-acf_true).*strue_pet,tmp),0);
lam_fbp_pwls = max(fbp2(exp(acf_shat_pwls-acf_true).*strue_pet,tmp),0);
lam_fbp_pl = max(fbp2(exp(acf_shat_pl - acf_true).*strue_pet,tmp),0);

figure;im(xtrue_pet,'True PET Image',c.soft);cbar([0 1]);
figure;im(lam_fbp_decomp,'PET Image with CTAC by Conv. Decomp',c.soft);cbar([0 1]);;
xlab = sprintf('NRMS = %0.3g',Nrms(lam_fbp_decomp,xtrue_pet));xlabel(xlab);
figure;im(lam_fbp_pwls,'PET Image with CTAC by PWLS',c.soft);cbar([0 1]);
xlab = sprintf('NRMS = %0.3g',Nrms(lam_fbp_pwls,xtrue_pet));xlabel(xlab);
figure;im(lam_fbp_pl,'PET Image with CTAC by PL',c.soft);cbar([0 1]);
xlab = sprintf('NRMS = %0.3g',Nrms(lam_fbp_pl,xtrue_pet));xlabel(xlab);

save de_sino_results_old.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%
%horizontal profile %
%%%%%%%%%%%%%%%%%%%%%
if isinf(f.scale)
    hori_index = 52;%52
    figure;subplot(1,2,1);imagesc(xtrue(:,:,1)'+xtrue(:,:,2)');colormap gray;colorbar;hold on;
    plot(hori_index*ones(size(xtrue,1)),'r--');hold off;title('True density');
    xtrue_profile = xtrue(:,hori_index,1)+xtrue(:,hori_index,2);
    xfbp_profile_decomp = xfbp_decomp(:,hori_index,1)+xfbp_decomp(:,hori_index,2);
    xfbp_profile_pwls = xfbp_pwls(:,hori_index,1)+xfbp_pwls(:,hori_index,2);
    xfbp_profile_pl = xfbp_pl(:,hori_index,1)+xfbp_pl(:,hori_index,2);
    subplot(122);plot(xtrue_profile,'b','linewidth',2);hold on;
    plot(xfbp_profile_decomp,'c-.','linewidth',2);hold on;
    plot(xfbp_profile_pwls,'y--','linewidth',2);hold on;
    plot(xfbp_profile_pl,'r','linewidth',2);hold off
    legend('True','Decomp-FBP','PWLS-FBP','PL-FBP','location','northwest');title('Noiseless Horizontal profile');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bilinear Transformation % SE-CT only
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pixel-wise transform / CT eff-energy 57 kVp
%From Burger's paper
rhocb = 1.919;
mucb = xray_read_atten('bone',se_ftab.xray.eff)*rhocb; %LAC of bone at CT energy
muw = xray_read_atten('water',se_ftab.xray.eff); %LAC of water at CT energy
muw_pet = xray_read_atten('water',511);
maccb_pet = xray_read_atten('bone',511);
mac_soft = xray_read_atten('soft',se_ftab.xray.eff);

%Conversion to CT number [HU]
se_hu(:,:,1) = (mac_soft*se_xfbp_decomp - muw)/muw * 1000;
se_hu(:,:,2) = (mac_soft*se_xfbp_pwls - muw)/muw *1000;
se_hu(:,:,3) = (mac_soft*se_xfbp_pl - muw)/muw *1000;
HUcb = (mucb - muw)/muw * 1000;

%Bilinear scaling and forward projection
for pp=1:3
    temp = se_hu(:,:,pp);
    se_mux_petn = (1 + temp/1000)*muw_pet;    
    se_mux_petp = (1 + (rhocb *maccb_pet/muw_pet - 1)*temp/HUcb)*muw_pet;
    ind = find(temp > 0);    
    se_mux_petn(ind) = se_mux_petp(ind);
    se_mux_pet(:,:,pp) = se_mux_petn; %mux at pet energy
    se_fb(:,:,pp) = reshape(G*col(se_mux_petn),[n.b n.a]); 
    %Think why bilinear scaling is needed??
    %Forward projection may not be necessary??
    %Compare with PET reconstructed Image !!
    %You are HERE !!
end

