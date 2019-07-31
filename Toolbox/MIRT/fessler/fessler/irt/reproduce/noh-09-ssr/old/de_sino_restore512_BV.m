%This is the code for DE-CT sinogram restoration
%Modified from spie02.m

%Iodine Contrast
%Bright Streak in DECT (w/ modified reg)
%Algorithm Check

clear all;
im off;
tic
%%%%%%%%%%
%%%%%%%%%%
% PART I % setup for data generation
%%%%%%%%%%
%%%%%%%%%%

%Geometry selection
geo_sel = 0; %0|1 (para|fan-beam)
%Iodine selection
Iodine_sel = 1

load xtrue512.mat; %512 by 512
[n.x n.y] = size(xtrue(:,:,1));
%mask
mask = sum(xtrue, 3) > 0;
mask = imdilate(mask, strel('disk', 5));
mask = logical(mask);

%%%%%%%%%%%%%%%%
%System matrix % for forward projection
%%%%%%%%%%%%%%%%
if geo_sel == 1 %fan-beam geometry
    f.dx = 0.1;%in cm
    n.bf = 888;n.af = 984;ds = 0.1;
    n.bp = 1024;n.ap = 800;dr = ds/2;
    %Image geometry
    ig = image_geom('nx', n.x, 'ny', n.y, 'dx', f.dx, 'mask', mask);
    %Fan-beam geometry
    sgf = sino_geom('fan', 'nb', n.bf, 'na', n.af, 'ds', ds, 'dsd', 949,...
                    'dod', 408, 'offset_s', 1.25); 
    %Parallel-beam geometry
    sgp = sino_geom('par', 'nb', n.bp, 'na', n.ap, 'dr', dr, ...
                    'offset_r', 0, 'orbit', 180); 
    %G = Gtomo2_strip(sg, ig); %small CT
    G = Gtomo2_dscmex(sgf, ig); %full CT (Fan-beam)
    gi = reshape(sum(G'), [n.bf n.af]); %simple projection

elseif geo_sel == 0 %parallel-beam geometry
    f.dx = 0.1; %in cm
    n.bf = 1024;n.af = 800;ds = 0.1;
    n.bp = n.bf;n.ap = n.af;dr = ds/2;
    %Image geometry
    ig = image_geom('nx', n.x, 'ny', n.y, 'dx', f.dx, 'mask', mask);
    %Parallel-beam geometry
    sgp = sino_geom('par', 'nb', n.bp, 'na', n.ap, 'dr', dr, ...
                    'offset_r', 0, 'orbit', 180); 
    G = Gtomo2_dscmex(sgp, ig); %full CT (Para-beam)
    gi = reshape(sum(G'), [n.bf n.af]); %simple projection
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f table and true images %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Iodine_sel == 0
    ftab = de_ftab_build({'ps1'}); %DE-Table
elseif Iodine_sel == 1
    ftab = de_ftab_build({'ps1t'}); %TE-Table
    %rhoI = 4.930; %Pure Iodine density (from Hubbell&Seltzer)
    rhoI = 0.09; %0.05 diluted Iodine density
    ell1 = [-10 0 0.3 0.5 0 rhoI];
    ell2 = [0 0 0.5 0.3 0 rhoI];
    xtrue(:,:,3) = ellipse_im(ig,ell1) + ellipse_im(ig,ell2);
end %iodine with soft tissue

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ideal and Noisy Measurement % y_{mi}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Noiseless sinogram and f
%
xtrueL = size(xtrue,3);
tmp = reshape(xtrue, n.x*n.y, xtrueL);
tic; %noiseless forward projection
strue = reshape(G * tmp(mask,:), [n.bf n.af xtrueL]);
printm('forward projection time %0.3g', toc)
ftrue = de_ftab_fm(strue, ftab.mac.mac, ftab.xray.Ide); %noiseless f
%ftrue = ftab.fm_fun(ftab, {strue(:,:,1),strue(:,:,2)}); %Fessler's
%ftrue = ftab.fm_fun(ftab, {strue(:,:,1),strue(:,:,2),strue(:,:,3)});

%Noisy measurement
ybi = zeros(n.bf, n.af, 2); %bar{y}_{mi}
ybi(:,:,1) = ftab.xray.I(1) * exp(-ftrue(:,:,1)); %true mean
ybi(:,:,2) = ftab.xray.I(2) * exp(-ftrue(:,:,2));

%Dose amount
%f.scale_high = inf; f.scale_low = inf;%for noiseless profile
f.scale_high = 5*10^5 / ybi(1,1,2); %high dose
f.scale_low = 2*10^5 / ybi(1,1,2); %low dose

%%%%%%%%%%%
%%%%%%%%%%%
% PART II % Multiple random measurements generation (Poissone distributed)
%%%%%%%%%%%
%%%%%%%%%%% 
leng_iter = 45; %# of iteration for multiple realization of measurements
down = 4;n.bp = n.bp/down;n.ap = n.ap/down; %downsampling
 
for ij=1:leng_iter %for bias & variance estimation
    ij

    if isinf(f.scale_high) || isinf(f.scale_low)
        ymi = ybi; se_ymi = ybi;	% noiseless
    else                                                      %%%%%%%%%%%%%%%%%%%%%%%
        ymi = poisson(f.scale_low*ybi, 0) / f.scale_low;      %corrupted by poisson %
        se_ymi = poisson(f.scale_high*ybi, 0) / f.scale_high; %%%%%%%%%%%%%%%%%%%%%%%
    end %measurements look like sinogram                 

    %Borrow neighbor(s) for any log(0) values
    off = 1; %DECT
    while any(ymi(:) == 0)
        ii = find(ymi(:) == 0);
        printm('fixing %d zeros in ymi with %d', length(ii), off)
        ymi(ii) = max(ymi(ii+off), ymi(ii-off));
        off = off + 1;
    end, clear off ii
    off = 1; %SECT
    while any(se_ymi(:) == 0)
        ii = find(se_ymi(:) == 0);
        printm('fixing %d zeros in se_ymi with %d', length(ii), off)
        se_ymi(ii) = max(se_ymi(ii+off), se_ymi(ii-off));
        off = off + 1;
    end, clear off ii

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Geometry Tx & downsampling %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if geo_sel == 1 %fan-beam
        ymipd(:,:,1) = downsample2(rebin_fan2par(ymi(:,:,1), sgf, sgp),down); %fan -> par
        ymipd(:,:,2) = downsample2(rebin_fan2par(ymi(:,:,2), sgf, sgp),down); %[n.bp n.ap]
        se_ymipd(:,:,1) = downsample2(rebin_fan2par(se_ymi(:,:,1), sgf, sgp),down); 
        se_ymipd(:,:,2) = downsample2(rebin_fan2par(se_ymi(:,:,2), sgf, sgp),down);
        for pp=1:xtrueL
            strue_pd(:,:,pp) = downsample2(rebin_fan2par(strue(:,:,pp), sgf, sgp),down);
        end
    elseif geo_sel == 0 %parallel-beam
        ymipd(:,:,1) = downsample2(ymi(:,:,1),down);  %downsampling
        ymipd(:,:,2) = downsample2(ymi(:,:,2),down); %[n.bp n.ap]
        se_ymipd(:,:,1) = downsample2(se_ymi(:,:,1),down); 
        se_ymipd(:,:,2) = downsample2(se_ymi(:,:,2),down); 
        for pp=1:xtrueL
            strue_pd(:,:,pp) = downsample2(strue(:,:,pp),down);
        end
    end


    %%%%%%%%%%%%
    % PART III %
    %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data Analysis from here % using down-sampled sinogram mesurements !!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fitting Table Generation for DE-CT Restoration %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Iodine_sel == 1, clear ftab;ftab = de_ftab_build({'ps1'});  end   %DE_Table recall
    LL = size(ftab.mac.mac,2); %# of decomposed materials

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % regularization design % reg_sel 0|1 (old|modified)
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %regbeta = 2^(-5)*ones(1,LL); reg_sel = 0 
    regbeta = 2^(-8)*ones(1,LL); reg_sel = 1
    %# of iterations
    if reg_sel == 0,    de_niter = 100; end
    if reg_sel == 1,    de_niter = 100; end
    %Kernel for smoothing
    f.fbp_kernel = [1 1 1]'; %Presmoothing
    f.fbp_kernel = f.fbp_kernel / sum(f.fbp_kernel);
    %ymipd(:,:,1) = conv2(ymipd(:,:,1), f.fbp_kernel, 'same'); 
    %ymipd(:,:,2) = conv2(ymipd(:,:,2), f.fbp_kernel, 'same');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DECT sinogram restoration %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fhat eq(19) 
    fhat.raw(:,:,1) = -log(ymipd(:,:,1) / ftab.xray.I(1)); 
    fhat.raw(:,:,2) = -log(ymipd(:,:,2) / ftab.xray.I(2));
    fhat.raw(isinf(fhat.raw)) = 0;

    %Conventional decomposition
    [shat_decomp cost_eval_decomp] = de_ftab_s_iter(ftab.fit,fhat.raw,'niter',10,'data',ymipd);%,'init',reshape(strue_pd,[],2));  %by (W)LS
    shat_decomp = reshape(shat_decomp,[n.bp n.ap LL]);
    %Post smoothing
    for ii=1:LL
        shat_decomp(:,:,ii) = conv2(shat_decomp(:,:,ii), f.fbp_kernel, 'same');
    end
    shat_decomp = reshape(shat_decomp, [], LL); %to initialize PWLS and PL

    %Two statistical restoration
    if Iodine_sel == 1, strue_pwls = [];strue_pl = [];  end
    if Iodine_sel == 0, strue_pwls = strue_pd;strue_pl = strue_pd;  end
    [shat_pwls cost_pwls1 cost_pwls2 nrms_pwls] = de_ftab_s_pwls(ftab.fit,fhat.raw,ymipd,'strue',strue_pwls,'init',shat_decomp,...
                                                                 'niter',de_niter,'beta',regbeta,'regsel',reg_sel); 
    [shat_pl cost_pl1 cost_pl2 nrms_pl] = de_ftab_s_pl(ftab.fit,ftab.xray,ftab.mac.mac,fhat.raw,ymipd,'strue',strue_pl,'init',shat_decomp,...
                                                       'niter',de_niter,'beta',regbeta,'curvtype','pc','regsel',reg_sel);

    shat_decomp = reshape(shat_decomp,[n.bp n.ap LL]);
    shat_pwls = reshape(shat_pwls,[n.bp n.ap LL]); 
    shat_pl = reshape(shat_pl,[n.bp n.ap LL]); 
    cost_eval_pwls = cost_pwls1+cost_pwls2;
    cost_eval_pl = cost_pl1+cost_pl2;

    %Post-smoothing (for DECT PWLS and PL w/ modified penalty)
    post_kernel = [0.5 2 0.5]';
    post_kernel = post_kernel / sum(post_kernel);
    if reg_sel == 1 
        for ii=1:LL
            shat_pwls(:,:,ii) = conv2(shat_pwls(:,:,ii), post_kernel, 'same');
            shat_pl(:,:,ii) = conv2(shat_pl(:,:,ii), post_kernel, 'same');
        end
    end

    %Individual terms of PL cost functions
    cost_pl1 = cost_pl1 - cost_pl1(1);
    cost_pl2 = cost_pl2 - cost_pl2(1);
    %subplot(121);plot(cost_pl1,'r-o');grid on;
    %subplot(122);plot(cost_pl2,'b-o');grid on;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SECT sinogram restoration % discard the 1st measurements / assuming soft tissue
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    se_ftab = de_ftab_build({'ps1s'}); %fitting for se-ct
    se_strue_pd = sum(strue_pd,3);
    se_fhat.raw(:,:,1) = -log(se_ymipd(:,:,1) / ftab.xray.I(1)); 
    se_fhat.raw(:,:,2) = -log(se_ymipd(:,:,2) / ftab.xray.I(2));
    se_fhat.raw(isinf(se_fhat.raw)) = 0;

    disp ('Single-kVp CT')
    [se_shat_decomp se_cost_eval_decomp] = de_ftab_s_iter(se_ftab.fit,se_fhat.raw(:,:,2),...
                                                          'niter',10,'data',se_ymipd(:,:,2)); %,'init',reshape(se_strue_pd,[],1)); %by (W)LS
    se_shat_decomp = reshape(se_shat_decomp,[n.bp n.ap]);
    se_shat_decomp = conv2(se_shat_decomp, f.fbp_kernel, 'same'); %post smoothing
    se_shat_decomp = col(se_shat_decomp);

    %Two statistical restorations
    disp ('Single-kVp CT')
    [se_shat_pwls se_cost_pwls1 se_cost_pwls2 se_nrms_pwls] = de_ftab_s_pwls(se_ftab.fit,se_fhat.raw(:,:,2),se_ymipd(:,:,2),'strue',se_strue_pd,'init',se_shat_decomp,...
                                                                             'niter',de_niter,'beta',regbeta(2),'regsel',reg_sel);
    disp ('Single-kVp CT')
    [se_shat_pl se_cost_pl1 se_cost_pl2 se_nrms_pl] = de_ftab_s_pl(se_ftab.fit,se_ftab.xray,se_ftab.mac.mac,se_fhat.raw(:,:,2),se_ymipd(:,:,2),...
                                                                   'strue',se_strue_pd,'init',se_shat_decomp,'niter',de_niter,'beta',regbeta(2),...
                                                                   'curvtype','pc','regsel',reg_sel); %by PL                                                                                                                                

    se_shat_decomp = reshape(se_shat_decomp,[n.bp n.ap]);
    se_shat_pwls = reshape(se_shat_pwls,[n.bp n.ap]);
    se_shat_pl = reshape(se_shat_pl,[n.bp n.ap]);
    se_cost_eval_pwls = se_cost_pwls1+se_cost_pwls2;
    se_cost_eval_pl = se_cost_pl1+se_cost_pl2;

    %Post-smoothing (for SECT)
    if reg_sel == 1 %modified penalty with post-smoothing
        se_shat_pwls = conv2(se_shat_pwls, post_kernel, 'same');
        se_shat_pl = conv2(se_shat_pl, post_kernel, 'same');
    end

    %Individual terms of PL cost functions
    se_cost_pl1 = se_cost_pl1 - se_cost_pl1(1);
    se_cost_pl2 = se_cost_pl2 - se_cost_pl2(1);
    %subplot(121);plot(se_cost_pl1,'r-o');grid on;
    %subplot(122);plot(se_cost_pl2,'b-o');grid on;

    %%%%%%%%%%%%%%%%%%%%%%
    % FBP reconstruction % to full size CT
    %%%%%%%%%%%%%%%%%%%%%%
    mask_d = logical(downsample2(mask,down));
    ig_d = image_geom('nx', n.x/down, 'ny', n.y/down, 'dx', f.dx*down,...
                      'mask', mask_d);
    sgp_d = sino_geom('par', 'nb', n.bp, 'na', n.ap,'dr', dr*down, ...
                      'offset_r', 0, 'orbit', 180); 
    fbp2tmp = fbp2(sgp_d, ig_d); %backward
    G_d = Gtomo2_dscmex(sgp_d, ig_d); %forward

    %DECT
    for ii=1:LL
        xfbp_decomp(:,:,ii) = fbp2(shat_decomp(:,:,ii), fbp2tmp);
        xfbp_pwls(:,:,ii) = fbp2(shat_pwls(:,:,ii), fbp2tmp);
        xfbp_pl(:,:,ii) = fbp2(shat_pl(:,:,ii), fbp2tmp);
    end
    xfbp_decomp = max(xfbp_decomp, 0); 
    xfbp_pwls = max(xfbp_pwls, 0); 
    xfbp_pl = max(xfbp_pl, 0); 

    %SECT
    se_xfbp_decomp = fbp2(se_shat_decomp, fbp2tmp);
    se_xfbp_decomp = max(se_xfbp_decomp, 0);
    se_xfbp_pwls = fbp2(se_shat_pwls, fbp2tmp);
    se_xfbp_pwls = max(se_xfbp_pwls, 0);
    se_xfbp_pl = fbp2(se_shat_pl, fbp2tmp);
    se_xfbp_pl = max(se_xfbp_pl,0);

    %%%%%%%%%%%%%%%%%%
    % ACF estimation % for DE-CT
    %%%%%%%%%%%%%%%%%%
    mac_511(1) = xray_read_atten('soft',511); %soft tissue
    mac_511(2) = xray_read_atten('bone',511); %bone
    mac_511(3) = xray_read_atten('iodine',511); %iodine
    acftrue = mac_511(1)*strue_pd(:,:,1) + mac_511(2)*strue_pd(:,:,2);
    if Iodine_sel == 1, acftrue = acftrue + mac_511(3)*strue_pd(:,:,3); end
    acf_decomp = mac_511(1)*shat_decomp(:,:,1) + mac_511(2)*shat_decomp(:,:,2);
    acf_pwls = mac_511(1)*shat_pwls(:,:,1) + mac_511(2)*shat_pwls(:,:,2);
    acf_pl = mac_511(1)*shat_pl(:,:,1) + mac_511(2)*shat_pl(:,:,2);

    %%%%%%%%%%%%%%%%%%%%
    % Bilinear Scaling % for SECT
    %%%%%%%%%%%%%%%%%%%%
    rhocb = 1.919;
    mucb = xray_read_atten('bone',se_ftab.xray.eff)*rhocb; %LAC of c-bone
    muw = xray_read_atten('water',se_ftab.xray.eff); %LAC of water
    muw_pet = xray_read_atten('water',511);
    maccb_pet = xray_read_atten('bone',511);
    mac_soft = xray_read_atten('soft',se_ftab.xray.eff);

    %Conversion to CT number [HU] (?? multiplied by soft mac ?? Check)
    HU(:,:,1) = (mac_soft*se_xfbp_decomp - muw)/muw * 1000;
    HU(:,:,2) = (mac_soft*se_xfbp_pwls - muw)/muw *1000;
    HU(:,:,3) = (mac_soft*se_xfbp_pl - muw)/muw *1000;
    HUcb = (mucb - muw)/muw * 1000;

    %Bilinear scaling and Forward projection
    for pp=1:3
        mux_petn = (1 + HU(:,:,pp)/1000)*muw_pet;    
        mux_petp = (1 + (rhocb *maccb_pet/muw_pet - 1)*HU(:,:,pp)/HUcb )*muw_pet; 
        ind = find(HU(:,:,pp) > 0);    
        mux_petn(ind) = mux_petp(ind);
        mux_pet(:,:,pp) = mux_petn; 
        tmp = reshape(mux_pet(:,:,pp), [], 1);
        se_acf(:,:,pp) = reshape(G_d*tmp(mask_d), [n.bp n.ap]); 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PET/CT image reconstruction %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xtrue_tmp = xtrue(:,:,1);
    xtrue_pet = xtrue_tmp + 0.1*xtrue(:,:,2); %original pet lambda
    if Iodine_sel == 1
        xtrue_pet = xtrue_pet + xtrue(:,:,3)/max(col(xtrue(:,:,3)));%.*mean(col(xtrue_tmp(xtrue_tmp~=0)));
    end
    xtrue_petd = downsample2(xtrue_pet,down); %downsampled pet lambda
    %PET sinogram
    temp = reshape(xtrue_pet, [], 1);
    strue_pet = reshape(G * temp(mask), [down*n.bp down*n.ap]);
    strue_pet = downsample2(strue_pet,down);
    %PET measurements (noisyless)
    ymi_pet = strue_pet .*exp(-acftrue); 

    %Attenuation correction (DECT)
    shat_pet_decomp = ymi_pet.*exp(acf_decomp);
    shat_pet_pwls = ymi_pet.*exp(acf_pwls);
    shat_pet_pl = ymi_pet.*exp(acf_pl);

    %Attenuation correction (SECT by bilinear scaling)
    se_shat_pet_decomp = ymi_pet.*exp(se_acf(:,:,1));
    se_shat_pet_pwls = ymi_pet.*exp(se_acf(:,:,2));
    se_shat_pet_pl = ymi_pet.*exp(se_acf(:,:,3));

    %FBP reconstruction for PET/CT
    fbp2tmp = fbp2(sgp_d, ig_d);
    xfbp_pet_decomp(:,:,ij) = max(fbp2(shat_pet_decomp,fbp2tmp),0);
    xfbp_pet_pwls(:,:,ij) = max(fbp2(shat_pet_pwls,fbp2tmp),0);
    xfbp_pet_pl(:,:,ij) = max(fbp2(shat_pet_pl,fbp2tmp),0);
    se_xfbp_pet_decomp(:,:,ij) = max(fbp2(se_shat_pet_decomp,fbp2tmp),0);
    se_xfbp_pet_pwls(:,:,ij) = max(fbp2(se_shat_pet_pwls,fbp2tmp),0);
    se_xfbp_pet_pl(:,:,ij) = max(fbp2(se_shat_pet_pl,fbp2tmp),0);
    
end %end of the loop for multiple measurement realizations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bias and Variance Estimation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bias
bias_pet_decomp = mean(xfbp_pet_decomp,3)-xtrue_petd;
bias_pet_pwls = mean(xfbp_pet_pwls,3)-xtrue_petd;
bias_pet_pl = mean(xfbp_pet_pl,3)-xtrue_petd;
bias_se_pet_decomp = mean(se_xfbp_pet_decomp,3)-xtrue_petd;
bias_se_pet_pwls = mean(se_xfbp_pet_pwls,3)-xtrue_petd;
bias_se_pet_pl = mean(se_xfbp_pet_pl,3)-xtrue_petd;

%variance
L = size(xfbp_pet_decomp,3);
var_pet_decomp = mean((xfbp_pet_decomp-repmat(mean(xfbp_pet_decomp,3),[1 1 L])).^2,3);
var_pet_pwls = mean((xfbp_pet_pwls-repmat(mean(xfbp_pet_pwls,3),[1 1 L])).^2,3);
var_pet_pl = mean((xfbp_pet_pl-repmat(mean(xfbp_pet_pl,3),[1 1 L])).^2,3);
var_se_pet_decomp = mean((se_xfbp_pet_decomp-repmat(mean(se_xfbp_pet_decomp,3),[1 1 L])).^2,3);
var_se_pet_pwls = mean((se_xfbp_pet_pwls-repmat(mean(se_xfbp_pet_pwls,3),[1 1 L])).^2,3);
var_se_pet_pl = mean((se_xfbp_pet_pl-repmat(mean(se_xfbp_pet_pl,3),[1 1 L])).^2,3);

%mean
mean_pet_decomp = mean(xfbp_pet_decomp,3);
mean_pet_pwls = mean(xfbp_pet_pwls,3);
mean_pet_pl = mean(xfbp_pet_pl,3);
mean_se_pet_decomp = mean(se_xfbp_pet_decomp,3);
mean_se_pet_pwls = mean(se_xfbp_pet_pwls,3);
mean_se_pet_pl = mean(se_xfbp_pet_pl,3);

%standard deviation
std_pet_decomp = sqrt(var_pet_decomp);
std_pet_pwls = sqrt(var_pet_pwls);
std_pet_pl = sqrt(var_pet_pl);
std_se_pet_decomp = sqrt(var_se_pet_decomp);
std_se_pet_pwls = sqrt(var_se_pet_pwls);
std_se_pet_pl = sqrt(var_se_pet_pl);

time = toc
% Data save
filename = sprintf('sino_restore_results512_BV_R%d_I%d.mat',reg_sel,Iodine_sel);
save(filename);

%%%%%%%%%%%%%%%%%%%%%%%
%Plotting and Display %
%%%%%%%%%%%%%%%%%%%%%%%
%Bias
im(231,bias_pet_decomp);cbar;
im(232,bias_pet_pwls);cbar;
im(233,bias_pet_pl);cbar;
im(234,bias_se_pet_decomp);cbar;
im(235,bias_se_pet_pwls);cbar;
im(236,bias_se_pet_pl);cbar;

% %NRSSB (normalized root sum of squared biases)
% sqrt(sum(col(bias_pet_decomp).^2)) / sqrt(sum(col(xtrue_petd).^2)) * 100
% sqrt(sum(col(bias_pet_pwls).^2)) / sqrt(sum(col(xtrue_petd))) * 100
% sqrt(sum(col(bias_pet_pl).^2)) / sqrt(sum(col(xtrue_petd))) * 100
% sqrt(sum(col(bias_se_pet_decomp).^2)) / sqrt(sum(col(xtrue_petd))) * 100
% sqrt(sum(col(bias_se_pet_pwls).^2)) / sqrt(sum(col(xtrue_petd))) * 100
% sqrt(sum(col(bias_se_pet_pl).^2)) / sqrt(sum(col(xtrue_petd))) * 100

%NRMSB (normalized root mean squared biase)
sqrt(mean(col(bias_pet_decomp).^2)) / mean(col(xtrue_petd)) * 100
sqrt(mean(col(bias_pet_pwls).^2)) / mean(col(xtrue_petd)) * 100
sqrt(mean(col(bias_pet_pl).^2)) / mean(col(xtrue_petd)) * 100
sqrt(mean(col(bias_se_pet_decomp).^2)) / mean(col(xtrue_petd)) * 100
sqrt(mean(col(bias_se_pet_pwls).^2)) / mean(col(xtrue_petd)) * 100
sqrt(mean(col(bias_se_pet_pl).^2)) / mean(col(xtrue_petd)) * 100

%Std
im(231,std_pet_decomp);cbar;
im(232,std_pet_pwls);cbar;
im(233,std_pet_pl);cbar;
im(234,std_se_pet_decomp);cbar;
im(235,std_se_pet_pwls);cbar;
im(236,std_se_pet_pl);cbar;

% %NRSV (normalized root sum of variances)
% sqrt(sum(col(var_pet_decomp))) / sqrt(sum(col(xtrue_petd).^2)) * 100
% sqrt(sum(col(var_pet_pwls))) / sqrt(sum(col(xtrue_petd).^2)) * 100
% sqrt(sum(col(var_pet_pl))) / sqrt(sum(col(xtrue_petd).^2)) * 100
% sqrt(sum(col(var_se_pet_decomp))) / sqrt(sum(col(xtrue_petd).^2)) * 100
% sqrt(sum(col(var_se_pet_pwls))) / sqrt(sum(col(xtrue_petd).^2)) * 100
% sqrt(sum(col(var_se_pet_pl))) / sqrt(sum(col(xtrue_petd).^2)) * 100

%NRMV (normalized root mean variance)
sqrt(mean(col(var_pet_decomp))) / mean(col(xtrue_petd)) * 100
sqrt(mean(col(var_pet_pwls))) / mean(col(xtrue_petd)) * 100
sqrt(mean(col(var_pet_pl))) / mean(col(xtrue_petd)) * 100
sqrt(mean(col(var_se_pet_decomp))) / mean(col(xtrue_petd)) * 100
sqrt(mean(col(var_se_pet_pwls))) / mean(col(xtrue_petd)) * 100
sqrt(mean(col(var_se_pet_pl))) / mean(col(xtrue_petd)) * 100







