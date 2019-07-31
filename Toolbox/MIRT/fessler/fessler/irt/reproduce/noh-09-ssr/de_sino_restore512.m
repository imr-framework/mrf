% de_sino_restore512.m
% This is the code for DE-CT sinogram restoration
% Modified from spie02.m by Joonki Noh

%Iodine Contrast
%Bright Streak in DECT (w/ modified reg)
%Algorithm Check

if ~isvar('ymi')
%	geo_sel = 1; %0|1 (para|fan-beam)
%	Iodine_sel = 0; %Iodine selection
	de_sino_restore512_setup
end

if Iodine_sel == 1
	ftab = de_ftab_build({'ps1'});
end % DE_Table recall
LL = size(ftab.mac.mac,2); % # of decomposed materials

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Analysis from here % using down-sampled sinogram mesurements !!
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% regularization design % reg_sel 0|1 (old|modified)
%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isvar('shat_decomp')
	% regbeta = 2^(-5)*ones(1,LL); reg_sel = 0
	regbeta = 2^(-8)*ones(1,LL); reg_sel = 1;
	% # of iterations
	if reg_sel == 0, de_niter = 100; end
	if reg_sel == 1, de_niter = 100; end
	% Kernel for smoothing
	f.fbp_kernel = [1 1 1]'; % Presmoothing
	f.fbp_kernel = f.fbp_kernel / sum(f.fbp_kernel);
	% ymipd(:,:,1) = conv2(ymipd(:,:,1), f.fbp_kernel, 'same');
	% ymipd(:,:,2) = conv2(ymipd(:,:,2), f.fbp_kernel, 'same');

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DECT sinogram restoration %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% fhat eq(19)
	fhat.raw(:,:,1) = -log(ymipd(:,:,1) / ftab.xrs.I(1));
	fhat.raw(:,:,2) = -log(ymipd(:,:,2) / ftab.xrs.I(2));
	fhat.raw(isinf(fhat.raw)) = 0;

	% Conventional decomposition
	% [shat_decomp cost_eval_decomp] = de_ftab_s_iter(ftab.fit,fhat.raw, 'niter', 10, 'data', ymipd);
	shat_decomp = ftab.inv2.fun(fhat.raw);
end

if 0 % Post smoothing
	for ii=1:LL
		shat_decomp(:,:,ii) = conv2(shat_decomp(:,:,ii), f.fbp_kernel, 'same');
	end
	shat_decomp = reshape(shat_decomp, [], LL); % to initialize PWLS and PL
end

% Two statistical restoration
if 1 % WH not needed
	if Iodine_sel == 1, strue_pwls = []; strue_pl = []; end
	if Iodine_sel == 0, strue_pwls = strue_pd; strue_pl = strue_pd; end
% todo: 2014-03-30 email about error
	[shat_pwls cost_pwls1 cost_pwls2 nrms_pwls] = ...
		de_ftab_s_pwls(ftab.fit, fhat.raw, ymipd, ...
			'strue', strue_pwls,'init', shat_decomp,...
			'niter', de_niter, 'beta', regbeta, 'regsel', reg_sel);
return
	[shat_pl cost_pl1 cost_pl2 nrms_pl] = ...
		de_ftab_s_pl(ftab.fit, ftab.xray, ftab.mac.mac, fhat.raw, ...
			ymipd, 'strue', strue_pl, 'init', shat_decomp, ...
			'niter', de_niter, 'beta', regbeta, ...
			'curvtype', 'pc', 'regsel', reg_sel);
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
end

if 0
%Cost functions (DECT) ?? a mismatch in cost function ??
figure;steptmp = 1;
subplot(211);plot((2:steptmp:de_niter/4),cost_eval_pwls(2:steptmp:end/4)-max(cost_eval_pwls),'r-o','linewidth',2);%grid on;
title('Cost function of PWLS');xlabel('Iteration #');ylabel('\Phi(x^n)-\Phi(x^0)');
subplot(212);plot((2:steptmp:de_niter/4),cost_eval_pl(2:steptmp:end/4)-max(cost_eval_pl),'b-o','linewidth',2);%grid on;
title('Cost function of PL');xlabel('Iteration #');ylabel('\Psi(x^n)-\Psi(x^0)');

%Individual terms of PL cost functions
cost_pl1 = cost_pl1 - cost_pl1(1);
cost_pl2 = cost_pl2 - cost_pl2(1);
%subplot(121);plot(cost_pl1,'r-o');grid on;
%subplot(122);plot(cost_pl2,'b-o');grid on;

%NRMS errors for each iteration (DECT)
% ?? strange behavior for modified regularizer (ask Fessler)
nrms_pwls_soft = nrms_pwls(:,1);nrms_pwls_bone = nrms_pwls(:,2);
nrms_pl_soft = nrms_pl(:,1);nrms_pl_bone = nrms_pl(:,2);
leng_tmp = size(nrms_pwls,1);
figure;subplot(211);plot((1:leng_tmp),nrms_pwls_soft(1:leng_tmp),'b-.','linewidth',2);
hold on;plot((1:leng_tmp),nrms_pl_soft(1:leng_tmp),'r-','linewidth',2);
legend('PWLS','PL');title('Restored sinogram (Soft tissue)');
xlabel('Number of iterations');ylabel('NRMS error');
subplot(212);plot((1:leng_tmp),nrms_pwls_bone(1:leng_tmp),'b-.','linewidth',2);
hold on;plot((1:leng_tmp),nrms_pl_bone(1:leng_tmp),'r-','linewidth',2);
legend('PWLS','PL');title('Restored sinogram (Bone)');
xlabel('Number of iterations');ylabel('NRMS error');
clear nrms_pwls_soft nrms_pwls_bone nrms_pl_soft nrms_pl_bone leng_tmp

%Restored sinogram plotting
figure;
if Iodine_sel == 0
	c.s1 = [0 30];c.s2 = [0 20];
	im(241,strue_pd(:,:,1),'strue1',c.s1),cbar(c.s1,'horiz'),ylabel('Soft Tissue');
	im(245,strue_pd(:,:,2),'strue2',c.s2),cbar(c.s2,'horiz'),ylabel('Bone');
	im(242,shat_decomp(:,:,1),'shat1 by decomp',c.s1),cbar(c.s1,'horiz')
	xlab = sprintf('NRMS = %0.3g',Nrms(shat_decomp(:,:,1), strue_pd(:,:,1)));xlabel(xlab);
	im(246,shat_decomp(:,:,2),'shat2 by decomp',c.s2),cbar(c.s2,'horiz')
	xlab = sprintf('NRMS = %0.3g',Nrms(shat_decomp(:,:,2), strue_pd(:,:,2)));xlabel(xlab);
	im(243,shat_pwls(:,:,1),'shat1 by pwls',c.s1),cbar(c.s1,'horiz')
	xlab = sprintf('NRMS = %0.3g',Nrms(shat_pwls(:,:,1), strue_pd(:,:,1)));xlabel(xlab);
	im(247,shat_pwls(:,:,2),'shat2 by pwls',c.s2),cbar(c.s2,'horiz')
	xlab = sprintf('NRMS = %0.3g',Nrms(shat_pwls(:,:,2), strue_pd(:,:,2)));xlabel(xlab);
	im(244,shat_pl(:,:,1),'shat1 by pl',c.s1),cbar(c.s1,'horiz')
	xlab = sprintf('NRMS = %0.3g',Nrms(shat_pl(:,:,1), strue_pd(:,:,1)));xlabel(xlab);
	im(248,shat_pl(:,:,2),'shat2 by pl',c.s2),cbar(c.s2,'horiz')
	xlab = sprintf('NRMS = %0.3g',Nrms(shat_pl(:,:,2), strue_pd(:,:,2)));xlabel(xlab);
elseif Iodine_sel == 1
	c.s1 = [0 55];
	im(221,sum(strue_pd,3),'strue',c.s1);cbar(c.s1);
	im(222,sum(shat_decomp,3),'shat by de-comp',c.s1);cbar(c.s1);
	xlab = sprintf('NRMS = %0.3g',Nrms(sum(shat_decomp,3), sum(strue_pd,3)));xlabel(xlab);
	im(223,sum(shat_pwls,3),'shat by de-pwls',c.s1);cbar(c.s1);
	xlab = sprintf('NRMS = %0.3g',Nrms(sum(shat_pwls,3), sum(strue_pd,3)));xlabel(xlab);
	im(224,sum(shat_pl,3),'shat by de-pl',c.s1);cbar(c.s1);
	xlab = sprintf('NRMS = %0.3g',Nrms(sum(shat_pl,3), sum(strue_pd,3)));xlabel(xlab);
end

end

if 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECT sinogram restoration % discard the 1st measurements / assuming soft tissue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
se_ftab = de_ftab_build({'ps1s'}); %fitting for se-ct
se_strue_pd = sum(strue_pd,3);
se_fhat.raw(:,:,1) = -log(se_ymipd(:,:,1) / ftab.xrs.I(1));
se_fhat.raw(:,:,2) = -log(se_ymipd(:,:,2) / ftab.xrs.I(2));
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

%Restored sinogram plotting
c.s1 = [0 55];figure;
im(221,se_strue_pd,'strue',c.s1),cbar(c.s1),ylabel('Strue');
im(222,se_shat_decomp,'shat by se-decomp',c.s1),cbar(c.s1)
xlab = sprintf('NRMS = %0.3g',Nrms(se_shat_decomp, se_strue_pd));xlabel(xlab);
im(223,se_shat_pwls,'shat by se-pwls',c.s1),cbar(c.s1)
xlab = sprintf('NRMS = %0.3g',Nrms(se_shat_pwls, se_strue_pd));xlabel(xlab);
im(224,se_shat_pl,'shat by se-pl',c.s1),cbar(c.s1)
xlab = sprintf('NRMS = %0.3g',Nrms(se_shat_pl, se_strue_pd));xlabel(xlab);

%Cost functions (SECT)
figure;tmp_leng = 20;
subplot(121);plot(se_cost_eval_pwls(2:tmp_leng)-max(se_cost_eval_pwls),'r-o','linewidth',2);grid on;
title('Cost function of SE-PWLS');xlabel('Iteration #');ylabel('\Phi(x^n)-\Phi(x^0)');
subplot(122);plot(se_cost_eval_pl(2:tmp_leng)-max(se_cost_eval_pl),'b-o','linewidth',2);grid on;
title('Cost function of SE-PL');xlabel('Iteration #');ylabel('\Psi(x^n)-\Psi(x^0)');
clear tmp_leng

%Individual terms of PL cost functions
se_cost_pl1 = se_cost_pl1 - se_cost_pl1(1);
se_cost_pl2 = se_cost_pl2 - se_cost_pl2(1);
%subplot(121);plot(se_cost_pl1,'r-o');grid on;
%subplot(122);plot(se_cost_pl2,'b-o');grid on;

%NRMS errors for each iteration (SECT)
leng_tmp = 20;
figure;plot((1:leng_tmp),se_nrms_pwls(1:leng_tmp),'b-o','linewidth',2);
hold on;plot((1:leng_tmp),se_nrms_pl(1:leng_tmp),'r-o','linewidth',2);
legend('PWLS','PL');title('Restored sinogram by SE');
xlabel('Number of iterations');ylabel('NRMS error');
clear leng_tmp

end % single energy


%%%%%%%%%%%%%%%%%%%%%%
% FBP reconstruction % to full size CT
%%%%%%%%%%%%%%%%%%%%%%
if 0
mask_d = logical(downsample2(mask,down));
ig_d = image_geom('nx', n.x/down, 'ny', n.y/down, 'dx', f.dx*down,...
			'mask', mask_d);
sgp_d = sino_geom('par', 'nb', n.bp, 'na', n.ap,'dr', dr*down, ...
		'offset_r', 0, 'orbit', 180);
fbp2tmp = fbp2(sgp_d, ig_d); %backward
G_d = Gtomo2_dscmex(sgp_d, ig_d); %forward
end

%DECT
if ~isvar('xfbp_decomp'), printm 'fbp'
	tmp = fbp2(sgsp, igs);
	xfbp_decomp = fbp2(shat_decomp, tmp);
	xfbp_decomp = max(xfbp_decomp, 0);
	im(xfbp_decomp)

	for ii=1:LL
	% xfbp_decomp(:,:,ii) = fbp2(shat_decomp(:,:,ii), fbp2tmp);
	% xfbp_pwls(:,:,ii) = fbp2(shat_pwls(:,:,ii), fbp2tmp);
	% xfbp_pl(:,:,ii) = fbp2(shat_pl(:,:,ii), fbp2tmp);
	end
	%xfbp_pwls = max(xfbp_pwls, 0);
	%xfbp_pl = max(xfbp_pl, 0);
end

%SECT
if 0
se_xfbp_decomp = fbp2(se_shat_decomp, fbp2tmp);
se_xfbp_decomp = max(se_xfbp_decomp, 0);
se_xfbp_pwls = fbp2(se_shat_pwls, fbp2tmp);
se_xfbp_pwls = max(se_xfbp_pwls, 0);
se_xfbp_pl = fbp2(se_shat_pl, fbp2tmp);
se_xfbp_pl = max(se_xfbp_pl,0);

figure;
im(331,xfbp_decomp(:,:,1));ylabel('Soft Tissue');xlabel('DE-Decomp');colorbar;
im(334,xfbp_decomp(:,:,2));ylabel('Bone');xlabel('DE-Decomp');colorbar;
im(332,xfbp_pwls(:,:,1));ylabel('Soft Tissue');xlabel('DE-PWLS');colorbar;
im(335,xfbp_pwls(:,:,2));ylabel('Bone');xlabel('DE-PWLS');colorbar;
im(333,xfbp_pl(:,:,1));ylabel('Soft Tissue');xlabel('DE-PL');colorbar;
im(336,xfbp_pl(:,:,2));ylabel('Bone');xlabel('DE-PL');colorbar;
im(337,se_xfbp_decomp);xlabel('SE-WLS');colorbar;
im(338,se_xfbp_pwls);xlabel('SE-PWLS');colorbar;
im(339,se_xfbp_pl);xlabel('SE-PL');colorbar;

figure;
im(se_xfbp_pl);colorbar;
figure;
im(xfbp_pl(:,:,2));colorbar;
end

%%%%%%%%%%%%%%%%%%
% ACF estimation % for DE-CT
%%%%%%%%%%%%%%%%%%
if ~isvar('acf_decomp'), printm 'acf'
	mac_511(1) = xray_read_atten('soft',511); %soft tissue
	mac_511(2) = xray_read_atten('bone',511); %bone
	mac_511(3) = xray_read_atten('iodine',511); %iodine
	acftrue = mac_511(1)*strue_pd(:,:,1) + mac_511(2)*strue_pd(:,:,2);
	if Iodine_sel == 1, acftrue = acftrue + mac_511(3)*strue_pd(:,:,3); end
	acf_decomp = mac_511(1)*shat_decomp(:,:,1) + mac_511(2)*shat_decomp(:,:,2);
	%acf_pwls = mac_511(1)*shat_pwls(:,:,1) + mac_511(2)*shat_pwls(:,:,2);
	%acf_pl = mac_511(1)*shat_pl(:,:,1) + mac_511(2)*shat_pl(:,:,2);

	im pl 1 2
	im(1, acftrue), im(2, acf_decomp)
end

if 0
%%%%%%%%%%%%%%%%%%%%
% Bilinear Scaling % for SECT
%%%%%%%%%%%%%%%%%%%%
rhocb = 1.919;
mucb = xray_read_atten('bone',se_ftab.xrs.eff)*rhocb; %LAC of c-bone
muw = xray_read_atten('water',se_ftab.xrs.eff); %LAC of water
muw_pet = xray_read_atten('water',511);
maccb_pet = xray_read_atten('bone',511);
mac_soft = xray_read_atten('soft',se_ftab.xrs.eff);

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

figure; %careful!!, taking exponential
im(241,(acftrue),'True ACF');colorbar horiz;
im(242,(acf_decomp),'ACF by DE-Decomp');colorbar horiz;
xlab = sprintf('NRMS = %0.3g',Nrms(exp(acf_decomp),exp(acftrue)));xlabel(xlab);
im(243,(acf_pwls),'ACF by DE-PWLS');colorbar horiz;
xlab = sprintf('NRMS = %0.3g',Nrms(exp(acf_pwls),exp(acftrue)));xlabel(xlab);
im(244,(acf_pl),'ACF by DE-PL');colorbar horiz;
xlab = sprintf('NRMS = %0.3g',Nrms(exp(acf_pl),exp(acftrue)));xlabel(xlab);
im(246,(se_acf(:,:,1)),'ACF by SE-WLS');colorbar horiz;
xlab = sprintf('NRMS = %0.3g',Nrms(exp(se_acf(:,:,1)),exp(acftrue)));xlabel(xlab);
im(247,(se_acf(:,:,2)),'ACF by SE-PWLS');colorbar horiz;
xlab = sprintf('NRMS = %0.3g',Nrms(exp(se_acf(:,:,2)),exp(acftrue)));xlabel(xlab);
im(248,(se_acf(:,:,3)),'ACF by SE-PL');colorbar horiz;
xlab = sprintf('NRMS = %0.3g',Nrms(exp(se_acf(:,:,3)),exp(acftrue)));xlabel(xlab);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PET/CT image reconstruction %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xtrue_tmp = xtrue(:,:,1);
xtrue_pet = xtrue_tmp + 0.1*xtrue(:,:,2); %original pet lambda
if Iodine_sel == 1
	xtrue3_tmp = ellipse_im(ig,ell1) + ellipse_im(ig,ell2); %removing 3rd iodine
	%xtrue3_tmp = ellipse_im(ig,ell1) + ellipse_im(ig,ell2) + ellipse_im(ig,ell3);
	xtrue_pet = xtrue_pet + xtrue3_tmp/max(col(xtrue3_tmp));%.*mean(col(xtrue_tmp(xtrue_tmp~=0)));
end
xtrue_petd = downsample2(xtrue_pet, f.down); %downsampled pet lambda

%PET sinogram
%temp = reshape(xtrue_pet, [], 1);
%strue_pet = reshape(G * temp(mask), [ f.down*n.bp down*n.ap]);
strue_pet = Gsp * xtrue_petd;

%strue_pet = downsample2(strue_pet,down);
%PET measurements (noisyless)
ymi_pet = strue_pet .*exp(-acftrue); % true attenuation

%Attenuation correction (DECT)
shat_pet_decomp = ymi_pet.*exp(acf_decomp);
%shat_pet_pwls = ymi_pet.*exp(acf_pwls);
%shat_pet_pl = ymi_pet.*exp(acf_pl);

return
% todo: done to here

%Attenuation correction (SECT by bilinear scaling)
se_shat_pet_decomp = ymi_pet.*exp(se_acf(:,:,1));
se_shat_pet_pwls = ymi_pet.*exp(se_acf(:,:,2));
se_shat_pet_pl = ymi_pet.*exp(se_acf(:,:,3));

%FBP reconstruction for PET/CT
fbp2tmp = fbp2(sgp_d, ig_d);
xfbp_pet_decomp = max(fbp2(shat_pet_decomp,fbp2tmp),0);
xfbp_pet_pwls = max(fbp2(shat_pet_pwls,fbp2tmp),0);
xfbp_pet_pl = max(fbp2(shat_pet_pl,fbp2tmp),0);
se_xfbp_pet_decomp = max(fbp2(se_shat_pet_decomp,fbp2tmp),0);
se_xfbp_pet_pwls = max(fbp2(se_shat_pet_pwls,fbp2tmp),0);
se_xfbp_pet_pl = max(fbp2(se_shat_pet_pl,fbp2tmp),0);

figure;
subplot(241);im(xtrue_petd,'True PET Image',c.soft);cbar([0 1],'h');
subplot(242);im(xfbp_pet_decomp,'PET Image with CTAC by DE Decomp',c.soft);cbar([0 1],'h');;
xlab = sprintf('NRMS = %0.2g',nrms(xfbp_pet_decomp,xtrue_petd));xlabel(xlab);
subplot(243);im(xfbp_pet_pwls,'PET Image with CTAC by DE PWLS',c.soft);cbar([0 1],'h');
xlab = sprintf('NRMS = %0.2g',nrms(xfbp_pet_pwls,xtrue_petd));xlabel(xlab);
subplot(244);im(xfbp_pet_pl,'PET Image with CTAC by DE PL',c.soft);cbar([0 1],'h');
xlab = sprintf('NRMS = %0.2g',nrms(xfbp_pet_pl,xtrue_petd));xlabel(xlab);
subplot(246);im(se_xfbp_pet_decomp,'PET Image with CTAC by SE WLS',c.soft);cbar([0 1],'h');;
xlab = sprintf('NRMS = %0.3g',nrms(se_xfbp_pet_decomp,xtrue_petd));xlabel(xlab);
subplot(247);im(se_xfbp_pet_pwls,'PET Image with CTAC by SE PWLS',c.soft);cbar([0 1],'h');
xlab = sprintf('NRMS = %0.3g',nrms(se_xfbp_pet_pwls,xtrue_petd));xlabel(xlab);
subplot(248);im(se_xfbp_pet_pl,'PET Image with CTAC by SE PL',c.soft);cbar([0 1],'h');
xlab = sprintf('NRMS = %0.3g',nrms(se_xfbp_pet_pl,xtrue_petd));xlabel(xlab);

%For TMI - no Iodine
figure;im(xtrue_petd,'True PET Image',c.soft);cbar([0 1]);
figure;im(xfbp_pet_decomp,'PET Image with CTAC by DE Decomp',c.soft);cbar([0 1]);
figure;im(xfbp_pet_pwls,'PET Image with CTAC by DE PWLS',c.soft);cbar([0 1]);
figure;im(xfbp_pet_pl,'PET Image with CTAC by DE PL',c.soft);cbar([0 1]);

%For TMI - with Iodine
figure;im(xtrue_petd,'True PET Image',c.soft);cbar([0 1]);
figure;im(xfbp_pet_decomp,'PET Image with CTAC by DE Decomp',c.soft);cbar([0 1]);
figure;im(xfbp_pet_pwls,'PET Image with CTAC by DE PWLS',c.soft);cbar([0 1]);
figure;im(xfbp_pet_pl,'PET Image with CTAC by DE PL',c.soft);cbar([0 1]);
figure;im(se_xfbp_pet_decomp,'PET Image with CTAC by SE WLS',c.soft);cbar([0 1]);
figure;im(se_xfbp_pet_pwls,'PET Image with CTAC by SE PWLS',c.soft);cbar([0 1]);
figure;im(se_xfbp_pet_pl,'PET Image with CTAC by SE PL',c.soft);cbar([0 1]);


%%%%%%%%%%%%%%%%%%%%
% Local RMS Errors % ROI selection issue
%%%%%%%%%%%%%%%%%%%%
pet_mask = zeros(size(xfbp_pet_decomp));
pet_mask((30:45),(50:70)) = 1; %region 1
pet_mask((61:70),(49:57)) = 1; %region 2
pet_mask((55:75),(60:80)) = 1; %region 3
pet_mask((50:80),(87:93)) = 1; %region 4

% pet_mask((35:44),(61:70)) = 1;
% pet_mask((61:70),(61:70)) = 1;
% pet_mask((55:74),(89:93)) = 1;

%masked true pet
mxtrue_petd = xtrue_petd.*pet_mask;
mxtrue_petd_comp = xtrue_petd.*(1-pet_mask);
%masked DE-decomp and SE-decomp
mxfbp_pet_decomp(:,:,1) = xfbp_pet_decomp.*pet_mask;
mxfbp_pet_decomp(:,:,2) = se_xfbp_pet_decomp.*pet_mask;
%masked DE-PWLS and SE-PWLS
mxfbp_pet_pwls(:,:,1) = xfbp_pet_pwls.*pet_mask;
mxfbp_pet_pwls(:,:,2) = se_xfbp_pet_pwls.*pet_mask;
%masked DE-PL and SE-PL
mxfbp_pet_pl(:,:,1) = xfbp_pet_pl.*pet_mask;
mxfbp_pet_pl(:,:,2) = se_xfbp_pet_pl.*pet_mask;

%nrms errors
nrms_mxfbp_pet_decomp = [nrms(mxfbp_pet_decomp(:,:,1),mxtrue_petd), nrms(mxfbp_pet_decomp(:,:,2),mxtrue_petd)]
nrms_mxfbp_pet_pwls = [nrms(mxfbp_pet_pwls(:,:,1),mxtrue_petd), nrms(mxfbp_pet_pwls(:,:,2),mxtrue_petd)]
nrms_mxfbp_pet_pl = [nrms(mxfbp_pet_pl(:,:,1),mxtrue_petd), nrms(mxfbp_pet_pl(:,:,2),mxtrue_petd)]

%masked region plot with red boxes
petmask_bound = pet_mask;
petmask_bound((31:44),(51:69)) = 0;
petmask_bound((62:69),(50:56)) = 0;
petmask_bound((56:74),(61:79)) = 0;
petmask_bound((51:79),(88:92)) = 0;
petmask_boundind = find(petmask_bound > 0);

mxtrue_pet_bound_r = xtrue_petd/max(max(xtrue_petd)); %red
mxtrue_pet_bound_g = xtrue_petd/max(max(xtrue_petd)); %green
mxtrue_pet_bound_b = xtrue_petd/max(max(xtrue_petd)); %blue

mxtrue_pet_bound_r(petmask_boundind) = 1;mxtrue_pet_bound(:,:,1) = mxtrue_pet_bound_r;
mxtrue_pet_bound_g(petmask_boundind) = 0;mxtrue_pet_bound(:,:,2) = mxtrue_pet_bound_g;
mxtrue_pet_bound_b(petmask_boundind) = 0;mxtrue_pet_bound(:,:,3) = mxtrue_pet_bound_b;

mxtrue_pet_bound = permute(mxtrue_pet_bound, [2 1 3]);
imagesc(mxtrue_pet_bound);axis square;axis off; %red boxes
mxtrue_pet_bound1 = permute(mxtrue_pet_bound, [2 1 3]);
im(mxtrue_pet_bound1(:,:,1)) %white boxes

%%%%%%%%%%%%%%%%%%%%%%
% Horizontal profile % for bias checking between DECT and SECT ?? (Ask Fessler)
%%%%%%%%%%%%%%%%%%%%%%
if isinf(f.scale_high) || isinf(f.scale_low)
	hori_index = 64;%52
	%DE-PL vs SE-PL
	figure(1);subplot(121);imagesc(xtrue_petd');axis square;colormap gray;colorbar;hold on;
	plot(hori_index*ones(size(xtrue_petd,1)),'r--');hold off;title('True PET Image');
	xtrue_pet_profile = xtrue_petd(:,hori_index);
	xfbp_pet_profile_pl = xfbp_pet_pl(:,hori_index);
	se_xfbp_pet_profile_pl = se_xfbp_pet_pl(:,hori_index);
	subplot(122);plot(xtrue_pet_profile,'c','linewidth',2);hold on;
	plot(xfbp_pet_profile_pl,'r','linewidth',2);hold on;
	plot(se_xfbp_pet_profile_pl,'b-.','linewidth',2);hold off;
	legend('True','DE-PL','SE-PL','location','northwest');
	title('PET Noiseless Horizontal profile');
	%DE-Decomp vs DE-PL
	figure(2);subplot(121);imagesc(xtrue_petd');axis square;colormap gray;colorbar;hold on;
	plot(hori_index*ones(size(xtrue_petd,1)),'r--');hold off;title('True PET Image');
	xtrue_pet_profile = xtrue_petd(:,hori_index);
	xfbp_pet_profile_decomp = xfbp_pet_decomp(:,hori_index);
	xfbp_pet_profile_pl = xfbp_pet_pl(:,hori_index);
	subplot(122);figure;plot(xtrue_pet_profile,'c','linewidth',2);hold on;
	plot(xfbp_pet_profile_pl,'r','linewidth',2);hold on;
	plot(xfbp_pet_profile_decomp,'b-.','linewidth',2);hold off;
	legend('True','DE-PL','DE-Decomp','location','northwest');
	title('PET Noiseless Horizontal profile');
end

% Data save
filename = sprintf('sino_restore_results512_R%d_I%d.mat',reg_sel,Iodine_sel);
%save(filename);
%load sino_restore_results512_R1_I1.mat;
