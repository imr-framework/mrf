%This is the code for DE-CT sinogram restoration
%Modified from spie02.m

clear all;
im off;
%%%%%%%%%%%%
%Load data %
%%%%%%%%%%%%
load xtrue.mat;
[n.x n.y] = size(xtrue(:,:,1));
n.b = 140; n.a = 128;f.dx = 0.16 * 2;
f.stype = 'ps1';
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
G = Gtomo2_strip(sg, ig);   %G = Gtomo2_table(sg, ig);
gi = reshape(sum(G'), [n.b n.a]); %simple projection
im clf, im(121, gi, 'gi'), im(122, gi>0, 'gi>0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ideal and Noisy Measurement % y_{mi}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Noiseless sinogram and f
%
tmp = reshape(xtrue, n.x*n.y, 2);
tic %noiseless forward projection
strue = reshape(G * tmp(mask,:), [n.b n.a 2]);
printm('forward projection time %g', toc)
%sinogram display
im clf, im(421, strue(:,:,1), 'noiseless s1'), cbar, im(422, strue(:,:,2), 'noiseless s2'), cbar
s1max = max(col(strue(:,:,1)));
s2max = max(col(strue(:,:,2)));

ftrue = ftab.fm_fun(ftab, {strue(:,:,1), strue(:,:,2)}); %noiseless f
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
end

f.title1 = sprintf('%2.0f kVp', ftab.xray.kvp(1));
f.title2 = sprintf('%2.0f kVp', ftab.xray.kvp(2));
if isinf(f.scale)
	f.title1 = [f.title1 ' (Noiseless)'];
	f.title2 = [f.title2 ' (Noiseless)'];
else
	printm('total counts %g', sum(ymi(:)) * f.scale)
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
printm('fhat.raw err %g%% %g%%', fhat.err.raw(1), fhat.err.raw(2))

%%%%%%%%%%%%%%%%%%%%%%%%
% Sinogram Restoration %
%%%%%%%%%%%%%%%%%%%%%%%%
%clear all;load de_ct_setup.mat; %high dose
%clear all;load de_ct_low4_setup.mat; %low dose : f.scale = 5*10^4 / ybi(1,1,2); 
%clear all;load de_ct_setup_noiseless.mat; %noiseless setup
sino_sel = 3; %1 for decomp, 2 for PWLS, 3 for PL
if sino_sel == 1
    shat = de_ftab_s_iter(ftab.fit, fhat.raw, 'niter', 300);  %by LS
elseif sino_sel == 2
    [shat cost_eval] = de_ftab_s_pwls(ftab.fit, fhat.raw, ymi,'niter', 300, 'beta', 2^(-5));  %by PWLS
    %[shat cost_eval] = de_ftab_s_pwls2(ftab.fit, fhat.raw, ymi,'niter', 300, 'beta', 2^(-7));  %by PWLS
    im clf;plot((0:20),cost_eval(1:21)-max(cost_eval),'r-o','linewidth',2);grid on;
    title('Cost function of PWLS');xlabel('Iteration n');ylabel('\Psi(x^n)-\Psi(x^0)');
elseif sino_sel == 3
    regbeta.soft = 2^(-5); regbeta.bone = regbeta.bone; %one beta for regularization
    shat_pwls = de_ftab_s_pwls(ftab.fit, fhat.raw, ymi,'niter',40,'beta',regbeta.soft);  %by PWLS
    [shat cost_eval] = de_ftab_s_pl(ftab.fit,ftab.xray, ftab.mac.mac,fhat.raw,ymi,...
                                    'init',shat_pwls,'niter',300,'beta',regbeta,'curvtype','pc'); %by PL
    im clf;plot((0:20),cost_eval(1:21)-max(cost_eval),'r-o','linewidth',2);grid on;
    title('Cost function of PL');xlabel('Iteration n');ylabel('\Psi(x^n)-\Psi(x^0)');
end
shat = reshape(shat,[n.b n.a 2]); 
shat_err = shat - strue;
c.s1 = [0 30];c.s2 = [0 20];c.s_err = [0 10];

im clf
im(231,strue(:,:,1),'strue1',c.s1),cbar(c.s1),ylabel('Soft Tissue');
im(232,shat(:,:,1),'shat1',c.s1), cbar(c.s1)
im(233,abs(shat_err(:,:,1)),'|shat err1|',c.s_err),cbar(c.s_err)
im(234,strue(:,:,2),'strue2',c.s2),cbar(c.s2),ylabel('Bone');
im(235,shat(:,:,2),'shat2',c.s2), cbar(c.s2)
im(236,abs(shat_err(:,:,2)),'|shat err2|',c.s_err),cbar(c.s_err)
printm('Sinogram restoration %g %g', ...
    Nrms(shat(:,:,1), strue(:,:,1)), Nrms(shat(:,:,2), strue(:,:,2)))
%s_iter Nrms : 0.074 / 0.19 (0.3162 / 0.8306  - low)

%s_pwls Nrms : 0.056 / 0.15, beta = 2^(-5) (0.1336 / 0.3445 - low)
%s_pwls Nrms : 0.048 / 0.12, beta = 2^(-7) 
%s_pwls Nrms : beta = 2^(-3) (0.1141/0.3000 - low)
%s_pwls Nrms : beta = 2^(-4) (0.1202/0.3123 - low)
%s_pwls Nrms : beta = 2^(-6) (0.1547/0.3981 - low)

%s_pl Nrms : 0.053 / 0.141, beta = 2^(-6) (0.1397/0.3592 - low)
%s_pl Nrms : 0.058 / 0.155, beta = 2^(-5) (0.1234/0.3182 - low)
%s_pl Nrms : beta = 2^(-4) (0.1144/0.2970 - low)
%s_pl Nrms : beta = 2^(-3) (0.1119 /0.2932 - low)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FBP reconstruction from pre-processed data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = fbp2(sg, ig); %for setup
if sino_sel == 1
    fbp_title = ['FBP [1 2 1]/4'];
    %f.fbp_kernel = 1;
    f.fbp_kernel = [1 2 1]';   %smoothing for the conventional decomp
else
    fbp_title = ['FBP'];
    f.fbp_kernel = 1; %no smoothing for pwls & pl
end 

f.fbp_kernel = f.fbp_kernel / sum(f.fbp_kernel);
sino1 = conv2(shat(:,:,1), f.fbp_kernel, 'same'); 
sino2 = conv2(shat(:,:,2), f.fbp_kernel, 'same');
xfbp(:,:,1) = fbp2(sino1, tmp);
xfbp(:,:,2) = fbp2(sino2, tmp);
xfbp = max(xfbp, 0); %negative out

im clf
im(231,xtrue(:,:,1),'Truth',c.soft),cbar([0 1]),ylabel('Soft Tissue');
im(234,xtrue(:,:,2),'',c.bone),cbar([0 2]),ylabel('Bone');
im(232,xfbp(:,:,1),'xhat1',c.soft),cbar([0 1])
im(235,xfbp(:,:,2),'xhat2',c.bone),cbar([0 2])
im(233,abs(xtrue(:,:,1)-xfbp(:,:,1)),'|Error|'),cbar
im(236,abs(xtrue(:,:,2)-xfbp(:,:,2))),cbar

%horizontal profile
if isinf(f.scale)
    hori_index = 52;%52
    im clf
    figure(1);subplot(1,2,1);imagesc(xtrue(:,:,1)'+xtrue(:,:,2)');colormap gray;colorbar;hold on;
    plot(hori_index*ones(size(xtrue,1)),'r--');hold off;title('True density');
    xtrue_profile = xtrue(:,hori_index,1)+xtrue(:,hori_index,2);
    xfbp_profile = xfbp(:,hori_index,1)+xfbp(:,hori_index,2);
    subplot(122);plot(xtrue_profile,'b','linewidth',2);hold on;
    plot(xfbp_profile,'r','linewidth',2);hold off
    if sino_sel == 1, legtit = ['Decomp-FBP'];  
    elseif sino_sel == 2, legtit = ['PWLS-FBP'];
    elseif sino_sel == 3, legtit = ['PL-FBP'];  end
    legend('True',legtit,'location','northwest');title('Noiseless Horizontal profile');
end

printm('FBP %g %g', Nrms(xfbp(:,:,1), xtrue(:,:,1)), Nrms(xfbp(:,:,2), xtrue(:,:,2)))
%FBP_s_iter Nrms : 0.187 / 0.313, (0.5408 / 0.6366 - low)

%FBP_s_pwls Nrms : 0.234 / 0.347, beta = 2^(-5) (0.3265 / 0.4137 - low)
%FBP_s_pwls Nrms : beta = 2^(-3) (0.3315/0.4573 - low)
%FBP_s_pwls Nrms : beta = 2^(-4) (0.3141/0.4242 - low)
%FBP_s_pwls Nrms : 0.204 / 0.307, beta = 2^(-6) (0.3917 / 0.4507 - low)
%FBP_s_pwls Nrms : 0.191 / 0.279, beta = 2^(-7) 

%FBP_s_pl Nrms : beta = 2^(-3) (0.3448/0.4705 - low)
%FBP_s_pl Nrms : beta = 2^(-4) (0.3203/0.4346 - low)
%FBP_s_pl Nrms : beta = 2^(-5) (0.3131/0.4106 - low)
%FBP_s_pl Nrms : 0.2162 / 0.3231, beta = 2^(-6) (0.3473/0.4205 - low)
%FBP_s_pl Nrms : 0.1959 / 0.2908, beta = 2^(-7) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attenuation Correction Factor %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mac at 511 kev for PET
mac_soft = xray_read_atten('soft',511);
mac_bone = xray_read_atten('bone',511);
%acfs
acf_true = mac_soft*strue(:,:,1) + mac_bone*strue(:,:,2);
%acf_fbp = mac_soft*shat(:,:,1) + mac_bone*shat(:,:,2); %?? Ask Fessler
acf_fbp = mac_soft*(G*xfbp(:,:,1)) + mac_bone*(G*xfbp(:,:,2));
if sino_sel==1, tit=sprintf('Estimated Atten Correc Fact by Decomp'); end
if sino_sel==2, tit=sprintf('Estimated Atten Correc Fact by PWLS');   end
if sino_sel==3, tit=sprintf('Estimated Atten Correc Fact by PL');     end
im(121,acf_true,'True Atten Correc Fact');cbar;
im(122,acf_fbp,tit);cbar;
printm('NRMS of ACF %g', Nrms(acf_fbp,acf_true))
%pl : 0.0286 / pwls : 0.0307 / decomp : 0.0655
%pl : 0.0455 / pwls : 0.0710 / decomp : 0.1914












