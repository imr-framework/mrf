%This is the code to test the designed penalizing regularizers 
%for spatially uniform resolution in DECT sinogram restoration.
% - modified from de_sino_restore.m

clear all;
load xtrue.mat;
[n.x n.y] = size(xtrue(:,:,1));
s_factor = 1;
n.b = 140*s_factor; n.a = 128*s_factor;f.dx = 0.16 * 2;
f.stype = 'ps1'; %setting dual spectra
ftab = de_ftab_build({f.stype});
%mask
mask = sum(xtrue,3) > 0;
mask = imdilate(mask, strel('disk', 5));
mask = logical(mask);

%%%%%%%%%%%%%%%%
%System matrix %
%%%%%%%%%%%%%%%%
ig = image_geom('nx', n.x, 'ny', n.y, 'dx', f.dx, 'mask', mask);
sg = sino_geom('par', 'nb', n.b, 'na', n.a, 'dr', f.dx);
G = Gtomo2_strip(sg, ig);   %G = Gtomo2_table(sg, ig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sinogram W/ and W/O perturbation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = reshape(xtrue,n.x*n.y,2);
%%%%%%%%%%%%%%%%%%%% 
%Original sinogram %
%%%%%%%%%%%%%%%%%%%%
strue = reshape(G * tmp(mask,:), [n.b n.a 2]);
% im(121,strue(:,:,1));cbar;
% im(122,strue(:,:,2));cbar;
%f
ftrue = ftab.fm_fun(ftab,{strue(:,:,1),strue(:,:,2)});
%noiseless measurement
ybi = zeros(n.b,n.a,2); %bar{y}_{mi}
ybi(:,:,1) = ftab.xray.I(1) * exp(-ftrue(:,:,1)); %true mean
ybi(:,:,2) = ftab.xray.I(2) * exp(-ftrue(:,:,2)); %assuming no background contribution
ymi = ybi;	% noiseless
%%%%%%%%%%%%%%%%%%%%%
%Perturbed sinogram %
%%%%%%%%%%%%%%%%%%%%%
Nd = size(strue);pert = zeros(Nd);
Nd = Nd(1:end-1);
delta = 0.001; %perturbation amplitude
locx = [(2:5)]/7;locy = [(1:7)]/8;
[locxx locyy] = ndgrid(locx,locy);
%assume each perturbation are far away from each other
pert(Nd(1)*locxx+1,Nd(2)*locyy+1,1) = delta; %only in 1st compo sino 
pert_strue = strue + pert;
%perturbed f
pert_ftrue = ftab.fm_fun(ftab,{pert_strue(:,:,1),pert_strue(:,:,2)});
%noiseless measurement
pert_ybi = zeros(n.b,n.a,2);
pert_ybi(:,:,1) = ftab.xray.I(1) * exp(-pert_ftrue(:,:,1));
pert_ybi(:,:,2) = ftab.xray.I(2) * exp(-pert_ftrue(:,:,2));
pert_ymi = pert_ybi;

%%%%%%%%%%%%%%%%%%%%%%%%
% Sinogram Restoration % 
%%%%%%%%%%%%%%%%%%%%%%%%
%penalty selection
reg_sel = 2; %0|1|2 (conventional/modified/hybrid)
if reg_sel == 1
    sniter = 30;regbeta = [2^(1) 2^(1)]; %
elseif reg_sel == 0
    sniter = 30;regbeta = [2^(1) 2^(1)]; %quick mode for check
elseif reg_sel == 2
    sniter = 30;regbeta = [2^(-5) 2^(-5)];
    %sniter = 90;regbeta = [0 0]; %ML-post smoothing
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Restoration for original sinogram %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fhat.raw(:,:,1) = -log(ymi(:,:,1) / ftab.xray.I(1)); 
fhat.raw(:,:,2) = -log(ymi(:,:,2) / ftab.xray.I(2));

%conventional decomposition
[shat_decomp cost_decomp] = de_ftab_s_iter(ftab.fit,fhat.raw,'data',ymi,'niter',10);  %by (W)LS
shat_decomp = reshape(shat_decomp,[n.b n.a 2]);
%Post smoothing
f.fbp_kernel = [1 2 1]';
f.fbp_kernel = f.fbp_kernel / sum(f.fbp_kernel);
shat_decomp(:,:,1) = conv2(shat_decomp(:,:,1), f.fbp_kernel, 'same'); 
shat_decomp(:,:,2) = conv2(shat_decomp(:,:,2), f.fbp_kernel, 'same');
shat_decomp = reshape(shat_decomp, [], 2); %to initialize PWLS and PL

%statistical restoration
[shat_pwls cost_eval_pwls nrms_pwls] = de_ftab_s_pwls(ftab.fit,fhat.raw,ymi,strue,'init',shat_decomp,...
                                                      'niter',sniter,'beta',regbeta,'regsel',reg_sel); 
[shat_pl cost_eval_pl nrms_pl] = de_ftab_s_pl(ftab.fit,ftab.xray,ftab.mac.mac,fhat.raw,ymi,strue,'init',shat_decomp,...
                                              'niter',sniter,'beta',regbeta,'curvtype','pc','regsel',reg_sel);

shat_decomp = reshape(shat_decomp,[n.b n.a 2]);
shat_pwls = reshape(shat_pwls,[n.b n.a 2]); 
shat_pl = reshape(shat_pl,[n.b n.a 2]); 
if reg_sel == 2
    post_kernel = [0.8 2 0.8]';
    post_kernel = post_kernel / sum(post_kernel);
    shat_pl(:,:,1) = conv2(shat_pl(:,:,1), post_kernel, 'same');
    shat_pl(:,:,2) = conv2(shat_pl(:,:,2), post_kernel, 'same');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Restoration for perturbed sinogram %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pert_fhat.raw(:,:,1) = -log(pert_ymi(:,:,1) / ftab.xray.I(1)); 
pert_fhat.raw(:,:,2) = -log(pert_ymi(:,:,2) / ftab.xray.I(2));

%conventional decomposition
[pert_shat_decomp pert_cost_decomp] = de_ftab_s_iter(ftab.fit,pert_fhat.raw,'data',pert_ymi,'niter',10);  %by (W)LS
pert_shat_decomp = reshape(pert_shat_decomp,[n.b n.a 2]);
%Post smoothing
pert_shat_decomp(:,:,1) = conv2(pert_shat_decomp(:,:,1), f.fbp_kernel, 'same'); 
pert_shat_decomp(:,:,2) = conv2(pert_shat_decomp(:,:,2), f.fbp_kernel, 'same');
pert_shat_decomp = reshape(pert_shat_decomp, [], 2); %to initialize PWLS and PL

%statistical restoration
[pert_shat_pwls pert_cost_eval_pwls pert_nrms_pwls] = de_ftab_s_pwls(ftab.fit,pert_fhat.raw,pert_ymi,pert_strue,'init',pert_shat_decomp,...
                                                      'niter',sniter,'beta',regbeta,'regsel',reg_sel); 
[pert_shat_pl pert_cost_eval_pl pert_nrms_pl] = de_ftab_s_pl(ftab.fit,ftab.xray,ftab.mac.mac,pert_fhat.raw,pert_ymi,pert_strue,'init',pert_shat_decomp,...
                                                'niter',sniter,'beta',regbeta,'curvtype','pc','regsel',reg_sel);

pert_shat_decomp = reshape(pert_shat_decomp,[n.b n.a 2]);
pert_shat_pwls = reshape(pert_shat_pwls,[n.b n.a 2]); 
pert_shat_pl = reshape(pert_shat_pl,[n.b n.a 2]); 
if reg_sel == 2
    pert_shat_pl(:,:,1) = conv2(pert_shat_pl(:,:,1), post_kernel, 'same');
    pert_shat_pl(:,:,2) = conv2(pert_shat_pl(:,:,2), post_kernel, 'same');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linearized local impulse response %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
llir_decomp = (pert_shat_decomp-shat_decomp)/delta;
llir_pwls = (pert_shat_pwls-shat_pwls)/delta;
llir_pl = (pert_shat_pl-shat_pl)/delta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contour plot of Linearized LIR %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;subplot(121);contour(llir_decomp(:,:,1),5);grid on;hold on;
plot(Nd(2)*locyy+1,flipud(Nd(1)*locxx+1),'rx');colorbar;
xlabel('Angle'):ylabel('Radius');title('Linearized LIR for soft (Decomp)');
subplot(122);contour(llir_decomp(:,:,2),5);grid on;hold on;
plot(Nd(2)*locyy+1,flipud(Nd(1)*locxx+1),'rx');colorbar;
xlabel('Angle'):ylabel('Radius');title('Linearized LIR for Bone (Decomp)');

llir_decomp1 = llir_decomp(:,:,1);figure;
xind = Nd(1)*locxx+1;xind = xind(:,1)';
yind = Nd(2)*locyy+1;yind = yind(1,:);

for ii=1:length(locx)
    for jj=1:length(locy)
        xx = (xind(ii)+(-5:5));
        %subplot(330+3*(ii-1)+jj);
        subplot(length(locx),length(locy),length(locy)*(ii-1)+jj);
        plot(-5:5,llir_decomp1(xx,yind(jj)),'r-o');grid on;
        ylim([min(llir_decomp1(:)) max(llir_decomp1(:))]);
        %tit = sprintf('Angle location = %d (Decomp)',yind(jj));
        %xlabel('Radial pixel location');ylabel('Linearized LIR');
        %title(tit);
    end
end

figure;subplot(121);contour(llir_pwls(:,:,1),10);grid on;hold on;
plot(Nd(2)*locyy+1,flipud(Nd(1)*locxx+1),'rx');colorbar;
xlabel('Angle'):ylabel('Radius');title('Linearized LIR for soft (PWLS)');
subplot(122);contour(llir_pwls(:,:,2),10);grid on;hold on;
plot(Nd(2)*locyy+1,flipud(Nd(1)*locxx+1),'rx');colorbar;
xlabel('Angle'):ylabel('Radius');title('Linearized LIR for bone (PWLS)');

llir_pwls1 = llir_pwls(:,:,1);figure;
for ii=1:length(locx)
    for jj=1:length(locy)
        xx = (xind(ii)+(-5:5));
        %subplot(330+3*(ii-1)+jj);
        subplot(length(locx),length(locy),length(locy)*(ii-1)+jj);
        plot(-5:5,llir_pwls1(xx,yind(jj)),'r-o');grid on;
        ylim([min(llir_pwls1(:)) max(llir_pwls1(:))]);
        %tit = sprintf('Angle location = %d (PWLS)',yind(jj));
        %xlabel('Radial pixel location');ylabel('Linearized LIR');
        %title(tit);
    end
end

figure;subplot(121);contour(llir_pl(:,:,1),10);grid on;hold on;
plot(Nd(2)*locyy+1,flipud(Nd(1)*locxx+1),'rx');colorbar;
xlabel('Angle'):ylabel('Radius');title('Linearized LIR for soft (PL)');
subplot(122);contour(llir_pl(:,:,2),10);grid on;hold on;
plot(Nd(2)*locyy+1,flipud(Nd(1)*locxx+1),'rx');colorbar;
xlabel('Angle'):ylabel('Radius');title('Linearized LIR for bone n (PL)');

llir_pl1 = llir_pl(:,:,1);figure;
for ii=1:length(locx)
    for jj=1:length(locy)
        xx = (xind(ii)+(-5:5));
        %subplot(330+3*(ii-1)+jj);
        subplot(length(locx),length(locy),length(locy)*(ii-1)+jj);
        plot(-5:5,llir_pl1(xx,yind(jj)),'r-o');grid on;
        ylim([min(llir_pl1(:)) max(llir_pl1(:))]);
        %tit = sprintf('(r,phi) =  (%d,%d) (PL)',xind(ii),yind(jj));
        %xlabel('Radial pixel location');ylabel('Linearized LIR');
        %title(tit);
    end
end

%Locations of impulse in the sinogram
figure;imagesc(strue(:,:,1));hold on;colormap gray;
xlabel('Angle');ylabel('Radius');title('Impulse response locations');
plot(Nd(2)*locyy+1,Nd(1)*locxx+1,'rx');colorbar;
%
figure;subplot(121);mesh(strue(:,:,1));
xlabel('Angle');ylabel('Radius');title('True sinogram');
subplot(122);plot(strue(:,65,1),'r-o')
xlabel('Radius');ylabel('Amplitude');grid on;
title('A profile of true sinogram at angle = 65');






