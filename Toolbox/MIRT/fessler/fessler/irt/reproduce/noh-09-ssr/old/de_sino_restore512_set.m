%This is the code for DE-CT sinogram restoration
%Modified from spie02.m

clear all;
im off;
%%%%%%%%%%%%
%Load data %
%%%%%%%%%%%%
load xtrue512.mat; %512 by 512
[n.x n.y] = size(xtrue(:,:,1));
f.stype = 'ps1'; %setting dual-kVp spectra
ftab = de_ftab_build({f.stype});
%mask
mask = sum(xtrue, 3) > 0;
mask = imdilate(mask, strel('disk', 5));
mask = logical(mask);
%true density map display
im clf
c.soft = [0 1.2];c.bone = [0 2.2];c.dens = [0 2.2];
im(221, xtrue(:,:,1), 'SoftTissue Density', c.soft), cbar([0 1])
im(222, xtrue(:,:,2), 'Bone Density', c.bone), cbar([0 2])
im(223, sum(xtrue,3), 'Density Map', c.dens), cbar([0 2])
im(224, (sum(xtrue,3)>0.5) + 20*double(mask), 'Reconstruction Support')

%%%%%%%%%%%%%%%%
%System matrix %
%%%%%%%%%%%%%%%%
%Geometry selection
geo_sel = 1; %0|1 (para|fan-beam)
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
    G = Gtomo2_dscmex(sgf, ig); %full CT (Para-beam)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ideal and Noisy Measurement % y_{mi}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Noiseless sinogram and f
%
tmp = reshape(xtrue, n.x*n.y, 2);
tic; %noiseless forward projection
strue = reshape(G * tmp(mask,:), [n.bf n.af 2]);
printm('forward projection time %0.3g', toc)
ftrue = ftab.fm_fun(ftab, {strue(:,:,1), strue(:,:,2)}); %noiseless f

%Noisy measurement
ybi = zeros(n.bf, n.af, 2); %bar{y}_{mi}
ybi(:,:,1) = ftab.xray.I(1) * exp(-ftrue(:,:,1)); %true mean
ybi(:,:,2) = ftab.xray.I(2) * exp(-ftrue(:,:,2));
%f.scale = inf; %for noiseless profile
%f.scale = 10^6 / ybi(1,1,2); %high dose
f.scale = 5*10^4 / ybi(1,1,2); %low dose
if isinf(f.scale)
	ymi = ybi;	% noiseless
else                                         %%%%%%%%%%%%%%%%%%%%%%%
	ymi = poisson(f.scale*ybi, 0) / f.scale; %corrupted by poisson %
end %measurements look like sinogram         %%%%%%%%%%%%%%%%%%%%%%%

% borrow neighbor(s) for any log(0) values
off = 1;
while any(ymi(:) == 0)
	ii = find(ymi(:) == 0);
	printm('fixing %d zeros in ymi with %d', length(ii), off)
	ymi(ii) = max(ymi(ii+off), ymi(ii-off));
	off = off + 1;
end, clear off

%%%%%%%%%%%%%%%
% Geometry Tx % & downsampling by averaging
%%%%%%%%%%%%%%%
down = 4;
if geo_sel == 1
    ymipd(:,:,1) = downsample2(rebin_fan2par(ymi(:,:,1), sgf, sgp),down); %fan -> par
    ymipd(:,:,2) = downsample2(rebin_fan2par(ymi(:,:,2), sgf, sgp),down); %[n.bp n.ap]
    strue_pd(:,:,1) = downsample2(rebin_fan2par(strue(:,:,1), sgf, sgp),down);
    strue_pd(:,:,2) = downsample2(rebin_fan2par(strue(:,:,2), sgf, sgp),down); 
elseif geo_sel == 0
    ymipd(:,:,1) = downsample2(ymi(:,:,1),down);  %downsampling
    ymipd(:,:,2) = downsample2(ymi(:,:,2),down); %[n.bp n.ap]
    strue_pd(:,:,1) = downsample2(strue(:,:,1),down);
    strue_pd(:,:,2) = downsample2(strue(:,:,2),down); 
end
n.bp = n.bp/down;n.ap = n.ap/down; %down-sampling by 4
%n.bp = 1024;n.ap = 800; (original)

save de_ct_low512_setup_q.mat;

