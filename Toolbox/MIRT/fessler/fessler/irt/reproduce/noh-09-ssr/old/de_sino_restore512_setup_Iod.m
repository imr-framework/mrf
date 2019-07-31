%This is the code for DE-CT sinogram restoration
%Modified from spie02.m

%%%%%%%%%%%%
%Load data %
%%%%%%%%%%%%
load xtrue512.mat; %512 by 512
    rhoI = 4.930; %Iodine density (from Hubbell&Seltzer)
    xrang = (241:260);yrang = xrang;
    xtrue(xrang,yrang,1) = 0;
    xtrue(xrang,yrang,2) = 0;
    xtrue(xrang,yrang,3) = rhoI;
[n.x n.y] = size(xtrue(:,:,1));
f.stype = 'ps1'; %setting dual-kVp spectra
ftab = de_ftab_build({f.stype});
%mask
mask = sum(xtrue, 3) > 0;
mask = imdilate(mask, strel('disk', 5));
mask = logical(mask);

%%%%%%%%%%%%%%%%
%System matrix %
%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ideal and Noisy Measurement % y_{mi}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Noiseless sinogram and f
%
tmp = reshape(xtrue, n.x*n.y, 3);
tic; %noiseless forward projection
strue = reshape(G * tmp(mask,:), [n.bf n.af 3]);
printm('forward projection time %0.3g', toc)
ftrue = ftab.fm_fun(ftab, {strue(:,:,1),strue(:,:,2),strue(:,:,3)}); %noiseless f

%Noisy measurement
ybi = zeros(n.bf, n.af, 2); %bar{y}_{mi}
ybi(:,:,1) = ftab.xray.I(1) * exp(-ftrue(:,:,1)); %true mean
ybi(:,:,2) = ftab.xray.I(2) * exp(-ftrue(:,:,2));

%Dose amount
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry Tx & downsampling %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
down = 4;
if geo_sel == 1 %fan-beam
    ymipd(:,:,1) = downsample2(rebin_fan2par(ymi(:,:,1), sgf, sgp),down); %fan -> par
    ymipd(:,:,2) = downsample2(rebin_fan2par(ymi(:,:,2), sgf, sgp),down); %[n.bp n.ap]
    strue_pd(:,:,1) = downsample2(rebin_fan2par(strue(:,:,1), sgf, sgp),down);
    strue_pd(:,:,2) = downsample2(rebin_fan2par(strue(:,:,2), sgf, sgp),down); 
elseif geo_sel == 0 %parallel-beam
    ymipd(:,:,1) = downsample2(ymi(:,:,1),down);  %downsampling
    ymipd(:,:,2) = downsample2(ymi(:,:,2),down); %[n.bp n.ap]
    strue_pd(:,:,1) = downsample2(strue(:,:,1),down);
    strue_pd(:,:,2) = downsample2(strue(:,:,2),down); 
end
n.bp = n.bp/down;n.ap = n.ap/down; %downsampling
%n.bp = 1024;n.ap = 800; (original)
%save de_ct_low512_setup.mat;
