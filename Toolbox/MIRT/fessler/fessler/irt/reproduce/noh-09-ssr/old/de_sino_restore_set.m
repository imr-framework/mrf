%This is the code for DE-CT sinogram restoration
%Modified from spie02.m

clear all;
im off;
%%%%%%%%%%%%
%Load data %
%%%%%%%%%%%%
%load xtrue.mat; %128 by 104
load xtrue1024.mat; %1024 by 1024
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
%G = Gtomo2_strip(sg, ig); %small size CT
G = Gtomo2_dscmex(sg, ig); %full size CT
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

save de_ct_low1024_setup.mat;
