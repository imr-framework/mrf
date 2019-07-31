dataprefix = path_find_dir('data');
datadir = [dataprefix '/ge/2013-05-14/'];
trajdir = [dataprefix '/trajectories/'];

f.fov = [24 24];
f.niter = 10;
f.updates = 50;
f.chat = 1;
f.plot = 1;
f.basis = 'dirac';

ang = pi/6;
R = [cos(ang) -sin(ang);sin(ang) cos(ang)];
cc = colormap('gray');
c1 = imadjust(cc,[],[],4);
c2 = imadjust(cc,[],[],.5);

%% Single coil axial plane
%
load([datadir 'mat/setup_sc_axial']);
N = size(xref_cg);
smap = ones(N);
% mask = imdilate(mask_obj, ones(5));
mask = true(N);
maskv = mask & xor(mask, circshift(mask, [0 N(2)/2]));

% 1) Single shot
% kss = (R*kss.').'; %#ok<SUSENS>

% 1D linear phase correction
mask_est = false(N);
mask_est(8:47,27:43) = true;

[ks_ss1d, y_ss1d, dkf_ss1d, phi_ss1d, pmap_ss1d, out_ss1d] = ghost_correct(yss, ...
	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly10');

% 2D phase map fitting correction
mask_est = false(N);
mask_est(8:47,27:43) = true;

[ks_ss2d, y_ss2d, dkf_ss2d, phi_ss2d, pmap_ss2d, out_ss2d] = ghost_correct(yss, ...
	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly11');

% Model based joint estimation 1st and 2nd order
[ks_ssj1, y_ssj1, dkf_ssj1, phi_ssj1, pmap_ssj1, out_ssj1] = ghost_correct(yss, ...
	kss, smap, 'mask', maskv, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'vproj2d', ...
	'updates', f.updates, 'fft_osmp', 1, 'itnrms', 1e-3);

% [ks_ssj2, y_ssj2, dkf_ssj2, phi_ssj2, pmap_ssj2, out_ssj2] = ghost_correct(yss, ...
% 	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
% 	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
% 	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2);

xres = zeros(64,64,3);
xres(:,:,1) = out_ss1d.xmod;
xres(:,:,2) = out_ss2d.xmod;
xres(:,:,3) = out_ssj1.xmod;

% xres(end,:,1:2) = 28;
figure(2),im('row',1,xres,' ','cbar',[0 26]),colormap(cc);
figure(3),im('row',1,xres,' ','cbar',[0 26]),colormap(c1);
figure(4),im('row',1,xres,' ','cbar',[0 26]),colormap(c2);


% 2) 2 shot

% 1D linear phase correction
% mask_est = false(N);
% mask_est(8:47,27:43) = true;

[ks_s21d, y_s21d, dkf_s21d, phi_s21d, pmap_s21d, out_s21d] = ghost_correct(ys2, ...
	ks2, smap, 'phi', phi_ss1d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly10');

% 2D phase map fitting correction

[ks_s22d, y_s22d, dkf_s22d, phi_s22d, pmap_s22d, out_s22d] = ghost_correct(ys2, ...
	ks2, smap, 'phi', phi_ss2d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly11');

% Model based joint estimation 1st and 2nd order (universal pmap)
% [ks_s2j1, y_s2j1, dkf_s2j1, phi_s2j1, pmap_s2j1, out_s2j1] = ghost_correct(ys2, ...
% 	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
% 	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
% 	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 0);

% [ks_s2j2, y_s2j2, dkf_s2j2, phi_s2j2, pmap_s2j2, out_s2j2] = ghost_correct(ys2, ...
% 	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
% 	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
% 	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 0);

% Model based joint estimation 1st and 2nd order (individual pmap)
[ks_s2j1i, y_s2j1i, dkf_s2j1i, phi_s2j1i, pmap_s2j1i, out_s2j1i] = ghost_correct(ys2, ...
	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 1, 'itnrms', 1e-4);

% [ks_s2j2i, y_s2j2i, dkf_s2j2i, phi_s2j2i, pmap_s2j2i, out_s2j2i] = ghost_correct(ys2, ...
% 	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
% 	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
% 	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 1);

xres = zeros(64,64,3);
xres(:,:,1) = out_s21d.xmod;
xres(:,:,2) = out_s22d.xmod;
xres(:,:,3) = out_s2j1i.xmod;

figure(2),im('row',1,xres,' ','cbar',[0 26]),colormap(cc);
figure(3),im('row',1,xres,' ','cbar',[0 26]),colormap(c1);
figure(4),im('row',1,xres,' ','cbar',[0 26]),colormap(c2);


% xres(end,:,1:2) = 28;

% im('notick',[xres(:,:,1);xres(:,:,2);xres(:,:,3)],' ',[0 28]);

% xind(:,:,1) = out_s2j1.xmod;
% xind(:,:,2) = out_s2j1i.xmod;
% 
% im('notick',abs(xind(:,:,2)),' ',[0 27]);


% 3) 4 shot

% 1D linear phase correction
% mask_est = false(N);
% mask_est(8:47,27:43) = true;

[ks_s41d, y_s41d, dkf_s41d, phi_s41d, pmap_s41d, out_s41d] = ghost_correct(ys4, ...
	ks4, smap, 'phi', phi_ss1d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly10');

% 2D phase map fitting correction

[ks_s42d, y_s42d, dkf_s42d, phi_s42d, pmap_s42d, out_s42d] = ghost_correct(ys4, ...
	ks4, smap, 'phi', phi_ss2d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly11');

% Model based joint estimation 1st and 2nd order (universal pmap)
% [ks_s4j1, y_s4j1, dkf_s4j1, phi_s4j1, pmap_s4j1, out_s4j1] = ghost_correct(ys4, ...
% 	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
% 	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
% 	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 0);
% 
% [ks_s4j2, y_s4j2, dkf_s4j2, phi_s4j2, pmap_s4j2, out_s4j2] = ghost_correct(ys4, ...
% 	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
% 	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
% 	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 0);

% Model based joint estimation 1st and 2nd order (individual pmap)
[ks_s4j1i, y_s4j1i, dkf_s4j1i, phi_s4j1i, pmap_s4j1i, out_s4j1i] = ghost_correct(ys4, ...
	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 1);
% 
% [ks_s4j2i, y_s4j2i, dkf_s4j2i, phi_s4j2i, pmap_s4j2i, out_s4j2i] = ghost_correct(ys4, ...
% 	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
% 	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
% 	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 1);

xres = zeros(64,64,3);
xres(:,:,1) = out_s41d.xmod;
xres(:,:,2) = out_s42d.xmod;
xres(:,:,3) = out_s4j1i.xmod;

figure(2),im('row',1,xres,' ','cbar',[0 26]),colormap(cc);
figure(3),im('row',1,xres,' ','cbar',[0 26]),colormap(c1);
figure(4),im('row',1,xres,' ','cbar',[0 26]),colormap(c2);


% xres(end,:,1:2) = 28;
% 
% im('notick',[xres(:,:,1);xres(:,:,2);xres(:,:,3)],' ',[0 28]);
% 
% xind(:,:,1) = out_s4j1.xmod;
% xind(:,:,2) = out_s4j1i.xmod;
% 
% im('notick',abs(xind(:,:,2)),' ',[0 27]);

%}



%% Single coil single oblique plane
%{
load([datadir 'mat/setup_sc_sobl']);
N = size(xref_cg);
smap = ones(N);
mask = imdilate(mask_obj, ones(1));

% 1) Single shot

% 1D linear phase correction
mask_est = false(N);
mask_est(14:50,29:44) = true;

[ks_ss1d, y_ss1d, dkf_ss1d, phi_ss1d, pmap_ss1d, out_ss1d] = ghost_correct(yss, ...
	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly10');

% 2D phase map fitting correction
mask_est = false(N);
mask_est(14:50,29:44) = true;

[ks_ss2d, y_ss2d, dkf_ss2d, phi_ss2d, pmap_ss2d, out_ss2d] = ghost_correct(yss, ...
	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly11');

% Model based joint estimation 1st and 2nd order
[ks_ssj1, y_ssj1, dkf_ssj1, phi_ssj1, pmap_ssj1, out_ssj1] = ghost_correct(yss, ...
	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1);

[ks_ssj2, y_ssj2, dkf_ssj2, phi_ssj2, pmap_ssj2, out_ssj2] = ghost_correct(yss, ...
	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2);


% 2) 2 shot

% 1D linear phase correction
% mask_est = false(N);
% mask_est(8:47,27:43) = true;

[ks_s21d, y_s21d, dkf_s21d, phi_s21d, pmap_s21d, out_s21d] = ghost_correct(ys2, ...
	ks2, smap, 'phi', phi_ss1d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly10');

% 2D phase map fitting correction

[ks_s22d, y_s22d, dkf_s22d, phi_s22d, pmap_s22d, out_s22d] = ghost_correct(ys2, ...
	ks2, smap, 'phi', phi_ss2d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly11');

% Model based joint estimation 1st and 2nd order (universal pmap)
[ks_s2j1, y_s2j1, dkf_s2j1, phi_s2j1, pmap_s2j1, out_s2j1] = ghost_correct(ys2, ...
	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 0);

[ks_s2j2, y_s2j2, dkf_s2j2, phi_s2j2, pmap_s2j2, out_s2j2] = ghost_correct(ys2, ...
	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 0);

% Model based joint estimation 1st and 2nd order (individual pmap)
[ks_s2j1i, y_s2j1i, dkf_s2j1i, phi_s2j1i, pmap_s2j1i, out_s2j1i] = ghost_correct(ys2, ...
	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 1);

[ks_s2j2i, y_s2j2i, dkf_s2j2i, phi_s2j2i, pmap_s2j2i, out_s2j2i] = ghost_correct(ys2, ...
	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 1);


% 3) 4 shot

% 1D linear phase correction
% mask_est = false(N);
% mask_est(8:47,27:43) = true;

[ks_s41d, y_s41d, dkf_s41d, phi_s41d, pmap_s41d, out_s41d] = ghost_correct(ys4, ...
	ks4, smap, 'phi', phi_ss1d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly10');

% 2D phase map fitting correction

[ks_s42d, y_s42d, dkf_s42d, phi_s42d, pmap_s42d, out_s42d] = ghost_correct(ys4, ...
	ks4, smap, 'phi', phi_ss2d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly11');

% Model based joint estimation 1st and 2nd order (universal pmap)
[ks_s4j1, y_s4j1, dkf_s4j1, phi_s4j1, pmap_s4j1, out_s4j1] = ghost_correct(ys4, ...
	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 0);

[ks_s4j2, y_s4j2, dkf_s4j2, phi_s4j2, pmap_s4j2, out_s4j2] = ghost_correct(ys4, ...
	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 0);

% Model based joint estimation 1st and 2nd order (individual pmap)
[ks_s4j1i, y_s4j1i, dkf_s4j1i, phi_s4j1i, pmap_s4j1i, out_s4j1i] = ghost_correct(ys4, ...
	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 1);

[ks_s4j2i, y_s4j2i, dkf_s4j2i, phi_s4j2i, pmap_s4j2i, out_s4j2i] = ghost_correct(ys4, ...
	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 1);
%}



%% Single coil double oblique plane
%{
load([datadir 'mat/setup_sc_dobl']);
N = size(xref_cg);
smap = ones(N);
mask = imdilate(mask_obj, ones(5));

% 1) Single shot

% 1D linear phase correction
mask_est = false(N);
mask_est(10:50,31:45) = true;

[ks_ss1d, y_ss1d, dkf_ss1d, phi_ss1d, pmap_ss1d, out_ss1d] = ghost_correct(yss, ...
	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly10');

% 2D phase map fitting correction
mask_est = false(N);
mask_est(10:50,31:45) = true;

[ks_ss2d, y_ss2d, dkf_ss2d, phi_ss2d, pmap_ss2d, out_ss2d] = ghost_correct(yss, ...
	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly11');

% Model based joint estimation 1st and 2nd order
[ks_ssj1, y_ssj1, dkf_ssj1, phi_ssj1, pmap_ssj1, out_ssj1] = ghost_correct(yss, ...
	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1);

[ks_ssj2, y_ssj2, dkf_ssj2, phi_ssj2, pmap_ssj2, out_ssj2] = ghost_correct(yss, ...
	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2);


% 2) 2 shot

% 1D linear phase correction
% mask_est = false(N);
% mask_est(8:47,27:43) = true;

[ks_s21d, y_s21d, dkf_s21d, phi_s21d, pmap_s21d, out_s21d] = ghost_correct(ys2, ...
	ks2, smap, 'phi', phi_ss1d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly10');

% 2D phase map fitting correction

[ks_s22d, y_s22d, dkf_s22d, phi_s22d, pmap_s22d, out_s22d] = ghost_correct(ys2, ...
	ks2, smap, 'phi', phi_ss2d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly11');

% Model based joint estimation 1st and 2nd order (universal pmap)
[ks_s2j1, y_s2j1, dkf_s2j1, phi_s2j1, pmap_s2j1, out_s2j1] = ghost_correct(ys2, ...
	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 0);

[ks_s2j2, y_s2j2, dkf_s2j2, phi_s2j2, pmap_s2j2, out_s2j2] = ghost_correct(ys2, ...
	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 0);

% Model based joint estimation 1st and 2nd order (individual pmap)
[ks_s2j1i, y_s2j1i, dkf_s2j1i, phi_s2j1i, pmap_s2j1i, out_s2j1i] = ghost_correct(ys2, ...
	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 1);

[ks_s2j2i, y_s2j2i, dkf_s2j2i, phi_s2j2i, pmap_s2j2i, out_s2j2i] = ghost_correct(ys2, ...
	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 1);


% 3) 4 shot

% 1D linear phase correction
% mask_est = false(N);
% mask_est(8:47,27:43) = true;

[ks_s41d, y_s41d, dkf_s41d, phi_s41d, pmap_s41d, out_s41d] = ghost_correct(ys4, ...
	ks4, smap, 'phi', phi_ss1d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly10');

% 2D phase map fitting correction

[ks_s42d, y_s42d, dkf_s42d, phi_s42d, pmap_s42d, out_s42d] = ghost_correct(ys4, ...
	ks4, smap, 'phi', phi_ss2d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly11');

% Model based joint estimation 1st and 2nd order (universal pmap)
[ks_s4j1, y_s4j1, dkf_s4j1, phi_s4j1, pmap_s4j1, out_s4j1] = ghost_correct(ys4, ...
	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 0);

[ks_s4j2, y_s4j2, dkf_s4j2, phi_s4j2, pmap_s4j2, out_s4j2] = ghost_correct(ys4, ...
	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 0);

% Model based joint estimation 1st and 2nd order (individual pmap)
[ks_s4j1i, y_s4j1i, dkf_s4j1i, phi_s4j1i, pmap_s4j1i, out_s4j1i] = ghost_correct(ys4, ...
	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 1);

[ks_s4j2i, y_s4j2i, dkf_s4j2i, phi_s4j2i, pmap_s4j2i, out_s4j2i] = ghost_correct(ys4, ...
	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 1);
%}



%% Mutli coil axial plane
%
load([datadir 'mat/setup_mc_axial']);
N = size(xref_sos);
% smap = ones(N);
mask = imdilate(mask_obj, ones(1));

% 1) Single shot

% 1D linear phase correction
mask_est = false(N);
mask_est(11:48,22:40) = true;

[ks_ss1d, y_ss1d, dkf_ss1d, phi_ss1d, pmap_ss1d, out_ss1d] = ghost_correct(yss, ...
	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly10','indcoil',0);

% 2D phase map fitting correction
mask_est = false(N);
mask_est(11:48,22:40) = true;

[ks_ss2d, y_ss2d, dkf_ss2d, phi_ss2d, pmap_ss2d, out_ss2d] = ghost_correct(yss, ...
	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly11','indcoil',0);

% Model based joint estimation 1st and 2nd order
[ks_ssj1, y_ssj1, dkf_ssj1, phi_ssj1, pmap_ssj1, out_ssj1] = ghost_correct(yss, ...
	kss, smap, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'vproj2d', ...
	'updates', f.updates, 'fft_osmp', 16, 'itnrms', 1e-3);

% [ks_ssj2, y_ssj2, dkf_ssj2, phi_ssj2, pmap_ssj2, out_ssj2] = ghost_correct(yss, ...
% 	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
% 	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
% 	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indcoil', 0);

xres = zeros(64,64,3);
xres(:,:,1) = out_ss1d.xmod;
xres(:,:,2) = out_ss2d.xmod;
xres(:,:,3) = out_ssj1.xmod;

figure(2),im('row',1,xres,' ','cbar',[0 1.5]),colormap(cc);
figure(3),im('row',1,xres,' ','cbar',[0 1.5]),colormap(c1);
figure(4),im('row',1,xres,' ','cbar',[0 1.5]),colormap(c2);


% xres(end,:,1:2) = 1.7;
% 
% im('notick',[xres(:,:,1);xres(:,:,2);xres(:,:,3)],' ',[0 1.7]);



% 2) 2 shot

% 1D linear phase correction
% mask_est = false(N);
% mask_est(8:47,27:43) = true;

[ks_s21d, y_s21d, dkf_s21d, phi_s21d, pmap_s21d, out_s21d] = ghost_correct(ys2, ...
	ks2, smap, 'phi', phi_ss1d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly10');

% 2D phase map fitting correction

[ks_s22d, y_s22d, dkf_s22d, phi_s22d, pmap_s22d, out_s22d] = ghost_correct(ys2, ...
	ks2, smap, 'phi', phi_ss2d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly11');

% Model based joint estimation 1st and 2nd order (universal pmap)
% [ks_s2j1, y_s2j1, dkf_s2j1, phi_s2j1, pmap_s2j1, out_s2j1] = ghost_correct(ys2, ...
% 	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
% 	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
% 	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 0);

% [ks_s2j2, y_s2j2, dkf_s2j2, phi_s2j2, pmap_s2j2, out_s2j2] = ghost_correct(ys2, ...
% 	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
% 	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
% 	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 0);

% Model based joint estimation 1st and 2nd order (individual pmap)
[ks_s2j1i, y_s2j1i, dkf_s2j1i, phi_s2j1i, pmap_s2j1i, out_s2j1i] = ghost_correct(ys2, ...
	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 1, 'itnrms', 1e-4);

% [ks_s2j2i, y_s2j2i, dkf_s2j2i, phi_s2j2i, pmap_s2j2i, out_s2j2i] = ghost_correct(ys2, ...
% 	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
% 	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
% 	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 1);

xres = zeros(64,64,3);
xres(:,:,1) = out_s21d.xmod;
xres(:,:,2) = out_s22d.xmod;
xres(:,:,3) = out_s2j1i.xmod;

figure(2),im('row',1,xres,' ','cbar',[0 1.5]),colormap(cc);
figure(3),im('row',1,xres,' ','cbar',[0 1.5]),colormap(c1);
figure(4),im('row',1,xres,' ','cbar',[0 1.5]),colormap(c2);

% xres(end,:,1:2) = 1.7;
% 
% im('notick',[xres(:,:,1);xres(:,:,2);xres(:,:,3)],' ',[0 1.7]);
% 
% xind(:,:,1) = out_s2j1.xmod;
% xind(:,:,2) = out_s2j1i.xmod;
% 
% im('notick',abs(xind(:,:,2)),' ',[0 1.6]);


% 3) 4 shot

% 1D linear phase correction
% mask_est = false(N);
% mask_est(8:47,27:43) = true;

[ks_s41d, y_s41d, dkf_s41d, phi_s41d, pmap_s41d, out_s41d] = ghost_correct(ys4, ...
	ks4, smap, 'phi', phi_ss1d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly10');

% 2D phase map fitting correction

[ks_s42d, y_s42d, dkf_s42d, phi_s42d, pmap_s42d, out_s42d] = ghost_correct(ys4, ...
	ks4, smap, 'phi', phi_ss2d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly11');

% Model based joint estimation 1st and 2nd order (universal pmap)
% [ks_s4j1, y_s4j1, dkf_s4j1, phi_s4j1, pmap_s4j1, out_s4j1] = ghost_correct(ys4, ...
% 	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
% 	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
% 	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 0);

% [ks_s4j2, y_s4j2, dkf_s4j2, phi_s4j2, pmap_s4j2, out_s4j2] = ghost_correct(ys4, ...
% 	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
% 	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
% 	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 0);

% Model based joint estimation 1st and 2nd order (individual pmap)
[ks_s4j1i, y_s4j1i, dkf_s4j1i, phi_s4j1i, pmap_s4j1i, out_s4j1i] = ghost_correct(ys4, ...
	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 1, 'itnrms', 1e-4);

% [ks_s4j2i, y_s4j2i, dkf_s4j2i, phi_s4j2i, pmap_s4j2i, out_s4j2i] = ghost_correct(ys4, ...
% 	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
% 	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
% 	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 1);

xres = zeros(64,64,3);
xres(:,:,1) = out_s41d.xmod;
xres(:,:,2) = out_s42d.xmod;
xres(:,:,3) = out_s4j1i.xmod;

figure(2),im('row',1,xres,' ','cbar',[0 1.5]),colormap(cc);
figure(3),im('row',1,xres,' ','cbar',[0 1.5]),colormap(c1);
figure(4),im('row',1,xres,' ','cbar',[0 1.5]),colormap(c2);


% xres(end,:,1:2) = 1.6;
% 
% im('notick',[xres(:,:,1);xres(:,:,2);xres(:,:,3)],' ',[0 1.6]);
% 
% xind(:,:,1) = out_s4j1.xmod;
% xind(:,:,2) = out_s4j1i.xmod;
% 
% im('notick',abs(xind(:,:,2)),' ',[0 1.6]);

%}



%% Mutli coil single oblique plane
%{
load([datadir 'mat/setup_mc_sobl']);
N = size(xref_sos);
% smap = ones(N);
mask = imdilate(mask_obj, ones(1));

% 1) Single shot

% 1D linear phase correction
mask_est = false(N);
mask_est(11:48,22:40) = true;

[ks_ss1d, y_ss1d, dkf_ss1d, phi_ss1d, pmap_ss1d, out_ss1d] = ghost_correct(yss, ...
	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly10','indcoil',0);

% 2D phase map fitting correction
mask_est = false(N);
mask_est(11:48,22:40) = true;

[ks_ss2d, y_ss2d, dkf_ss2d, phi_ss2d, pmap_ss2d, out_ss2d] = ghost_correct(yss, ...
	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly11','indcoil',0);

% Model based joint estimation 1st and 2nd order
[ks_ssj1, y_ssj1, dkf_ssj1, phi_ssj1, pmap_ssj1, out_ssj1] = ghost_correct(yss, ...
	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d(:,1), 'morder', 1, 'indcoil', 0);

[ks_ssj2, y_ssj2, dkf_ssj2, phi_ssj2, pmap_ssj2, out_ssj2] = ghost_correct(yss, ...
	kss, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indcoil', 0);


% 2) 2 shot

% 1D linear phase correction
% mask_est = false(N);
% mask_est(8:47,27:43) = true;

[ks_s21d, y_s21d, dkf_s21d, phi_s21d, pmap_s21d, out_s21d] = ghost_correct(ys2, ...
	ks2, smap, 'phi', phi_ss1d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly10');

% 2D phase map fitting correction

[ks_s22d, y_s22d, dkf_s22d, phi_s22d, pmap_s22d, out_s22d] = ghost_correct(ys2, ...
	ks2, smap, 'phi', phi_ss2d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly11');

% Model based joint estimation 1st and 2nd order (universal pmap)
[ks_s2j1, y_s2j1, dkf_s2j1, phi_s2j1, pmap_s2j1, out_s2j1] = ghost_correct(ys2, ...
	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 0);

[ks_s2j2, y_s2j2, dkf_s2j2, phi_s2j2, pmap_s2j2, out_s2j2] = ghost_correct(ys2, ...
	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 0);

% Model based joint estimation 1st and 2nd order (individual pmap)
[ks_s2j1i, y_s2j1i, dkf_s2j1i, phi_s2j1i, pmap_s2j1i, out_s2j1i] = ghost_correct(ys2, ...
	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 1);

[ks_s2j2i, y_s2j2i, dkf_s2j2i, phi_s2j2i, pmap_s2j2i, out_s2j2i] = ghost_correct(ys2, ...
	ks2, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 1);


% 3) 4 shot

% 1D linear phase correction
% mask_est = false(N);
% mask_est(8:47,27:43) = true;

[ks_s41d, y_s41d, dkf_s41d, phi_s41d, pmap_s41d, out_s41d] = ghost_correct(ys4, ...
	ks4, smap, 'phi', phi_ss1d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly10');

% 2D phase map fitting correction

[ks_s42d, y_s42d, dkf_s42d, phi_s42d, pmap_s42d, out_s42d] = ghost_correct(ys4, ...
	ks4, smap, 'phi', phi_ss2d, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'poly11');

% Model based joint estimation 1st and 2nd order (universal pmap)
[ks_s4j1, y_s4j1, dkf_s4j1, phi_s4j1, pmap_s4j1, out_s4j1] = ghost_correct(ys4, ...
	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 0);

[ks_s4j2, y_s4j2, dkf_s4j2, phi_s4j2, pmap_s4j2, out_s4j2] = ghost_correct(ys4, ...
	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 0);

% Model based joint estimation 1st and 2nd order (individual pmap)
[ks_s4j1i, y_s4j1i, dkf_s4j1i, phi_s4j1i, pmap_s4j1i, out_s4j1i] = ghost_correct(ys4, ...
	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', phi_ss1d, 'morder', 1, 'indshot', 1);

[ks_s4j2i, y_s4j2i, dkf_s4j2i, phi_s4j2i, pmap_s4j2i, out_s4j2i] = ghost_correct(ys4, ...
	ks4, smap, 'emask', mask_est, 'mask', mask, 'fov', f.fov, 'niter', f.niter, ...
	'chat', f.chat, 'plot', f.plot, 'basis', f.basis, 'etype', 'joint', 'updates', ...
	f.updates, 'phiinit', [phi_ss1d;0;0;0], 'morder', 2, 'indshot', 1);
%}

