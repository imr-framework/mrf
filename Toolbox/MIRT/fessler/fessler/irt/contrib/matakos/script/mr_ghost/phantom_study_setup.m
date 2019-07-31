dataprefix = path_find_dir('data');
datadir = [dataprefix '/ge/2013-05-14/'];
trajdir = [dataprefix '/trajectories/'];

f.nx = 64;
f.ny = 64;
f.nc = 8;

N = [f.nx f.ny];

f.fov = [24 24];
f.nufft = {N, [1 1], N, N/2, 'linear'};
f.niter = 20; % for all coils
f.plot = 1;
f.chat = 1;
f.mask = true(N);

f.ep = 0;				% default parameter
f.ishosp = 0;			% default parameter
f.rhbline = 1;		% default parameter
f.ne = 1;			% # of echoes

%% Multi-coil reference images and mask

% 1) Axial plane

% cartesian acquisition
load([trajdir 'epi_s64_std/input_traj_epi_s64_std']);

f.fn = [datadir 'pfiles/P_epi_s64_std_axial_8ch.7'];
f.fnb = [datadir 'pfiles/P_epi_s64_std_axial_bd.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 1;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1:f.nc,1:f.ns,f.ne); % ind. recon
y = mean(y,4);

yb = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1,1:f.ns,f.ne); % ind. recon
yb = mean(yb,4);


Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xref_fft = ifftshift(ifft2(ifftshift(y)));
xref_conj = embed(Gd'*reshape(y,[],f.nc),f.mask)/prod(N);

xref_cg = zeros(f.nx, f.ny, f.nc);
for ii=1:f.nc
	xtmp = qpwls_pcg1(xref_fft(:,:,ii), Gd, 1, vec(y(:,:,ii)), 0,'niter', f.niter);
	xref_cg(:,:,ii) = embed(xtmp,f.mask);
end

xref_sos = sqrt(sum(abs(xref_cg).^2,3));

xbd_fft = ifftshift(ifft2(ifftshift(yb)));
xbd_conj = embed(Gd'*yb(:),f.mask)/prod(N);
xbd_cg = qpwls_pcg1(xbd_fft(f.mask), Gd, 1, yb(:), 0,'niter', f.niter);
xbd_cg = embed(xbd_cg,f.mask);

mask_obj = abs(xref_sos) > 2;

[smap,sinit] = ge_create_smaps(f.fn,N,'bodycoil',f.fnb,'l2b',-1);
[smap_nb,sinit_nb] = ge_create_smaps(f.fn,N,'l2b',-1);

ycart = y;
kcart = f.ks;

% Single shot acquisition
load([trajdir 'epi_ss_std/input_traj_epi_ss_std']);

f.fn = [datadir 'pfiles/P_epi_ss_std_axial_8ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 10;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1:f.nc,1:f.ns,f.ne); % ind. recon
y = mean(y,4);
y = y(f.idx,:);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xss_conj = embed(Gd'*reshape(y,[],f.nc),f.mask)/prod(N);
xss_cg = zeros(f.nx, f.ny, f.nc);
for ii=1:f.nc
	xtmp = qpwls_pcg1(xss_conj(:,:,ii), Gd, 1, vec(y(:,ii)), 0,'niter', f.niter);
	xss_cg(:,:,ii) = embed(xtmp,f.mask);
end

yss = y;
kss = f.ks;

% 2 shot acquisition
load([trajdir 'epi_s2_std/input_traj_epi_s2_std']);

f.fn = [datadir 'pfiles/P_epi_s2_std_axial_8ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 10;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1:f.nc,1:f.ns,f.ne); % ind. recon
y = mean(y,4);
y = y(f.idx,:,:);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xs2_conj = embed(Gd'*reshape(y,[],f.nc),f.mask)/prod(N);
xs2_cg = zeros(f.nx, f.ny, f.nc);
for ii=1:f.nc
	xtmp = qpwls_pcg1(xs2_conj(:,:,ii), Gd, 1, vec(y(:,:,ii)), 0,'niter', f.niter);
	xs2_cg(:,:,ii) = embed(xtmp,f.mask);
end

ys2 = y;
ks2 = f.ks;


% 4 shot acquisition
load([trajdir 'epi_s4_std/input_traj_epi_s4_std']);

f.fn = [datadir 'pfiles/P_epi_s4_std_axial_8ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 10;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1:f.nc,1:f.ns,f.ne); % ind. recon
y = mean(y,4);
y = y(f.idx,:,:);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xs4_conj = embed(Gd'*reshape(y,[],f.nc),f.mask)/prod(N);
xs4_cg = zeros(f.nx, f.ny, f.nc);
for ii=1:f.nc
	xtmp = qpwls_pcg1(xs4_conj(:,:,ii), Gd, 1, vec(y(:,:,ii)), 0,'niter', f.niter);
	xs4_cg(:,:,ii) = embed(xtmp,f.mask);
end

ys4 = y;
ks4 = f.ks;

save setup_mc_axial mask_obj smap sinit smap_nb sinit_nb xref_sos xref_conj ...
	xref_cg ycart kcart xss_conj xss_cg yss kss xs2_conj xs2_cg ys2 ks2 ...
	xs4_conj xs4_cg ys4 ks4 xbd_conj xbd_cg



% 2) Single oblique plane

% cartesian acquisition
load([trajdir 'epi_s64_std/input_traj_epi_s64_std']);

f.fn = [datadir 'pfiles/P_epi_s64_std_sobl_8ch.7'];
f.fnb = [datadir 'pfiles/P_epi_s64_std_sobl_bd.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 1;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1:f.nc,1:f.ns,f.ne); % ind. recon
y = mean(y,4);

yb = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1,1:f.ns,f.ne); % ind. recon
yb = mean(yb,4);


Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xref_fft = ifftshift(ifft2(ifftshift(y)));
xref_conj = embed(Gd'*reshape(y,[],f.nc),f.mask)/prod(N);

xref_cg = zeros(f.nx, f.ny, f.nc);
for ii=1:f.nc
	xtmp = qpwls_pcg1(xref_fft(:,:,ii), Gd, 1, vec(y(:,:,ii)), 0,'niter', f.niter);
	xref_cg(:,:,ii) = embed(xtmp,f.mask);
end

xref_sos = sqrt(sum(abs(xref_cg).^2,3));

xbd_fft = ifftshift(ifft2(ifftshift(yb)));
xbd_conj = embed(Gd'*yb(:),f.mask)/prod(N);
xbd_cg = qpwls_pcg1(xbd_fft(f.mask), Gd, 1, yb(:), 0,'niter', f.niter);
xbd_cg = embed(xbd_cg,f.mask);

mask_obj = abs(xref_sos) > 2;

[smap,sinit] = ge_create_smaps(f.fn,N,'bodycoil',f.fnb,'l2b',-1);
[smap_nb,sinit_nb] = ge_create_smaps(f.fn,N,'l2b',-1);

ycart = y;
kcart = f.ks;

% Single shot acquisition
load([trajdir 'epi_ss_std/input_traj_epi_ss_std']);

f.fn = [datadir 'pfiles/P_epi_ss_std_sobl_8ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 10;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1:f.nc,1:f.ns,f.ne); % ind. recon
y = mean(y,4);
y = y(f.idx,:);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xss_conj = embed(Gd'*reshape(y,[],f.nc),f.mask)/prod(N);
xss_cg = zeros(f.nx, f.ny, f.nc);
for ii=1:f.nc
	xtmp = qpwls_pcg1(xss_conj(:,:,ii), Gd, 1, vec(y(:,ii)), 0,'niter', f.niter);
	xss_cg(:,:,ii) = embed(xtmp,f.mask);
end

yss = y;
kss = f.ks;

% 2 shot acquisition
load([trajdir 'epi_s2_std/input_traj_epi_s2_std']);

f.fn = [datadir 'pfiles/P_epi_s2_std_sobl_8ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 10;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1:f.nc,1:f.ns,f.ne); % ind. recon
y = mean(y,4);
y = y(f.idx,:,:);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xs2_conj = embed(Gd'*reshape(y,[],f.nc),f.mask)/prod(N);
xs2_cg = zeros(f.nx, f.ny, f.nc);
for ii=1:f.nc
	xtmp = qpwls_pcg1(xs2_conj(:,:,ii), Gd, 1, vec(y(:,:,ii)), 0,'niter', f.niter);
	xs2_cg(:,:,ii) = embed(xtmp,f.mask);
end

ys2 = y;
ks2 = f.ks;


% 4 shot acquisition
load([trajdir 'epi_s4_std/input_traj_epi_s4_std']);

f.fn = [datadir 'pfiles/P_epi_s4_std_sobl_8ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 10;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1:f.nc,1:f.ns,f.ne); % ind. recon
y = mean(y,4);
y = y(f.idx,:,:);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xs4_conj = embed(Gd'*reshape(y,[],f.nc),f.mask)/prod(N);
xs4_cg = zeros(f.nx, f.ny, f.nc);
for ii=1:f.nc
	xtmp = qpwls_pcg1(xs4_conj(:,:,ii), Gd, 1, vec(y(:,:,ii)), 0,'niter', f.niter);
	xs4_cg(:,:,ii) = embed(xtmp,f.mask);
end

ys4 = y;
ks4 = f.ks;

save setup_mc_sobl mask_obj smap sinit smap_nb sinit_nb xref_sos xref_conj ...
	xref_cg ycart kcart xss_conj xss_cg yss kss xs2_conj xs2_cg ys2 ks2 ...
	xs4_conj xs4_cg ys4 ks4 xbd_conj xbd_cg


%% Single coil reference images and mask

% 1) Axial plane

% cartesian acquisition
load([trajdir 'epi_s64_std/input_traj_epi_s64_std']);

f.fn = [datadir 'pfiles/P_epi_s64_std_axial_1ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 1;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1,1:f.ns,f.ne); % ind. recon
y = mean(y,4);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xref_fft = ifftshift(ifft2(ifftshift(y)));
xref_conj = embed(Gd'*y(:),f.mask)/prod(N);
xref_cg = qpwls_pcg1(xref_fft(f.mask), Gd, 1, y(:), 0,'niter', f.niter);
xref_cg = embed(xref_cg,f.mask);
ycart = y;
kcart = f.ks;

mask_obj = abs(xref_cg) > 2;

xax = xref_cg;

% Single shot acquisition
load([trajdir 'epi_ss_std/input_traj_epi_ss_std']);

f.fn = [datadir 'pfiles/P_epi_ss_std_axial_1ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 10;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1,1:f.ns,f.ne); % ind. recon
y = mean(y,4);
y = y(f.idx,:);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xss_conj = embed(Gd'*y(:),f.mask)/prod(N);
xss_cg = qpwls_pcg1(xss_conj(f.mask), Gd, 1, y(:), 0,'niter', f.niter);
xss_cg = embed(xss_cg,f.mask);
yss = y;
kss = f.ks;


% 2 shot acquisition
load([trajdir 'epi_s2_std/input_traj_epi_s2_std']);

f.fn = [datadir 'pfiles/P_epi_s2_std_axial_1ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 10;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1,1:f.ns,f.ne); % ind. recon
y = mean(y,4);
y = y(f.idx,:);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xs2_conj = embed(Gd'*y(:),f.mask)/prod(N);
xs2_cg = qpwls_pcg1(xs2_conj(f.mask), Gd, 1, y(:), 0,'niter', f.niter);
xs2_cg = embed(xs2_cg,f.mask);
ys2 = y;
ks2 = f.ks;


% 4 shot acquisition
load([trajdir 'epi_s4_std/input_traj_epi_s4_std']);

f.fn = [datadir 'pfiles/P_epi_s4_std_axial_1ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 10;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1,1:f.ns,f.ne); % ind. recon
y = mean(y,4);
y = y(f.idx,:);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xs4_conj = embed(Gd'*y(:),f.mask)/prod(N);
xs4_cg = qpwls_pcg1(xs4_conj(f.mask), Gd, 1, y(:), 0,'niter', f.niter);
xs4_cg = embed(xs4_cg,f.mask);
ys4 = y;
ks4 = f.ks;

save setup_sc_axial mask_obj xref_conj xref_cg ycart kcart xss_conj xss_cg yss ...
	kss xs2_conj xs2_cg ys2 ks2 xs4_conj xs4_cg ys4 ks4



% 2) Single oblique plane

% cartesian acquisition
load([trajdir 'epi_s64_std/input_traj_epi_s64_std']);

f.fn = [datadir 'pfiles/P_epi_s64_std_sobl_1ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 1;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1,1:f.ns,f.ne); % ind. recon
y = mean(y,4);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xref_fft = ifftshift(ifft2(ifftshift(y)));
xref_conj = embed(Gd'*y(:),f.mask)/prod(N);
xref_cg = qpwls_pcg1(xref_fft(f.mask), Gd, 1, y(:), 0,'niter', f.niter);
xref_cg = embed(xref_cg,f.mask);
ycart = y;
kcart = f.ks;

mask_obj = abs(xref_cg) > 2;
xso = xref_cg;


% Single shot acquisition
load([trajdir 'epi_ss_std/input_traj_epi_ss_std']);

f.fn = [datadir 'pfiles/P_epi_ss_std_sobl_1ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 10;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1,1:f.ns,f.ne); % ind. recon
y = mean(y,4);
y = y(f.idx,:);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xss_conj = embed(Gd'*y(:),f.mask)/prod(N);
xss_cg = qpwls_pcg1(xss_conj(f.mask), Gd, 1, y(:), 0,'niter', f.niter);
xss_cg = embed(xss_cg,f.mask);
yss = y;
kss = f.ks;


% 2 shot acquisition
load([trajdir 'epi_s2_std/input_traj_epi_s2_std']);

f.fn = [datadir 'pfiles/P_epi_s2_std_sobl_1ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 10;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1,1:f.ns,f.ne); % ind. recon
y = mean(y,4);
y = y(f.idx,:);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xs2_conj = embed(Gd'*y(:),f.mask)/prod(N);
xs2_cg = qpwls_pcg1(xs2_conj(f.mask), Gd, 1, y(:), 0,'niter', f.niter);
xs2_cg = embed(xs2_cg,f.mask);
ys2 = y;
ks2 = f.ks;


% 4 shot acquisition
load([trajdir 'epi_s4_std/input_traj_epi_s4_std']);

f.fn = [datadir 'pfiles/P_epi_s4_std_sobl_1ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 10;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1,1:f.ns,f.ne); % ind. recon
y = mean(y,4);
y = y(f.idx,:);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xs4_conj = embed(Gd'*y(:),f.mask)/prod(N);
xs4_cg = qpwls_pcg1(xs4_conj(f.mask), Gd, 1, y(:), 0,'niter', f.niter);
xs4_cg = embed(xs4_cg,f.mask);
ys4 = y;
ks4 = f.ks;

save setup_sc_sobl mask_obj xref_conj xref_cg ycart kcart xss_conj xss_cg yss ...
	kss xs2_conj xs2_cg ys2 ks2 xs4_conj xs4_cg ys4 ks4



% 3) Double oblique plane

% cartesian acquisition
load([trajdir 'epi_s64_std/input_traj_epi_s64_std']);

f.fn = [datadir 'pfiles/P_epi_s64_std_dobl_1ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 1;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1,1:f.ns,f.ne); % ind. recon
y = mean(y,4);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xref_fft = ifftshift(ifft2(ifftshift(y)));
xref_conj = embed(Gd'*y(:),f.mask)/prod(N);
xref_cg = qpwls_pcg1(xref_fft(f.mask), Gd, 1, y(:), 0,'niter', f.niter);
xref_cg = embed(xref_cg,f.mask);
ycart = y;
kcart = f.ks;

mask_obj = abs(xref_cg) > 2;
xdo = xref_cg;

% Single shot acquisition
load([trajdir 'epi_ss_std/input_traj_epi_ss_std']);

f.fn = [datadir 'pfiles/P_epi_ss_std_dobl_1ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 10;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1,1:f.ns,f.ne); % ind. recon
y = mean(y,4);
y = y(f.idx,:);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xss_conj = embed(Gd'*y(:),f.mask)/prod(N);
xss_cg = qpwls_pcg1(xss_conj(f.mask), Gd, 1, y(:), 0,'niter', f.niter);
xss_cg = embed(xss_cg,f.mask);
yss = y;
kss = f.ks;


% 2 shot acquisition
load([trajdir 'epi_s2_std/input_traj_epi_s2_std']);

f.fn = [datadir 'pfiles/P_epi_s2_std_dobl_1ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 10;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1,1:f.ns,f.ne); % ind. recon
y = mean(y,4);
y = y(f.idx,:);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xs2_conj = embed(Gd'*y(:),f.mask)/prod(N);
xs2_cg = qpwls_pcg1(xs2_conj(f.mask), Gd, 1, y(:), 0,'niter', f.niter);
xs2_cg = embed(xs2_cg,f.mask);
ys2 = y;
ks2 = f.ks;


% 4 shot acquisition
load([trajdir 'epi_s4_std/input_traj_epi_s4_std']);

f.fn = [datadir 'pfiles/P_epi_s4_std_dobl_1ch.7'];
f.ks = in_traj.trajectory{1};

f.fr = in_traj.frsize;
f.ns = 10;			% # of slices
f.idx = in_traj.idx - min(in_traj.idx) + 1;

y = rawload_jfn(f.fn,f.fr,f.ep,f.ishosp,f.rhbline,1,1:f.ns,f.ne); % ind. recon
y = mean(y,4);
y = y(f.idx,:);

Gd = Gmri(f.ks, f.mask, 'fov', f.fov, 'basis', {'dirac'},'nufft', f.nufft);

xs4_conj = embed(Gd'*y(:),f.mask)/prod(N);
xs4_cg = qpwls_pcg1(xs4_conj(f.mask), Gd, 1, y(:), 0,'niter', f.niter);
xs4_cg = embed(xs4_cg,f.mask);
ys4 = y;
ks4 = f.ks;

save setup_sc_dobl mask_obj xref_conj xref_cg ycart kcart xss_conj xss_cg yss ...
	kss xs2_conj xs2_cg ys2 ks2 xs4_conj xs4_cg ys4 ks4

