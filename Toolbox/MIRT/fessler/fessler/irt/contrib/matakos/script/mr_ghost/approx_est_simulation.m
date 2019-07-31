load ../../data/mri_ghost_correction/data_img_256x256x160
dataprefix = path_find_dir('data');
datadir = [dataprefix '/ge/2013-05-14/'];
load([datadir 'mat/setup_mc_axial']);

N = [192 256];
nc = 4;
f.fov = [18 24];

ig = image_geom('nx',N(1),'dx',f.fov(1)/N(1),'ny',N(2),'dy', f.fov(2)/N(2),...
	'offset_x', 0.5, 'offset_y', 0.5);

smap = mri_sensemap_sim('nx', ig.nx, 'ny', ig.ny, 'dx', 10*ig.dx, ...
	'dy', 10*ig.dy, 'rcoil', 120, 'ncoil', nc, 'chat',0);

% nc = 1;
% f.smap = ones(N);
f.smap = smap;

f.mask = true(N);
f.nufft = {N, [1 1], N, N/2, 'linear'};
% f.dkx = 1.2/f.fov(1);
f.dkx = 1.21726354/f.fov(1);
f.dky = 0.217263549/f.fov(2);
f.phi = [0.1229846371; 2*pi*f.fov(1)/N(1)*f.dkx; 2*pi*f.fov(2)/N(2)*f.dky];

[ks,om,wi] = mri_trajectory('epi',{1},N,f.fov);

ksx = reshape(ks(:,1),N(1),[]);
ksxp = ksx(:,1:2:end); % data going right
ksxn = ksx(:,2:2:end); % data goint left

ksy = reshape(ks(:,2),N(1),[]);
ksyp = ksy(:,1:2:end); % data going right
ksyn = ksy(:,2:2:end); % data goint left

ksp = [ksxp(:) ksyp(:)];
ksn = [ksxn(:) ksyn(:)];

xx = (1:N(1))'-N(1)/2-1;
yy = (1:N(2))'-N(2)/2-1;
[X,Y] = ndgrid(xx,yy);

T = [ones(sum(f.mask(:)),1) X(f.mask) Y(f.mask)];
pmap = T*f.phi;


G = Gmri(ks,f.mask,'fov',f.fov,'basis',{'dirac'},'nufft',f.nufft);
Gp = Gmri(ksp,f.mask,'fov',f.fov,'basis',{'dirac'},'nufft',f.nufft);
Gn = Gmri(ksn,f.mask,'fov',f.fov,'basis',{'dirac'},'nufft',f.nufft);

Gpc = cell(nc,1);
Gnc = cell(nc,1);

for ic=1:nc
	Gpc{ic} = Gp * Gdiag(f.smap(:,:,ic));
	Gnc{ic} = Gn * Gdiag(f.smap(:,:,ic));
end
Gp = block_fatrix(Gpc,'type','col');
Gn = block_fatrix(Gnc,'type','col');

Gm = block_fatrix({Gp,Gn*Gdiag(exp(-1i*pmap(f.mask)),'mask',f.mask)},'type','col');

xt = img(:,:,80);
xt = xt(33:end-32,:);

mask = xt > 10;
mask = imdilate(mask, ones(5));

maskp = xor(mask,circshift(mask,[0 128])) & mask;

yt = Gm*xt(:);
y = reshape(yt, [], 2);
yp = reshape(y(:,1), N(1), N(2)/2, nc);
yn = reshape(y(:,2), N(1), N(2)/2, nc);
y = zeros([N nc]);
y(:,1:2:end,:) = yp;
y(:,2:2:end,:) = yn;
y = reshape(y,[],1,nc);

xres = zeros(192,256,8);
idx = 1;

for mm=[1 2 4 8 12 18 24]
	[~, ~, ~, phi, pmap] = ghost_correct(y,ks,f.smap,'etype','approx2d',...
		'mask',maskp,'fov',f.fov,'chat',1,'plot',0,'basis','dirac','fft_osmp',mm); 

	if mm==1
		phi = 0*phi;
		pmap = 0*pmap;
	end
	
	printm('%dx oversampling\n',mm);
	nrms(phi(1),f.phi(1));
	nrms(phi(2),f.phi(2));
	nrms(phi(3),f.phi(3));
	
	Gm = block_fatrix({Gp,Gn*Gdiag(exp(-1i*pmap(f.mask)),'mask',f.mask)},'type','col');

	ticker reset, ticker print, ticker(4);
	xi = qpwls_pcg1(0*xt(:), Gm, 1, yt, 0, 'niter', 40);
	xres(:,:,idx) = abs(embed(xi,f.mask));
	idx = idx+1;
end

[~, ~, ~, phi, pmap] = ghost_correct(y,ks,f.smap,'etype','vproj2d',...
		'mask',maskp,'fov',f.fov,'chat',1,'plot',0,'basis','dirac','fft_osmp',4); 
printm('QS estimate oversampling\n',mm);
nrms(phi(1),f.phi(1));
nrms(phi(2),f.phi(2));
nrms(phi(3),f.phi(3));

Gm = block_fatrix({Gp,Gn*Gdiag(exp(-1i*pmap(f.mask)),'mask',f.mask)},'type','col');

ticker reset, ticker print, ticker(4);
xi = qpwls_pcg1(0*xt(:), Gm, 1, yt, 0, 'niter', 40);
xres(:,:,idx) = abs(embed(xi,f.mask));

nrms(xres,xt);

cc = colormap('gray');

c1 = imadjust(cc,[],[],3);
c2 = imadjust(cc,[],[],.5);

figure(1), im('notick',xres,' ','cbar',[0 255]), colormap(cc);
figure(2), im('notick',xres,' ','cbar',[0 255]), colormap(c1);
figure(3), im('notick',xres,' ','cbar',[0 255]), colormap(c2);
figure(4), im('notick',log10(xres),' ','cbar',[0 2.5]), colormap(cc);



