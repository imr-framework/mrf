load ../../data/mri_ghost_correction/data_img_256x256x160

N = [192 192];
N2 = [192 216];

f.smap = ones(N);
f.fov = [24 24];
f.fovc = [24 48];

f.mask = true(N);
f.nufft = {N, [1 1], N, N/2, 'linear'};
% f.dkx = 1.2/f.fov(1);
f.dkx = 1.21726354/f.fov(1);
% f.dkxi = 0.8/f.fov(1);
f.dkxi = 0;
f.dky = 0.217263549/f.fov(2);
% f.dky = 0.051957234/f.fov(2);
% f.dky = 0;
% f.dkyi = 0.1/f.fov(2);
f.dkyi = 0;
f.phii = [0; 2*pi*f.fov(1)/N(1)*f.dkxi; 2*pi*f.fov(2)/N(2)*f.dkyi];
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

xx = (1:N2(1))'-N2(1)/2-1;
yy = (1:N2(2))'-N2(2)/2-1;
[X,Y] = ndgrid(xx,yy);

T2 = [ones(prod(N2),1) X(:) Y(:)];
pmap2 = T2*f.phi;
% f.mask = msk2;


G = Gmri(ks,f.mask,'fov',f.fov,'basis',{'dirac'},'nufft',f.nufft);
Gp = Gmri(ksp,f.mask,'fov',f.fov,'basis',{'dirac'},'nufft',f.nufft);
Gn = Gmri(ksn,f.mask,'fov',f.fov,'basis',{'dirac'},'nufft',f.nufft);
Gn1 = Gn * Gdiag(exp(-1i*pmap(f.mask)),'mask',f.mask);

I1 = Gdiag(ones(12*192,1));
I2 = Gdiag(ones(168*192,1));
Z1 = zero_fatrix([12 12]*192);
Z2 = zero_fatrix([12 168]*192);
Z3 = zero_fatrix([168 12]*192);


R1 = block_fatrix({Z1,I1,Z2,Z1,I1},'type','row');
R2 = block_fatrix({Z3,Z3,I2,Z3,Z3},'type','row');
R3 = block_fatrix({I1,Z1,Z2,I1,Z1},'type','row');
Pm = block_fatrix({R1,R2,R3},'type','col');
Gr = block_fatrix({Gp*Pm,Gn*Pm*Gdiag(exp(-1i*pmap2))},'type','col');
Gm = block_fatrix({Gp,Gn1},'type','col');

xt = img(:,:,80);
xt = xt(33:end-32,:);
xt = padarray(xt,[0 64],0,'both');
maske = xt > 10;
% maske = imdilate(maske,ones(5));
maske = xor(maske,circshift(maske,[0 192])) & maske;
maske = maske(:,85:300);
maskn = maske(:, 13:end-12);
xt2 = xt(:,85:300);
xt = reshape(Pm*xt2(:), N);

mask = xt > 10;
mask = imdilate(mask, ones(5));

maskp = padarray(xor(mask,circshift(mask,[0 96])) & mask, [0 12],'replicate','both');
% maskc = imdilate(xt2>10, ones(5));
maskc = padarray(mask, [0 12],'replicate','both');
% maskr = padarray(mask, [0 12],false,'both');
% maskr = ~maskr;
maskr = ~maskc;
maskf = false(size(xt2));
maskf(:,13:end-12) = true;

y = Gr*xt2(:);


M = Gdiag2(ones(sum(maskc(:)),1),'mask_out',maskc(:));
R = Gdiag(ones(sum(maskr(:)),1),'mask',maskr(:));
Grc = Gr*M;

ticker reset, ticker print, ticker(4);
xi = qpwls_pcg1(0*xt(:), Gm, 1, y, 0, 'niter', 40);
xhm = abs(embed(xi, maskf));
xtm = abs(embed(xi, f.mask));

xi = qpwls_pcg1(0*xt2(maskc), Grc, 1, y, 0, 'niter', 40);
xhc = abs(embed(xi, maskc));
xtc = embed(Pm*xhc(:), f.mask);

beta = sqrt(10000);
xi = qpwls_pcg1(0*xt2(:), Gr, 1, y, beta*R,'niter', 40);
xhr = abs(embed(xi, true(N2)));
xtr = embed(Pm*xhr(:), f.mask);


nrms(xt(maskn),xtm(maskn))
nrms(xt(maskn),xtc(maskn))
nrms(xt(maskn),xtr(maskn))





