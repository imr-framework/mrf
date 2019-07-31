load ../results/mat/res_quality

nb = 4;
ne = 18;

% xra = xr_ap_sep_wave(nb+1:end-nb,nb+1:end-nb);
% xre = xr_ex_sep_wave(nb+1:end-nb,nb+1:end-nb);

% xra = xr_ap_sep_tviso(nb+1:end-nb,nb+1:end-nb);
% xre = xr_ex_sep_tviso(nb+1:end-nb,nb+1:end-nb);

% xra = xr_ap_sep_tvaniso(nb+1:end-nb,nb+1:end-nb);
% xre = xr_ex_sep_tvaniso(nb+1:end-nb,nb+1:end-nb);

% xra = xr_ap_uni_wave(nb+1:end-nb,nb+1:end-nb);
% xre = xr_ex_uni_wave(nb+1:end-nb,nb+1:end-nb);

% xra = xr_ap_uni_tviso(nb+1:end-nb,nb+1:end-nb);
% xre = xr_ex_uni_tviso(nb+1:end-nb,nb+1:end-nb);

% xra = xr_ap_uni_tvaniso(nb+1:end-nb,nb+1:end-nb);
% xre = xr_ex_uni_tvaniso(nb+1:end-nb,nb+1:end-nb);

xra = xres(nb+1:end-nb,nb+1:end-nb,1);
xre = xres(nb+1:end-nb,nb+1:end-nb,2);

xt = double(imread('cameraman.tif')).';
xt = xt/max(xt(:));

xt = xt(nb+1:end-nb,nb+1:end-nb);
[nx ny] = size(xt);

mask = false(nx,ny);
mask(ne+1:end-ne,ne+1:end-ne) = true;

clim = [0 1.1];
climd = [-0.5 0.5];

nrms(xra,xt)
nrms(xre,xt)

nrms(xra(mask),xt(mask))
nrms(xre(mask),xt(mask))

nrms(xra(~mask),xt(~mask))
nrms(xre(~mask),xt(~mask))

figure(1), im(xra,'Approximate image',clim);
figure(2), im(xre,'Exact image',clim);

figure(3), im(xra-xt,'Difference of approximate image',climd,'cbar');
figure(4), im(xre-xt,'Difference of exact image',climd,'cbar');
