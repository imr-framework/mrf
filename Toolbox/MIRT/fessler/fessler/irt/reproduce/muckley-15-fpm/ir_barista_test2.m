% ir_barista_test2
% by Matt Muckley
% for the real brain data results
% todo: needs the data to be read from web or such

nx = 144;
ny = 256;
nc = 8;

load braindat;
load kmask;
% load brainmaskxinf;
load brainxinf;
% load braindb4xinf;

A = (1/sqrt(nx*ny))*Gdft('mask', true(nx,ny), 'fftshift', 1, ...
    'ifftshift', 1);

dat = A*reshape(coilim, [nx*ny nc]); clear coilim;

% mask = imdilate(mask, true(10,10));
% % mask = repmat(mask, [1 1 nc]);
% smap(~mask) = 0.1;
mask = true(nx,ny);

kmask = kmask((256-nx)/2+1:256-(256-nx)/2,:);
dat = col(dat(col(kmask),:));

n = kmask; b = 1:nx*ny;
n = b(n); clear b;

Q = (1/sqrt(nx*ny))*Gdft('ifftshift', 1, 'fftshift', 1, ...
    'samp', kmask);
S = cell(nc,1);
F = cell(nc,1);
for i=1:nc
    S{i} = Gdiag(smap(:,:,i));
    F{i} = Q;
end
S = block_fatrix(S, 'type', 'col');
F = block_fatrix(F, 'type', 'diag');
A = F*S;
% A = Apsm('knownfn', 'time', 'v', 1, 'n', n(:), 'smap', smap, 'immask', ...
%     true(nx,ny), 'nk', nx*ny);

x = A'*dat; x = x./ col(sum(abs(smap).^2,3));

% mask = abs(xinf) > 0.1*max(col(abs(xinf)));
% mask = true(nx,ny);
% x = x(mask);

xinit = x;

[~, info] = ir_mri_barista_wavelet(dat, mask, smap, kmask, 'niter', 500, ...
    'x0', xinit, 'beta', 2^8, 'xinf', xinf, 'wtype', 'haar', 'kmask', kmask, ...
    'timeonly', 0, 'restart', 0, 'obarista', 0);
nrbaristadist = info.dist;

[~, info] = ir_mri_barista_wavelet(dat, mask, smap, kmask, 'niter', 500, ...
    'x0', xinit, 'beta', 2^8, 'xinf', xinf, 'wtype', 'haar', 'kmask', kmask, ...
    'timeonly', 0, 'restart', 1, 'obarista', 0);
baristadist = info.dist;

[~, info] = ir_mri_barista_wavelet(dat, mask, smap, kmask, 'niter', 500, ...
    'x0', xinit, 'beta', 2^8, 'xinf', xinf, 'wtype', 'haar', 'kmask', kmask, ...
    'timeonly', 0, 'restart', 1, 'obarista', 0, 'fista', 1);
fistadist = info.dist;

[~, info] = ir_mri_barista_wavelet(dat, mask, smap, kmask, 'niter', 500, ...
    'x0', xinit, 'beta', 2^8, 'xinf', xinf, 'wtype', 'haar', 'kmask', kmask, ...
    'timeonly', 0, 'restart', 0, 'obarista', 0, 'fista', 1);
nrfistadist = info.dist;

xinf = col(xinf);
plot(20*log10(nrfistadist/norm(xinf)), 'Color', [0.3 0.8 0.3], 'LineWidth', 2);
hold on;
plot(20*log10(nrbaristadist/norm(xinf)), 'Color', [0.5 0.5 1], 'LineWidth', 2);
hold on;
plot(20*log10(fistadist/norm(xinf)), 'Color', [0 0.4 0], 'LineWidth', 2); 
hold on;
plot(20*log10(baristadist/norm(xinf)), 'b', 'LineWidth', 2);
axis([0 450 -250 0]);
legend('FISTA', 'RFISTA', 'NRBARISTA', 'BARISTA');
xlabel('Iteration Number');
ylabel('$\xi(k)$ (dB)', 'interpreter', 'latex');
