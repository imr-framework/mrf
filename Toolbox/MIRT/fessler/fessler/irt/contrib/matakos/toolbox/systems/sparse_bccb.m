function S = sparse_bccb(A,varargin)
%| function sparse_bccb(A,varargin)
%|
%| Create sparse matlab matrix for a BCCB matrix. Input can be a fatrix object
%| or a psf. In case of psf the "image" size must be specified
%|
%| in
%|	A		[np nt]		fatrix or
%|			[nbx nby]	psf
%|
%| options
%|	nx		scalar		x size of "image"
%|	ny		scalar		y size of "image"
%|	tol		scalar		tolerance for finding non-zero elements of psf

if nargin < 1
	sparse_bccb_test;
	return;
end

arg.nx = 16;
arg.ny = 16;
arg.tol = 1e-7;

arg = vararg_pair(arg, varargin);

if strcmpi(class(A),'fatrix2')
	if isfield(A.arg, 'psf')
		arg.tol = min(abs(A.arg.psf(:)))/2;
	end
else
	arg.tol = min(abs(A(:)))/2;
	
	mask = true(arg.nx,arg.ny);
	A = Gblur(mask,'psf',A,'type','fft,same');
end

nx = A.odim(1);
ny = A.odim(2);

if length(A.odim) > 2
	nz = A.odim(3);
else
	nz = 1;
end

psf_full = reshape(A(:,1),A.odim);
psf_full(abs(psf_full) < arg.tol) = 0;

% Find actual length of psf
nbx = zeros(nz,1);
nby = zeros(nz,1);

for iz = 1:nz
	[I J] = find(psf_full(:,:,iz));
	nbx(iz) = length(unique(I));
	nby(iz) = length(unique(J));
end

iiz = zeros(sum(nbx.*nby)*nx*ny,1);
jjz = zeros(sum(nbx.*nby)*nx*ny,1);
bbz = zeros(sum(nbx.*nby)*nx*ny,1);

for iz = 1:nz
	psf = psf_full(:,:,iz);
	
	% Keep only the nbx x nby non-zero entries of psf
	indx = (find(psf(:,1))-1).';
	indy = find(psf(1,:))-1;
	
	[ix jx] = ndgrid((0:nx-1).',(0:nx-1).');
	nij = mod(ix-jx,nx);
	
	inc = vec(repmat((0:ny-1)*nx,[nbx(iz)*nx 1]));
	
	mkx = false(nx);
	for kk = indx
		mkx = mkx | (nij == kk);
	end
	
	% Initialize sparse matrix components
	len = nbx(iz)*nby(iz)*nx*ny;
	lenx = nbx(iz)*nx*ny;
	
	ii = zeros(len,1);
	jj = zeros(len,1);
	bb = zeros(len,1);
	
	for kk = 0:length(indy)-1
		
		bt = reshape(psf(nij+1,indy(kk+1)+1),nx,nx);
		bb(kk*lenx+1:(kk+1)*lenx) = repmat(bt(mkx),[ny 1]);
		
		it = ix(mkx) + indy(kk+1)*nx;
		ii(kk*lenx+1:(kk+1)*lenx) = mod(repmat(it,[ny 1]) + inc,nx*ny);
		
		jt = jx(mkx);
		jj(kk*lenx+1:(kk+1)*lenx) = repmat(jt,[ny 1]) + inc;
		
	end
	
	iiz((iz-1)*len+1:iz*len) = ii + (iz-1)*nx*ny;
	jjz((iz-1)*len+1:iz*len) = jj;
	bbz((iz-1)*len+1:iz*len) = bb;
end

S = sparse(iiz+1,jjz+1,bbz);

end


function sparse_bccb_test

nx = 2048;
ny = 2048;

mask = true(nx,ny);

psf = [1 2 3 4 5;
	2 3 4 5 6;
	3 4 5 6 7;
	2 3 4 5 6;
	1 2 3 4 5];

Ab = Gblur(mask,'psf',psf,'type','fft,same');
Ac = Gwave2('mask',mask,'nlevel',1);

As1 = sparse_bccb(Ab);
As2 = sparse_bccb(psf,'nx',nx,'ny',ny);
As3 = sparse_bccb(Ac,'nx',nx,'ny',ny);

% figure(1), im(As1);
% figure(2), im(As2);
% figure(3), im(As3);

end