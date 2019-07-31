function ob = BlurMatrix(psf, imsize, flag)
% Copyright Jan. 25, 2008, Dr.WEN You-Wei
% email: wenyouwei@graduate.hku.hk
if nargin < 3, flag = 'cir'; end
boundarycond = flag;

if nargin == 0	% required by Mathworks
	ob = class(ob, 'BlurMatrix'); return
end

if numel(imsize) ==2,imsize = [imsize 1]; end

h = psf;
if size(psf,3)==1 && imsize(3)>1
    for j = 1:imsize(3)
        h(:,:,j) = psf;
    end
end

psf = h;

ob.psf = psf;
ob.boundarycond = boundarycond;
ob.eigblurmatrix = zeros(imsize);
s = fix(([size(psf,1) size(psf,2)]+1)/2);

switch ob.boundarycond
    case 'cir'
        ob.eigblurmatrix(1:size(psf,1), 1:size(psf,2),:) = psf;
        ob.eigblurmatrix = circshift(ob.eigblurmatrix, 1-s);
        ob.eigblurmatrix = fft2(ob.eigblurmatrix);
    case 'refl'
        h1 = psf(s(1):end, s(2):end,:);
        h2 = h1(2:end,:,:); h2(s(1),:,:) = 0;
        h3 = h1 + h2;
        h4 = h3(:,2:end,:); h4(:,s(2),:) = 0;
        h5 = h4 + h3; 
        
        ob.eigblurmatrix(1:s(1),1:s(2),:) = h5;
        e1 = zeros(imsize); e1(1,1,:) = 1;
        for k = 1:imsize(3)
            ob.eigblurmatrix(:,:,k) = dct2(ob.eigblurmatrix(:,:,k))./dct2(e1(:,:,k));
        end
        %ob.ForTF         = @(x)dct2(x);
        %ob.BackTF        = @(x)idct2(x);
    case 'imfilter'
        ob.eigblurmatrix = psf;
    otherwise
        error('Wrong. Flag should be cir or refl or imfilter');
end

ob = class(ob, 'BlurMatrix');
