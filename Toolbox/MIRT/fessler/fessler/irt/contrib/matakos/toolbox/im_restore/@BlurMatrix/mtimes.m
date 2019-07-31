function y = mtimes(a,x)
%function y = mtimes(a, x)
% y = G * x	 
% Copyright Jan. 25, 2008, Dr.WEN You-Wei
% email: wenyouwei@graduate.hku.hk
if ~isa(a, 'BlurMatrix') 
    if numel(a) ~= 1
        error('Wrong, the first paramter must be a scalar or a class of BlurMatrix');
    else
        y = x;
        y.eigblurmatrix = a * x.eigblurmatrix;
    end
else
    if isa(x,'BlurMatrix')
        y = a;
        y.eigblurmatrix = y.eigblurmatrix .* x.eigblurmatrix;
    else
        if numel(x) == 1
            y = a;
            y.eigblurmatrix = y.eigblurmatrix  * x;
        else
        switch a.boundarycond
            case 'cir'
                y = ifft2((a.eigblurmatrix) .* fft2(x));
                y = real(y);
            case 'refl'
                for k = 1:size(a.psf,3)
                    y(:,:,k) = idct2((a.eigblurmatrix(:,:,k)) .* dct2(x(:,:,k)));
                end
            case 'imfilter'    
                [n, m, r] = size(a.eigblurmatrix);
                for j = 1:r
                    h = a.eigblurmatrix;
                    y(:,:,j) = imfilter(x(:,:,j), h(:,:,j),'symmetric');
                end
                %s = fix((size(a.eigblurmatrix)+1)/2);
                %[n, m] = size(x);
                %bx = [x(:,s(2):-1:1) x x(:,m:-1:m-s(2)+1)];
                %bx = [bx(s(1):-1:1,:); bx; bx(n:-1:n-s(1)+1,:)];
                %temp = conv2(bx,a.eigblurmatrix);
                %y = temp(2*s(1)+1:n+2*s(1),2*s(2)+1:m+2*s(2));
        end
        end
        
    end
end