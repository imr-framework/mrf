function y = mldivide(a,b)
% To calculate inv(a)*b.
% Copyright Jan. 25, 2008, Dr.WEN You-Wei
% email: wenyouwei@graduate.hku.hk

if isa(a,'BlurMatrix')
    error('the first parameter must be the class of BlurMatrix');
end

if isa(b,'BlurMatrix')
	y = b;
	y.eigblurmatrix = a.eigblurmatrix .\ b.eigblurmatrix;
else
	switch a.boundarycond
		case 'cir'
			y = ifft2((a.eigblurmatrix) .\ fft2(b));
	    case 'refl'
			y = idct2((a.eigblurmatrix) .\ dct2(b));
		case 'psf'
			error('Wrong');
	end

	y = real(y);
end
