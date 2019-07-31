function x = cconvn(a,b,N)
%| function x = cconvn(a,b,N)
%|
%| Function to implement N-dimensional circular convolution of matrices
%| using FFT.
%| 
%| If N is scalar an N point convolution is calculated in all dimensions
%| If N is a vector it specifies how many point will be the convolution output
%| in each dimension. N must be 1xM where M = ndims(x)
%| If N is ommited the output has a size of size(a,i) + size(b,i) - 1 in
%| each dimension
%| 
%| Also matrices a and b must have the same dimensions
%|
%| input:
%|	a: the first matrix of the convolution
%|	b: The second matrix of the convolution
%|	M: The scalar or vector specifying the output size of the convolution
%|		in each dimension
%|
%| output:
%|	x: convolution output. It is a matrix of the same dimension as a and b
%|		with size specified by N
%|
%| Copyright 2008-02-17, Antonis Matakos, The University of Michigan
%| 
%| Caution: This is a testing version. Use at your own risk.

if nargin < 2
    error('cconvn:NoInputs','No input arguments specified. There should be at least two input arguments.')
end

if ndims(a) ~= ndims(b)
    error('cconvn:InvalidInput','Invalid input matrices. Matrices a and b must have the same dimensions.')
end

is_vec = 0;

if length(a) == numel(a) && length(b) == numel(b) % a and b are vectors
    if (size(a,1) == 1 && size(b,1) == 1) || (size(a,2) == 1 && size(b,2) == 1)
        siz_a = length(a);
        siz_b = length(b);
        n = 1;
        is_vec = 1;
    else
        error('cconvn:InvalidInput','Invalid input vectors. Vectors a and b must be both column vectors or row vectors')
    end
elseif (length(a) == numel(a) && length(b) ~= numel(b)) || (length(a) ~= numel(a) && length(b) == numel(b))
    error('cconvn:InvalidInput','Invalid input matrices. Matrices a and b must have the same dimensions.')
else
    siz_a = size(a);
    siz_b = size(b);
    n = ndims(a);
end

if nargin < 3 % N is not specified. Create a default
    N = siz_a + siz_b - 1;
else
    if is_vec == 1
        N = max(N);
    elseif isscalar(N)
        N = N*ones(n,1);
    elseif all(size(N) > 1)
        error('cconvn:InvalidInput','Invalid input N. Must be a vector or scalar')
    elseif length(N) ~= n
        error('cconvn:InvalidInput','Invalid input N. length(N) must be equal to ndims(a)')
    elseif ~all(N)
        error('cconvn:InvalidInput','Invalid input N. Vector M must contain no zeros')
    end
end

% The actual calculation of circular convolution using the fftn
x = ifftn(fftn(datawrapn(a,N)).*fftn(datawrapn(b,N)));

