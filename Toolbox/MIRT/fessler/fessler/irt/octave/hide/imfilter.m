function out = imfilter(f, h, varargin)
if numel(varargin)
	arg = varargin{1};
	if streq(arg, 'circular')
		fail 'circular not implemented in fake imfilter for octave'
	end
%	if streq(arg, 'same') % todo
end
out = conv2(f, h, 'same');
