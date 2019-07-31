 function n = nrms(x, xtrue, arg)
%function n = nrms(x, xtrue, arg)
% normalized rms error
% in
%	x	[[Nd], Nrep]
%	xtrue	[[Nd]]
%	arg			use '2' to do along 2nd dim for same size
% out
%		scalar
%		vector if x and xtrue have same size and arg = 2
% multiple x's, one xtrue
if nargin < 2, help(mfilename), error(mfilename), end

if nargin > 2
	if any(size(x) ~= size(xtrue))
		error 'need same size'
	end
	n = mean(abs(xtrue - x).^2);
return
end

xtrue = xtrue(:);
np = length(xtrue);
nrep = prod(size(x)) / np;
x = reshape(x, np, nrep);
n = abs(x - xtrue(:,ones(1,nrep)));
n = sqrt(sum(n.^2)') / norm(xtrue);
