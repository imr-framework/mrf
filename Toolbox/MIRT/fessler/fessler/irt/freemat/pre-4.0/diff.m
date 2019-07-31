function y = diff(x, id)

if nargin < 2
	if size(x,1) > 1
		id = 1; % col vector
	else
		id = 2; % row vector
	end
end

dim = size(x);
tmp = [prod(dim(1:id-1)) dim(id) prod(dim(id+1:end))];
x = reshape(x, tmp);
y = x(:,2:end,:) - x(:,1:end-1,:);
dim(id) = dim(id) - 1;
y = reshape(y, dim);
