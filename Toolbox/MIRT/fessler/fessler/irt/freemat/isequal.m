 function out = isequal(a,b)
%function out = isequal(a,b)
% needed for freemat 4.0

if nargin < 2, help(mfilename), error(mfilename), end

out = false;

% check that class type matches
if ~isnumeric(a) || ~isnumeric(b)
	if ~streq(class(a), class(b)), return, end
end

if ischar(a)
	out = streq(a, b);
return
end

if iscell(a)
	out = isequal_cell(a,b);
return
end

if ~isnumeric(a) || ~isnumeric(b)
	fail('this class not done')
end

% if here, then numeric!

% check dimension
sa = size(a);
sb = size(b);
if length(sa) ~= length(sb), return, end

% check size
if any(sa ~= sb), return, end

% check values
out = all(a == b);

%
% isequal_cell()
%
function out = isequal_cell(a, b)
out = false;

% check dimension
sa = size(a);
sb = size(b);
if length(sa) ~= length(sb), return, end

% check size
if any(sa ~= sb), return, end

% check values
for ii=1:prod(sa)
	if ~isequal(a{ii}, b{ii}), return, end
end
out = true;
