 function out = strmatch(str, strs, type)
%function out = strmatch(str, strs, type)
% return matching row indices

if nargin == 1 && streq(str, 'test'), strmatch_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end
if ~streq(type, 'exact'), error('only "exact" done'), end

if iscell(strs)
	nstr = length(strs);
	out = zeros(nstr, 1, 'logical');
	for ii=1:nstr
		out(ii) = isequal(strs{ii}, str);
	end
	out = find(out);

elseif ischar(strs)
	nstr = size(strs,1);
	out = zeros(nstr, 1, 'logical');
	for ii=1:nstr
		% ignore trailing spaces!
		tmp = strs(ii,:);
		space = ' ';
		if tmp(end) == ' ' && length(tmp) > length(str)
			space = space(1,ones(1,length(tmp) - length(str)));
			str_try = [str space];
		else
			str_try = str;
		end
		out(ii) = isequal(strs(ii,:), str_try);
	end
	out = find(out);

else
	error('bad class')
end

function strmatch_test
out = strmatch('no', {'yes', 'no', 'maybe', 'no'}, 'exact');
jf_equal(out, [2 4]')
tmp = strvcat('yes', 'no', 'maybe', 'no');
out = strmatch('no', tmp, 'exact');
jf_equal(out, [2 4]')
