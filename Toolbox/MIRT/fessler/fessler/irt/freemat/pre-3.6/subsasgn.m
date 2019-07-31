 function out = subsasgn(a, s, b)
%function out = subsasgn(a, s, b)

if nargin == 1 && streq(a, 'test'), subsasgn_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

out = a;
% out = builtin('subsasgn', varargin{:}); % seems to recurse

% a.s = b;
if length(s) == 1 && streq(s(1).type, '.')
	out.(s(1).subs) = b;

% a.s1.s2 = b;
elseif length(s) == 2 && streq(s(1).type, '.') && streq(s(2).type, '.')
	s1 = s(1).subs;
	s2 = s(2).subs;
	if isfield(out, s1)
		tmp = out.(s1);
	else
		tmp = struct;
	end
	tmp.(s2) = b;
	out.(s1) = tmp;
else
	fail('not done')
end

function subsasgn_test
a.b = 1;
s(1).type = '.';
s(1).subs = 'c';
o1 = subsasgn(a, s, 2);
jf_equal(o1.c, 2)

s = struct('type', {'.', '.'}, 'subs', {'d', 'e'});
o2 = subsasgn(a, s, 3);
jf_equal(o2.d.e, 3)
