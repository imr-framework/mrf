 function out = subsref(a, s)
%function out = subsref(a, s)

%if nargin < 2, help(mfilename), error(mfilename), end

s1 = s(1).subs;

switch s(1).type
case '.'
	out = a.(s1);

case '()'
	out = a(s1);

case '{}'
	out = a{s1{:}};

otherwise
	fail('type "%s" not done', s(1).type)
end

% recurse if needed
if length(s) > 2
	out = subsref(out, s(2:end));
end
