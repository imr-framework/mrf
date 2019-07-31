 function xtick(arg)
%function xtick(arg)
%	set axis xticks to just end points

if ~nargin
	lim = get(gca, 'xlim');
	if lim(1) == -lim(2)
		lim = [lim(1) 0 lim(2)];
	end
	set(gca, 'xtick', lim)

elseif nargin == 1
	if ischar(arg) & streq(arg, 'off')
		set(gca, 'xtick', [])
		set(gca, 'xticklabel', [])
	else
		set(gca, 'xtick', arg)
	end
else
	error arg
end
