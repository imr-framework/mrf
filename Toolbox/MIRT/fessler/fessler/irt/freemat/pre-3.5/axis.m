 function axis(varargin)
% placeholder

persistent warned

if isempty(warned)
	warned = 0;
end

if ~warned
	warning(['function ' mfilename ' not yet implemented, ignoring!'])
	warned = 1;
end
