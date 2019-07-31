%| function out = strvcat(varargin)
%|
%| (for freemat)

function out = strvcat(varargin)

% for compatability with matlab, remove empty arguments
tmp = {};
for ii=1:length(varargin)
	if ~isempty(varargin{ii})
		tmp = {tmp{:}, varargin{ii}};
	end
end

out = char(tmp{:});
