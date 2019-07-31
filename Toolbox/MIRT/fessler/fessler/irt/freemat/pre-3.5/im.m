 function h = im(varargin)
%function h = im([options,] [xx,] [yy,] zz, [scale|clim,] [title])
% show matrix zz as an image, possibly with (x,y) axes labeled by xx,yy
% options:
%	subplot		an integer for subplot
%	'notick'	no axis tick marks
%	'colorneg'	negatives as red
%	'black0'	make sure image value of 0 is black (max white)
%
% after options:
%	scale		scale image by this factor
%	clim		limits for colorbar
%	title		axis title
%
% Copyright 1997, Jeff Fessler, The University of Michigan

%
% handle state of display or not
%
persistent Display
persistent DoColorNeg
if isempty(Display)
	Display = true;
end
if isempty(DoColorNeg)
	DoColorNeg = false;
end

% default is to give help
if ~nargin & ~nargout
	help(mfilename)
	if im, disp('im enabled'), else, disp('im disabled'), end
	error(mfilename)
end

% for conditional of the form 'if im, ..., end'
if ~nargin & nargout
	h = Display;
	return
end

%
% process single string command arguments
%
if nargin == 1 & ischar(varargin{1})
	arg = varargin{1};

	if streq(arg, 'on')
		Display = true;
		disp 'enabling images'
	elseif streq(arg, 'off')
		Display = false;
		disp 'disabling images'
	elseif streq(arg, 'clf')
		if Display
			clf
		end
	elseif streq(arg, 'ison')	% query
		if Display
			h = true;
		else
			h = false;
		end
	elseif streq(arg, 'colorneg')
		DoColorNeg = true;
	elseif streq(arg, 'reset')
		DoColorNeg = false;
	else
		error 'unknown argument'
	end
return
end

scale = 1;
titlearg = {};
axisstr = 'image';
clim = [];
xx = [];
yy = [];
colorneg = false;	% put negatives in color and 0=blue ?
isxy = false;
isplot = false;
tick = true;
isblack0 = false;

%
% optional arguments
%
zz_arg_index = 1;
while notempty(varargin)
	arg = varargin{1};
	if isempty(arg)
		0; % do nothing
	elseif max(size(arg)) == 1
		if Display, subplot(varargin{1}), end
	elseif streq(arg, 'notick')
		tick = false;
	elseif streq(arg, 'colorneg')
		colorneg = true;
	elseif streq(arg, 'black0')
		isblack0 = true;
	else
		break
	end

	varargin = {varargin{2:end}};
	zz_arg_index = 1 + zz_arg_index;
end

% xx, yy
if isempty(varargin), help(mfilename), error('args'), end

if ndims(varargin{1}) <= 2 & min(size(varargin{1})) == 1
	if length(varargin) < 3
		isplot = true;
		plot_data = varargin{1};
	else
		xx = varargin{1};
		if length(varargin) < 2, help(mfilename), error('need both xx,yy'), end
		if (min(size(varargin{2})) ~= 1), error('both xx,yy need to be 1D'), end
		yy = varargin{2};
		varargin = {varargin{3:end}};
		isxy = 1;
		zz_arg_index = 1 + zz_arg_index;
	end
end

if ~Display, disp(['im disabled: ' inputname(zz_arg_index)]), return, end

% zz
if notempty(varargin)
	zz = double(varargin{1});
	if ndims(zz) == 3 & min(size(zz)) == 1
		zz = squeeze(zz);	% handle [1,n2,n3] case as [n2,n3]
	end
	varargin = {varargin{2:end}};
else
	error 'no image?'
end

% title, scale
while notempty(varargin)
	arg = varargin{1};
	if isempty(arg)
		0; % do nothing
	elseif ischar(arg)
		if strcmp(arg, 'equal')
			axisstr = arg;
		else
			titlearg = {arg};
			if ~isempty(strfind(arg, '\tex'))
				titlearg = {arg, 'interpreter', 'none'};
			end
		end
	elseif isa(arg, 'double')
		if max(size(arg)) == 1
			scale = arg;
		elseif all(size(arg) == [1 2])
			clim = arg;
		else
			error 'nonscalar scale / nonpair clim?'
		end
%		disp(['scale = ' num2str(scale)])
	else
		error 'unknown arg'
	end
	varargin = {varargin{2:end}};
end

if isplot
	plot(plot_data)
	if ~isempty(titlearg), title(titlearg{:}), end
return
end

if issparse(zz), zz = full(zz); end
if ~isreal(zz)
	zz = abs(zz);
	printf('warn %s: magnitude of complex image', mfilename)
end 

zmax = max(zz(:));
zmin = min(zz(:));

if isblack0
	if ~isempty(clim)
		warning 'black0 overrules clim'
	end
	clim = [0 zmax];
end

if (scale ~= 1)
	if scale == 0
		zmin = 0;
		scale = 1;
	elseif scale < 0
		zmin = 0;
		scale = -scale;
	end
end

if colorneg | DoColorNeg
	cmap = hot(512);	cmap = flipud(cmap(1:256,:));
	cmap = [cmap; [0 0 1]; gray(256)];
	colormap(cmap)
	zt = zz;
	zz(zt > 0) = 257+1+floor(255 * zz(zt > 0) / (abs(max(zt(:))) + eps));
	zz(zt == 0) = 257;
	zz(zt < 0) = 1+floor(-255 * zz(zt < 0) / (abs(min(zt(:))) + eps));
else
	colormap(gray)
end

if (scale ~= 1)	% fix: use Clim?
	n = nrow(colormap);
	if (zmax ~= zmin)
		zz = (n - 1) * (zz - zmin) / (zmax - zmin);
	else
		if zmin == 0
			zz(:) = 0;
			clim = [0 1];
		else
			zz(:) = n - 1;
		end
	end
	zz = 1 + round(scale * zz);
	zz = min(zz,n);
	zz = max(zz,1);
elseif zmin == zmax
	if zmin == 0
		clim = [0 1];
	else
		zz(:) = 1;
		clim = [0 1];
	end
end

if ndims(zz) < 3
	zz = zz';
	if isxy
		hh = image(xx, yy, zz, 'CDataMapping', 'scaled');
	else
		hh = image(zz, 'CDataMapping', 'scaled');
%		set(gca,'CLimMode','auto')
		% unclutter axes by only showing end limits
		n1 = size(zz,2);
		n2 = size(zz,1);
%		set(gca, 'xtick', [1 n1], 'ytick', [1 n2])

		if tick
			set(gca, 'xtick', [], 'ytick', [1 n2])
			text(1, 1.04*n2, '1', ...
				'horizontalalign', 'left', ...
				'verticalalign', 'top')
			text(n1, 1.04*n2, num2str(n1), ...
				'horizontalalign', 'right', ...
				'verticalalign', 'top')
		end
		% problem with this is it fails to register
		% space used for 'ylabel'
		% so i stick with manual xtick since that is what overlaps
		% with the colorbar
		if 0 & tick
			text(-0.01*n1, n2, '1', ...
				'verticalalign', 'bottom', ...
				'horizontalalign', 'right')
			text(-0.01*n1, 1, num2str(n2), ...
				'verticalalign', 'top', ...
				'horizontalalign', 'right')
		end
	end
else
	n1 = size(zz,1);
	n2 = size(zz,2);
	if 1
		zz(:,end+1,:) = max(zz(:));	% white border
		zz(end+1,:,:) = max(zz(:));
	end
	zz = montager(zz);
	hh = image(zz', 'CDataMapping', 'scaled');
	if tick & n1 > 1 & n2 > 1	% unclutter
		set(gca, 'xtick', [1 n1], 'ytick', [1 n2])
	end
end

if (zmax == zmin)
%	fprintf('Uniform image %g [', zmin)
%	fprintf(' %g', size(zz))
%	fprintf(' ]\n')
	texts(0.5, 0.5, sprintf('Uniform: %g', zmin), 'center', 'color', 'blue')
end

if colorneg | DoColorNeg
	set(hh, 'CDataMapping', 'direct')
else
	if ~isempty(clim),
		set(gca, 'CLim', clim)
	elseif ~ishold
		set(gca, 'CLimMode', 'auto')
	end
end

set(gca, 'TickDir', 'out')

if nargout > 0
	h = hh;
end

axis(axisstr)
% flip if y goes from -ymin to ymax
if isxy & yy(1) < 0, axis('xy'), end

if ~isempty(titlearg)
	title(titlearg{:})
%	drawnow
else
	title(sprintf('Range: [%g %g]', zmin, zmax))
end

if ~tick
	set(gca, 'xtick', [], 'ytick', [])
end
