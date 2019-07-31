 function [fw, angle, rad] = fwhm2(image, pixelsize, chat)
%function [fw, angle, rad] = fwhm2(image, pixelsize, chat)
% compute 2d fwhm of point-spread function
% uses image maximum as center and the contourc function

if nargin < 2
	pixelsize = 1;
end
if ~isvar('chat') | isempty(chat), chat = logical(0); end

% image maximum
if min(size(image)) < 11
	image = padn(image, max(size(image), 23));
end
[nx, ny] = size(image);
ii = imax(image, 2);

% find better center estimate by local centroid
cx = ii(1);
cy = ii(2);
if 1
	ix = [-5:5]';
	iy = [-5:5]';
	t = image(cx + ix, cy + iy);
	if chat
		im(ix, iy, t)
	end
	o.x = 0; o.y = 0;
	o.x = sum(t,2)'/sum(t(:)) * ix;
	o.y = sum(t,1)/sum(t(:)) * iy;
	cx = cx + o.x;
	cy = cy + o.y;
	if 0
		hold on
		plot(o.x, o.y, 'rx')
		hold off
	return
	end
end

cc = contourc(image, [1e30 0.5 * max(image(:))]);
if isempty(cc), error 'empty contour?  check minimum!', end
cc = cc(:,2:length(cc))';
cc = cc(:,[2 1]);	% swap row,col or x,y

% check center pixel found
if chat
	clf
	im(121, image)
	hold on
	plot(cc(:,1), cc(:,2), '+')
	plot(cx, cy, 'rx')
	title(sprintf('length(cc)=%d', length(cc)))
	hold off
%	axis([40 60 50 70])
	prompt
end

x = cc(:,1) - cx;
y = cc(:,2) - cy;
rsamp = sqrt(x.^2+y.^2) * pixelsize;
tsamp = 180/pi*atan2(y,x);

angle = [0:180]';
t2 = [0:180]';
r1 = interp1x(tsamp, rsamp, angle);
r2 = interp1x(tsamp, rsamp, angle-180);

% plot(tsamp, rsamp, 'o', angle, r1, '-', angle-180, r2, '-')

rad = r1 + r2;
fw = mean(rad);

if chat
	subplot(122)
	plot(angle, rad, '-o')
	xlabel 'Angle [degrees]'
	ylabel 'FWHM [pixels]'
	zoom on
end
