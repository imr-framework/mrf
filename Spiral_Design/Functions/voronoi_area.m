function dcf_out = voronoi_area(k)
%dcf = voronoi_area(k)
%  Calculate the area (density) of a k-space trajectory via Voronoi
%  triangulation
%
%   k     k-space trajectory  ([-0.5 0.5]*res)
%   dcf  area assosiated with each k-space point
%
% Created Florian Wiesinger
% Modified 7/2006 Rolf Schulte
%
% See also VORONOIN.
if (nargin<1), help(mfilename); dcf = NaN; return; end;
if size(k,1)~=1, error('size(k,1)~=1'); end

fufa = 0.7; % fudge factor to reduce outside areas

tmp = version;
if str2double(tmp(1))>7
    [ku,ki] = unique(k,'stable');
else
    [ku,ki] = unique(k);
end
nupts = length(k)-length(ki);
nsimi = sum(abs(diff(sort(k)))<0.02);
if nsimi>nupts, 
    warning('voronoi_area:nsimi',...
        'Detected %g similar points -> check dcf',nsimi);
end
if nupts>0, 
    fprintf('%g non-unique points excluded\n',nupts);
end
kx  = real(ku).';
ky  = imag(ku).';
kxy = [kx,ky];
[V,C] = voronoin(kxy);      % calculate voronoi vertices
dcf  = zeros(1,length(kxy));
makx  = max(abs(kx));
maky  = max(abs(ky));
makxy = min([makx maky]);
tmp1 = 0;

for i1 = 1:length(C)%length(kxy),     % loop through all triangles
    x = V(C{i1},1);
    y = V(C{i1},2);
    % check for points outside the sampled area
    if  any(isinf(x)) || any(isinf(y)) || ...
            any(abs(x)>makx) || any(abs(y)>maky) || ...
            any(sqrt(x.^2+y.^2)>makxy), 
        % outside: radius -> approximate area
        tmp1 = tmp1+1;
        radii = sort(sqrt( (x-kx(i1)).^2 + (y-ky(i1)).^2 ));
        dcf(i1) = radii(1)*radii(2)*fufa;
        % dcf(i1) = mean(radii(1:2))^2;
        % dcf(i1) = mean(radii(1:2))^2*pi;
    else
        % inside: correct area
        dcf(i1) = polyarea(x,y);
    end
end

if any(isinf(dcf(:))),
    dcf(isinf(dcf)) = 0;
    warning('voronoi_area:inf','dcf contains inf -> setting to zero');
end
dcf_out = zeros(1,length(k));
dcf_out(1,ki) = dcf;

% 
mean_dcf_unsamp = mean(dcf_out(dcf_out>1));
if mean_dcf_unsamp>1, 
    warning('voronoi_area:unsamp','Numerically stable?');
    fprintf('mean(dcf(dcf>1))=%g\n',mean_dcf_unsamp);
    figure,plot(dcf_out);
end

