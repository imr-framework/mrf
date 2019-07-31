function [kspace, omega, wi] = mri_trajectory(ktype, arg_traj, N, fov, arg_wi)
%|function [kspace, omega, wi] = mri_trajectory(ktype, arg_traj, N, fov, arg_wi)
%| generate kspace trajectory samples and density compensation functions.
%|
%| in
%|	ktype		string	k-space trajectory type.  see choices below.
%|	arg_traj	cell	arguments for a specific trajectory
%|	N	[1,2-3]		target image size
%|	fov	[1,2-3]		field of view in x and y (and z)
%|	arg_wi		cell	optional arguments to pass to mri_density_comp
%| out
%|	kspace	[Nk,2-3]	kspace samples in units 1/fov
%|	omega	[Nk,2-3]	trajectory samples over [-pi,pi]
%|	wi	[Nk,1]		(optional) density compensation factors
%|
%| trajectory types:
%| 'cartesian' 'radial' 'cart:y/2' 'random'
%| 'half+8' 'epi-sin'
%| 'spiral0' 'spiral1' 'spiral3'
%| new types:
%| 'racetrack'	EPI with 'racetrack' oversampling of center of kspace. The
%|				racetrack pattern is described in techavipoo:08:gdc.
%|				arg_traj{1} specifies the oversampled lines at the center of
%|				kspace and must be an integer
%| 'rct'			EPI with oversampling of center of kspace. More straightforward
%|				oversampling than the racetrack option. arg_traj{1} specifies
%|				the oversampled lines and must be an integer.
%| 'epi'			EPI trajectory with undersampling option when the arg_traj input
%|				is used. The undersampling factor is agr_traj{1} and should be
%|				integer
%| 'epimod'		Interleaved EPI with 2 interleaves. First is bottom up and
%|				second is top bottom
%| 'epinter'		Interleaved EPI with variable interleaves. All go in a bottom
%|				up fashion. The arg_traj input can be used as a cell array
%|				specifying the number of interleaves (arg_traj{1}) and the sign
%|				(arg_traj{2} - optional). The sign value specifies the starting
%|				location of the EPI. 0 means left and 1 means right. Mostly used
%|				to accommodate real data from the Philips scanner.
%| 'cartmod1'	Interleaved cartesian with 2 interleaves. Both interleaves go
%|				in a bottom up way
%| 'cartmod2'	Interleaved cartesian with 2 interleaves. First is bottom up and
%|				second is top bottom
%|
%| Copyright 2004-4-21, Jeff Fessler, The University of Michigan
%| Edited 2008-06-02, Antonis Matakos, The University of Michigan

if nargin < 1 || (nargin == 1 && streq(ktype, 'test'))
	mri_trajectory_test
	return
end
if nargin < 4, help(mfilename), error 'args', end
if ~isvar('arg_wi'), arg_wi = {}; end

if length(N) == 1, N = [N N]; end
if length(fov) == 1, fov = fov * ones(size(N)); end


% default density compensation, which works for ordinary Cartesian (only)
wi = [];
if isempty('arg_wi')
	wi = 1 / prod(fov);
end

%
% trajectory choices
%

% ideal cartesian
switch ktype
case 'cartesian'
	% bug fix - now it works properly with odd dimensions
	Noff = mod(N,2);
	o1 = ((0:(N(1)-1))/N(1) - 0.5*(N(1)-Noff(1))/N(1))*2*pi;
	o2 = ((0:(N(2)-1))/N(2) - 0.5*(N(2)-Noff(2))/N(2))*2*pi;
% 	o1 = ([0:(N(1)-1)]/N(1) - 0.5)*2*pi;
% 	o2 = ([0:(N(2)-1)]/N(2) - 0.5)*2*pi;
	if length(N) == 2
		[o1 o2] = ndgrid(o1, o2);
		omega = [o1(:) o2(:)];
	elseif length(N) == 3
		o3 = ((0:(N(3)-1))/N(3) - 0.5)*2*pi;
		[o1 o2 o3] = ndgrid(o1, o2, o3);
		omega = [o1(:) o2(:) o3(:)];
	else
		error 'only 2d and 3d done'
	end
	wi = 1 / prod(fov);

case 'cartmod1'
	% bug fix - now it works properly with odd dimensions
	Noff = mod(N,2);
	o1 = ((0:(N(1)-1))/N(1) - 0.5*(N(1)-Noff(1))/N(1))*2*pi;
	o2 = ((0:(N(2)-1))/N(2) - 0.5*(N(2)-Noff(2))/N(2))*2*pi;
% 	o1 = ([0:(N(1)-1)]/N(1) - 0.5)*2*pi;
% 	o2 = ([0:(N(2)-1)]/N(2) - 0.5)*2*pi;
	if length(N) == 2
		[o1 o2] = ndgrid(o1, o2);
		o2 = [o2(:,1:2:end) o2(:,2:2:end)];
		omega = [o1(:) o2(:)];
	elseif length(N) == 3
		o3 = ((0:(N(3)-1))/N(3) - 0.5)*2*pi;
		[o1 o2 o3] = ndgrid(o1, o2, o3);
		omega = [o1(:) o2(:) o3(:)];
	else
		error 'only 2d and 3d done'
	end
	wi = 1 / prod(fov);

case 'cartmod2'
	% bug fix - now it works properly with odd dimensions
	Noff = mod(N,2);
	o1 = ((0:(N(1)-1))/N(1) - 0.5*(N(1)-Noff(1))/N(1))*2*pi;
	o2 = ((0:(N(2)-1))/N(2) - 0.5*(N(2)-Noff(2))/N(2))*2*pi;
% 	o1 = ([0:(N(1)-1)]/N(1) - 0.5)*2*pi;
% 	o2 = ([0:(N(2)-1)]/N(2) - 0.5)*2*pi;
	if length(N) == 2
		[o1 o2] = ndgrid(o1, o2);
		o2 = [o2(:,1:2:end) fliplr(o2(:,2:2:end))];
		omega = [o1(:) o2(:)];
	elseif length(N) == 3
		o3 = ((0:(N(3)-1))/N(3) - 0.5)*2*pi;
		[o1 o2 o3] = ndgrid(o1, o2, o3);
		omega = [o1(:) o2(:) o3(:)];
	else
		error 'only 2d and 3d done'
	end
	wi = 1 / prod(fov);

case 'radial'
	[omega wi] = mri_trajectory_radial(N(1:2), fov(1:2), arg_traj{:});
	if length(N) == 3 % 'stack of radials'
		omega = mri_trajectory_stack(omega, N(3));
		wi = repmat(wi, N(3), 1);
		wi = wi(:) / fov(3);
	end
	if ~isempty(arg_wi)
		wi = []; % discard the default "analytical" DCFs
	end

% half cartesian + 8 rows
case 'half+8'
	o1 = ((0:(N(1)-1))/N(1) - 0.5)*2*pi;
	o2 = (-N(2)/2:8)/N(2) * 2*pi;
	[oo1 oo2] = ndgrid(o1, o2);
	omega = [oo1(:), oo2(:)];

% echo-planar with sinusoid:
case 'epi-sin'
	if isempty(arg_traj)
		oversample = 1;
	elseif iscell(arg_traj) && length(arg_traj) == 1
		oversample = arg_traj{1};
	else
		error 'bad trajectory argument'
	end
	Npt = oversample*prod(N);
	t = (0:(Npt-1))'/Npt;
	omega = [pi*sin(2*pi*t*N(2)/2) t*2*pi-pi];

% bad spiral:
case 'spiral0'
	Nspiral = 400;
	omega = linspace(0, 10*2*pi, Nspiral)';
	omega = pi*[cos(omega) sin(omega)] .* omega(:,[1 1])/max(omega);
	if isempty(arg_wi)
		wi = abs(omega(:,1) + 1i * omega(:,2)); % simple |r| weighting
	end

% crude spiral:
case 'spiral1'
	Nspiral = round(prod(N) * pi/4);
	omega = linspace(0, N(1)*2*pi, Nspiral)';
	omega = pi*[cos(omega) sin(omega)] .* omega(:,[1 1])/max(omega);
	if isempty(arg_wi)
		wi = abs(omega(:,1) + 1i * omega(:,2)); % simple |r| weighting
	end

% 3T spiral:
case 'spiral3'
	if fov(1) ~= fov(2) || N(1) ~= N(2), error 'only square done', end
	[kspace omega] = mri_kspace_spiral('N', max(N(1:2)), ...
				'fov', max(fov(1:2)));

	if length(N) == 3, % stack of these spirals
		omega = mri_trajectory_stack(omega, N(3));
	end

	if isempty(arg_wi)
		wi = abs(omega(:,1) + 1i * omega(:,2)); % simple |r| weighting
%		wi = wi / (2 * pi * prod(N) * prod(fov));
		if length(N) == 3
			wi = wi / fov(3); % cartesian-style weighting in z
		end
	end
% Rosette
case 'rosette'
	if length(N) == 2
		if isempty(arg_traj)
			nshot = 1;
			Nt = 2*N(1)^2;
			dt = 4e-6;
		elseif iscell(arg_traj) && length(arg_traj) == 1
			nshot = arg_traj{1};
			Nt = 2*N(1)^2;
			dt = 4e-6;
		elseif iscell(arg_traj) && length(arg_traj) == 2
			nshot = arg_traj{1};
			Nt = arg_traj{2};
			dt = 4e-6;
		elseif iscell(arg_traj) && length(arg_traj) == 3
			nshot = arg_traj{1};
			Nt = arg_traj{2};
			dt = arg_traj{3};
		else
			error 'bad trajectory argument'
		end
	
		[ks, omega, gr, sl, ts, idx, wi] = make_rosette(fov(1), N(1), dt, ...
			'Nt', Nt, 'shot', nshot);
		omega = omega(idx,:);
% 		wi = [];
	elseif length(N) == 3
		[omega wi] = mri_trajectory_rosette3(N, fov);
	end
% random
case 'random'
	rand('twister', 0)
	omega = (rand(N(1)*N(2)*2, 2)-0.5)*2*pi;

% 2D FT, undersampled in "y" (phase encode) direction
case 'cart:y/2'
	o1 = ((0:(N(1)/1-1))/(N(1)/1) - 0.5)*2*pi;
	o2 = ((0:(N(2)/2-1))/(N(2)/2) - 0.5)*2*pi;
	[oo1 oo2] = ndgrid(o1, o2);
	omega = [oo1(:), oo2(:)];

case 'racetrack'
	
	% account for center of kspace oversampling
	if isempty(arg_traj)
		L = 2;
	elseif iscell(arg_traj) && length(arg_traj) == 1
		L = arg_traj{1};
	else
		error 'bad trajectory argument'
	end
	
	[omega, wi] = mri_trajectory_racetrack(N(1:2),fov(1:2),L);
	
	if length(N) == 3, % stack of these racetracks
		omega = mri_trajectory_stack(omega, N(3));
	end

case 'rct'
    
    % account for center of kspace oversampling
    if isempty(arg_traj)
		L = 2;
	elseif iscell(arg_traj) && length(arg_traj) == 1
		L = arg_traj{1};
	else
		error 'bad trajectory argument'
    end
    
    [omega, wi] = mri_trajectory_rct(N(1:2),fov(1:2),L);
    
	if length(N) == 3, % stack of these racetracks
        omega = mri_trajectory_stack(omega, N(3));
	end

case 'epi'
	
	% account for center of kspace oversampling
    if isempty(arg_traj)
		L = 1;
	elseif iscell(arg_traj) && length(arg_traj) == 1
		L = arg_traj{1};
	else
		error 'bad trajectory argument'
    end
    
    [omega, wi] = mri_trajectory_epi(N(1:2),fov(1:2),L);
    
	if length(N) == 3, % stack of these racetracks
        omega = mri_trajectory_stack(omega, N(3));
	end

case 'epimod' % new addition - testing version - seems to be OK now!
    
    % account for center of kspace oversampling
%     if isempty(arg_traj)
% 		L = 2;
% 	elseif iscell(arg_traj) && length(arg_traj) == 1
% 		L = arg_traj{1};
% 	else
% 		error 'bad trajectory argument'
%     end
    
    [omega, wi] = mri_trajectory_epimod(N(1:2),fov(1:2));
    
	if length(N) == 3, % stack of these racetracks
        omega = mri_trajectory_stack(omega, N(3));
	end

case 'epinter' % new addition - testing version - seems to be OK now!
    
    % account for center of kspace oversampling
    if isempty(arg_traj)
		L = 2;
		sgn = 1;
	elseif iscell(arg_traj) && length(arg_traj) == 1
		L = arg_traj{1};
		sgn = 1;
	elseif iscell(arg_traj) && length(arg_traj) == 2
		L = arg_traj{1};
		sgn = arg_traj{2};
	else
		error 'bad trajectory argument'
    end
    
    [omega, wi] = mri_trajectory_epinter(N(1:2),fov(1:2),L,sgn);
    
	if length(N) == 3, % stack of these racetracks
        omega = mri_trajectory_stack(omega, N(3));
	end

otherwise
	error('unknown trajectory "%s"', ktype)
end

% convert to physical units
kspace = zeros(size(omega));
for id=1:length(N)
	dx = fov(id) / N(id);
	kspace(:,id) = omega(:,id) / (2*pi) / dx;
end

if ~isempty(arg_wi) && nargout > 2 && isempty(wi)
	wi = mri_density_comp(kspace, arg_wi{:});
end

end


%
% mri_trajectory_stack()
% make 3D "stack of ..." trajectory from a 2D trajectory
%
function omega = mri_trajectory_stack(omega2, N3)
o3 = ((0:(N3-1))/N3 - 0.5)*2*pi;
o3 = repmat(o3, nrow(omega2), 1); % [N12,N3]
omega = repmat(omega2, N3, 1); % [N12*N3,2]
omega = [omega o3(:)]; % [N12*N3,3]

end


%
% mri_trajectory_radial()
% todo: generalize to 3D using barger:02:trc
%
function [omega wi] = mri_trajectory_radial(N, fov, varargin)
arg.na_nr = 2*pi;	% default ensures proper sampling at edge of k-space
arg.na = [];		% angular spokes (default: na_nr * nr)
arg.nr = max(N)/2;	% radial samples per spoke
arg.ir = [];		% default: 0:nr
arg.omax = pi;		% maximum omega
arg = vararg_pair(arg, varargin);
if isempty(arg.ir), arg.ir = (0:arg.nr); end
if isempty(arg.na), arg.na = 4*ceil(arg.na_nr * arg.nr/4); end % mult of 4
om = arg.ir/arg.nr * pi;
ang = (0:arg.na-1)/arg.na * 2*pi;
[om ang] = ndgrid(om, ang); % [nr+1, na]
omega = [col(om.*cos(ang)) col(om.*sin(ang))];

% density compensation factors based on "analytical" voronoi
if any(fov ~= fov(1)), fail('only square FOV implemented for radial'), end
du = 1/fov(1); % assume this radial sample spacing
wi = pi * du^2 / arg.na * 2 * arg.ir(:); % see lauzon:96:eop, joseph:98:sei
wi(arg.ir == 0) = pi * (du/2)^2 / arg.na; % area of center disk
wi = repmat(wi, [1 arg.na]);
wi = wi(:);

end


%
% NEW
%
% mri_trajectory_racetrack()
% note: Testing version, use at your own risk - fixed weights
%
function [omega, wi] = mri_trajectory_racetrack(N, fov, L)
Noff = mod(N,2);
o1 = ((0:(N(1)-1))/N(1) - 0.5*(N(1)-Noff(1))/N(1))*2*pi;
o2 = ((0:(N(2)-1))/N(2) - 0.5*(N(2)-Noff(2))/N(2))*2*pi;
% o1 = ([0:(N(1)-1)]/N(1) - 0.5)*2*pi;
% o2 = ([0:(N(2)-1)]/N(2) - 0.5)*2*pi;

[o1 o2] = ndgrid(o1, o2);
    
% add oversampling columns for o1
o1 = [o1 repmat(o1(:,1),[1 L])];

% create initial racetrack
o1(:,1:2:end) = flipud(o1(:,1:2:end)); % comment out for testing

mid = ceil((N(2)+1)/2) - floor((L+1)/2) + 1;

% add oversampling columns for o2
o2 = [o2(:,1:mid) repmat(o2(:,mid),[1 L]) o2(:,mid+1:end)];


for i=0:L-1
    % first replicate appropriate columns
    o2(:,mid+2*i:mid+2*i+1) = repmat(o2(:,mid+L+i),[1 2]);
    % now swap columns to correct the trajectory
    o2(:,mid+2*i-1:mid+2*i) = fliplr(o2(:,mid+2*i-1:mid+2*i));
end
% final swap
o2(:,mid+2*L-1:mid+2*L) = fliplr(o2(:,mid+2*L-1:mid+2*L));

% create initial weights
wi = 1 / prod(fov) * ones(size(o1));
% reduce duplicate weights by 1/2
wi(:,mid:mid+2*L-1) = wi(:,mid:mid+2*L-1)/2;
% swap first and last duplicate columns to match trajectory
wi(:,mid-1:mid) = fliplr(wi(:,mid-1:mid));
wi(:,mid+2*L-1:mid+2*L) = fliplr(wi(:,mid+2*L-1:mid+2*L));

omega = [o1(:) o2(:)];
wi = wi(:);


end


%
% NEW
%
% mri_trajectory_rct() - artificial testing trajectory
% note: Testing version, use at your own risk - fixed weights
%
function [omega, wi] = mri_trajectory_rct(N, fov, L)
Noff = mod(N,2);
o1 = ((0:(N(1)-1))/N(1) - 0.5*(N(1)-Noff(1))/N(1))*2*pi;
o2 = ((0:(N(2)-1))/N(2) - 0.5*(N(2)-Noff(2))/N(2))*2*pi;
% o1 = ([0:(N(1)-1)]/N(1) - 0.5)*2*pi;
% o2 = ([0:(N(2)-1)]/N(2) - 0.5)*2*pi;

[o1 o2] = ndgrid(o1, o2);

% create initial racetrack
o1(:,1:2:end) = flipud(o1(:,1:2:end)); % comment out for testing

mid = ceil((N(2)+1)/2) - floor((L+1)/2) + 1;

% add oversampling columns for o2
o1 = [o1(:,1:mid-1) o1(:,mid:mid+L-1) o1(:,mid:mid+L-1) o1(:,mid+L:end)];
o2 = [o2(:,1:mid-1) o2(:,mid:mid+L-1) o2(:,mid:mid+L-1) o2(:,mid+L:end)];

% create initial weights
wi = 1 / prod(fov) * ones(size(o1));
% reduce duplicate weights by 1/2
wi(:,mid:mid+2*L-1) = wi(:,mid:mid+2*L-1)/2;

omega = [o1(:) o2(:)];
wi = wi(:);


end



%
% NEW
%
% mri_trajectory_epi() - epi trajectory w/ undersampling
% note: Testing version, use at your own risk - fixed weights
%
function [omega, wi] = mri_trajectory_epi(N, fov, L)
Noff = mod(N,2);
o1 = ((0:(N(1)-1))/N(1) - 0.5*(N(1)-Noff(1))/N(1))*2*pi;
o2 = ((0:(N(2)-1))/N(2) - 0.5*(N(2)-Noff(2))/N(2))*2*pi;
% o1 = ([0:(N(1)-1)]/N(1) - 0.5)*2*pi;
% o2 = ([0:(N(2)-1)]/N(2) - 0.5)*2*pi;

[o1 o2] = ndgrid(o1, o2);

% undersample if necessary
o1 = o1(:,1:L:end);
o2 = o2(:,1:L:end);

% create epi
o1(:,1:2:end) = flipud(o1(:,1:2:end)); % comment out for testing

% create initial weights
wi = 1 / prod(fov);

omega = [o1(:) o2(:)];
wi = wi(:);


end


%
% NEW
%
% mri_trajectory_rct() - artificial testing trajectory
% note: Testing version, use at your own risk - fixed weights
%
function [omega, wi] = mri_trajectory_epimod(N, fov)
Noff = mod(N,2);
o1 = ((0:(N(1)-1))/N(1) - 0.5*(N(1)-Noff(1))/N(1))*2*pi;
o2 = ((0:(N(2)-1))/N(2) - 0.5*(N(2)-Noff(2))/N(2))*2*pi;
% o1 = ([0:(N(1)-1)]/N(1) - 0.5)*2*pi;
% o2 = ([0:(N(2)-1)]/N(2) - 0.5)*2*pi;

[o1 o2] = ndgrid(o1, o2);

% create initial racetrack
o1(:,1:2:end) = flipud(o1(:,1:2:end)); % comment out for testing
o2 = [o2(:,1:2:end) o2(:,2:2:end)];

% create initial weights
wi = 1 / prod(fov);

omega = [o1(:) o2(:)];
wi = wi(:);

end


function [omega, wi] = mri_trajectory_epinter(N, fov, L, sgn)
Noff = mod(N,2);
o1 = ((0:(N(1)-1))/N(1) - 0.5*(N(1)-Noff(1))/N(1))*2*pi;
o2 = ((0:(N(2)-1))/N(2) - 0.5*(N(2)-Noff(2))/N(2))*2*pi;
% o1 = ([0:(N(1)-1)]/N(1) - 0.5)*2*pi;
% o2 = ([0:(N(2)-1)]/N(2) - 0.5)*2*pi;

[o1 o2] = ndgrid(o1, o2);

% o2 = fliplr(o2);
w1 = [];
w2 = [];

for i=1:L
	tmp = o2(:,i:L:end);
	w2 = [tmp w2];
	tmp = o1(:,i:L:end);
	if sgn == 1
		tmp(:,2:2:end) = flipud(tmp(:,2:2:end));
	else
		tmp(:,1:2:end) = flipud(tmp(:,1:2:end));
	end
	w1 = [w1 tmp];
end

% create initial weights
wi = 1 / prod(fov);

% w2 = fliplr(w2);
omega = [w1(:) w2(:)];
wi = wi(:);

end



%
% mri_trajectory_rosette3()
% 3d rosette, with default parameters from bucholz:08:miw
%
function [omega wi] = mri_trajectory_rosette3(N, fov, varargin)
arg.f1 = 211;
arg.f2 = 117.13;
arg.f3 = 73.65;
arg.nshot = 32; % todo: shots vs arms
arg.omax = pi; % maximum omega
arg.nt = 12300; % # time samples (65.536 ms for 4 usec dt)
arg.dt = 4e-6; % # time sample spacing (4 usec)
arg.ti = []; % # time samples
arg = vararg_pair(arg, varargin);
if isempty(arg.ti)
	arg.ti = (0:arg.nt-1)' * arg.dt;
end
tmp = 2 * pi * arg.ti;
p1 = f1 * tmp;
p2 = f2 * tmp;
p3 = f3 * tmp;
kx = arg.omax * sin(p1) .* cos(p2) .* cos(p3);
ky = arg.omax * sin(p1) .* sin(p2) .* cos(p3); 
kz = arg.omax * sin(a1) .* sin(a3);
omega = [kx ky kz];
for is=1:(arg.nshot-1) % n-shot, rotate kx,ky by 2 pi / N
	ang = is * 2 * pi / arg.nshot;
	c = cos(ang);
	s = sin(ang);
	ox = c * kx + s * ky;
	oy = -s * kx + c * ky;
	omega = [omega; [ox oy kz]];
end

wi = omax^3 * abs( sin(p1)^2 .* cos(p1) .* cos(p3) ); % from bucholz:08:miw
end

%
% test routine
%
function mri_trajectory_test

nx = 64;
ny = 64;
ig = image_geom_mri('nx', nx, 'ny', nx, 'fov', 24);
N = ig.dim;
%N = [32 30];
%fov = 250;	% 250 mm FOV
st = mri_objects('rect2', [0 0 ig.fov/2 ig.fov/2 1]);
%keyboard

arg = {1,1*ny^2};
% arg = {};

% ktype = 'cartesian';
% ktype = 'spiral3';
% ktype = 'spiral0';
ktype = 'rosette';
% ktype = 'random';
% ktype = 'epi-sin'; arg = {2};
% ktype = 'radial'; arg = {'na_nr', pi/2};
arg_wi = {'voronoi'};
% arg_wi = {};
nufft = {N, [6 6], 2*N, N/2, 'table', 2^10, 'minmax:kb'};

[kspace omega wi] = mri_trajectory(ktype, arg, N, ig.fov, arg_wi);

printm 'setup Gnufft object'
G = Gnufft(ig.mask, {omega, N, [6 6], 2*N, N/2, 'table', 2^10, 'minmax:kb'});
% G = Gnufft(ig.mask, {omega, N, [1 1], N, N/2, 'linear'});

Gm = Gmri(kspace, ig.mask, 'fov', ig.fov, 'basis', {'rect'}, 'nufft', nufft);

R = Robject(ig.mask, 'beta', 2^1, 'potential', 'quad', 'order', 1);

printm 'setup data'
yi = st.kspace(kspace(:,1), kspace(:,2));

printm 'conj. phase reconstruction'
xcp = G' * (wi .* yi);
xcp = ig.embed(xcp);

xrec = qpwls_pcg1(xcp(ig.mask), Gm, 1, yi(:), R.C, 'niter', 40);
xrec = ig.embed(xrec(:,end));

xt = st.image(ig.xg, ig.yg);

figure,
im pl 2 2
im subplot 1, plot(omega(:,1), omega(:,2), '.-')
title(sprintf('"%s" with %d k-space samples', ktype, size(omega,1)))
axis(pi*[-1 1 -1 1]), axis_pipi, axis square

im(2, ig.x, ig.y, xt, 'f true'), cbar
im(3, ig.x, ig.y, abs(xrec), 'Iterative Recon.'), cbar

im subplot 4
ix = 1:ig.nx; iy = ig.ny/2+1;
plot(	ig.x, xt(ix,iy), '-', ...
	ig.x, real(xrec(ix,iy)), 'g.-', ...
	ig.x, imag(xrec(ix,iy)), 'y.-')
axis tight


figure,
im pl 2 2
im subplot 1, plot(omega(:,1), omega(:,2), '.-')
title(sprintf('"%s" with %d k-space samples', ktype, size(omega,1)))
axis(pi*[-1 1 -1 1]), axis_pipi, axis square

im(2, ig.x, ig.y, xt, 'f true'), cbar
im(3, ig.x, ig.y, abs(xcp), 'Conj. Phase Recon.'), cbar

im subplot 4
ix = 1:ig.nx; iy = ig.ny/2+1;
plot(	ig.x, xt(ix,iy), '-', ...
	ig.x, real(xcp(ix,iy)), 'g.-', ...
	ig.x, imag(xcp(ix,iy)), 'y.-')
axis tight


nrms(xt(:),abs(xcp(:)))
nrms(xt(:),abs(xrec(:)))

end
