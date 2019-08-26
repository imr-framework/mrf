function sphere_par = NIST_ROI_ana(map_data, map_type, map_size, fov)
%Input map mat file and return ROI analysis of a ISMRM/NIST phantom, all 
%14 spheres are detected. Their locations, radius, mean, and std are
%returned in correct order. 
% INPUT
%   map_data  parameter map, unit should be in second 
%   map_type  type of map (T1 or T2)
%   map_size  size of map, for example, 128 means the map size is 128x128
%        fov  field of view                                  [mm]
% 
% OUTPUT
%   sphere_par  circle parameters. It has four fields which contain sphere
%               locations, radius, mean, and std in correct order. 
%
% 7/2019 Enlin Qian
% # Copyright of the Board of Trustees of Columbia University in the City of New York

%% Check inputs
if ~exist('map_data','var'), error('No map data is found'); end
if ~exist('map_type','var'), error('No map type is found'); end
if map_size<1, warning('dcs:map_size','map_size(=%g)<1',map_size); input(''); end
if fov<1, warning('dcs:fov','fov(=%g)<1',fov); input(''); end

%% Load standard parameters for ISMRM/NIST phantom
NIST_std_par = gen_NIST_std_par(map_type, map_size, fov);

%% Threshold Sphere 1 based on intensity
map_binary = map_data>=NIST_std_par.intensity_range(1) & map_data<=NIST_std_par.intensity_range(2);

%% Initial guess of circles to find Sphere 1
[guess_center,guess_radii] = imfindcircles(map_binary,NIST_std_par.sphere_radius,...
    'Sensitivity',0.95); 

%% Visualize initial guess of Sphere 1
figure;imagesc(map_data); colormap hot; axis equal tight;
viscircles(guess_center, guess_radii);

%% If only 1 sphere is detected, if yes, use the circle, if no, pick the one closest to true value
num_guess = size(guess_center, 1);
if num_guess == 1
    Sphere1_center = guess_center;
    Sphere1_radii = guess_radii;
else
    Mean_all = zeros(num_guess, 1);
    for n = 1: num_guess
        [~, PointsIndex] = Circle_Analysis(guess_center(n, :), guess_radii(n), map_size, map_size);
        Mean_all(n) = mean(map_data(PointsIndex));
    end
    [~, Index] = min(abs(Mean_all - NIST_std_par.true_value(1)));
    Sphere1_center = guess_center(Index, :);
    Sphere1_radii = guess_radii(Index);
end

%% Find all template sphere locations based on Sphere 1 location
NIST_sphere_loc = gen_NIST_sphere_loc(map_type, map_size, fov, Sphere1_center);

%% Visualize the template (without rotation)
% figure;imagesc(map_data); colormap hot; axis equal tight;
% viscircles(NIST_sphere_loc.template_loc, repmat(Sphere1_radii, NIST_sphere_loc.num_spheres, 1));

%% Calculate transformation matrix for rotation
rot_center = NIST_sphere_loc.temp_rot_center;
template_vec = [NIST_sphere_loc.template_Sphere1_loc(1)-NIST_sphere_loc.temp_rot_center(1),...
    NIST_sphere_loc.template_Sphere1_loc(2)-NIST_sphere_loc.temp_rot_center(2)]; % Vector for template 
map_vec = [Sphere1_center(1)-rot_center(1), Sphere1_center(2)-rot_center(2)]; % Vector for image data
rot_angle = acosd(dot(template_vec,map_vec)/(norm(template_vec)*norm(map_vec))); % Find the angle between two vectors
trans_matrix = [cosd(rot_angle), -sind(rot_angle);
                sind(rot_angle), cosd(rot_angle)];
            
%% Multiply the template with the transformation matirx
map_loc = (NIST_sphere_loc.template_loc - rot_center) * trans_matrix + rot_center; % rotation is respect to origin
map_radii = repmat(Sphere1_radii, NIST_sphere_loc.num_spheres, 1);

%% Visualize the new locations with rotation fixed
figure;imagesc(map_data); colormap hot; axis equal tight;
viscircles(map_loc, map_radii);

%% Calculate mean and SD for each sphere
for nx = 1: NIST_std_par.num_spheres
    [~, PointsIndex] = Circle_Analysis(map_loc(nx, :), map_radii(nx), map_size, map_size);
    Mean_all_sphere(nx) = mean(map_data(PointsIndex));
    std_all_sphere(nx) = std(map_data(PointsIndex));
end

%% Construct output variable circle_par
sphere_par.all_sphere_loc = map_loc;
sphere_par.all_sphere_rad = map_radii;
sphere_par.all_sphere_mean = Mean_all_sphere';
sphere_par.all_sphere_std = std_all_sphere';
