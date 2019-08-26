function NIST_sphere_loc = gen_NIST_sphere_loc(map_type, map_size, fov, Sphere1_center)
%Input map parameters and Sphere 1 loc to generate a map that finds all
%locations of other spheres.
% INPUT
%   map_type  type of map (T1 or T2)
%   map_size  size of map, for example, 128 means the map size is 128x128
%        fov  field of view                                  [mm]
%   Sphere1_loc  centers of sphere 1 (x, y)
%
% OUTPUT
%   NIST_sphere_loc  sphere locations for all spheres
%
% 7/2019 Enlin Qian
% # Copyright of the Board of Trustees of Columbia University in the City of New York
NIST_sphere_loc.num_spheres = 14;
NIST_sphere_loc.voxel_size = fov/map_size; % unit in mm/voxel
NIST_sphere_loc.temp_rot_center = [map_size/2, map_size/2]; % unit in pixel
switch(map_type)
    case 'T1'
        NIST_sphere_loc.dist_map = [0, 31.324, 58.893, 81.810, 95.898, 100.317,...
            95.991, 81.077, 58.860, 31.451, 35.380, 35.931, 73.124, 73.157]; % unit in mm
        NIST_sphere_loc.ang_map = [0, 70.0850, 53.3538, 35.3270, 16.9426, -0.6324,...
            -18.7705, -36.9326, -56.4287, -72.7011, 33.6653, -34.0284, 14.9729, -16.4091]; % unit in degree
            
    case 'T2'
        NIST_sphere_loc.dist_map = [0, 31.547, 58.642, 80.697, 95.713, 99.812,...
            96.050, 80.936, 58.824, 30.036, 36.204, 35.882, 73.002, 72.552]; % unit in mm
        NIST_sphere_loc.ang_map = [0, 71.4899, 54.0302, 35.7043, 17.5603, -0.1549,...
            -18.3830, -35.8010, -54.2831, -70.3638, 33.1854, -33.0284, 15.2130, -16.2861]; % unit in degree
end

%% convert mm to voxel based on map_size and fov
NIST_sphere_loc.dist_map = NIST_sphere_loc.dist_map/NIST_sphere_loc.voxel_size; % convert to unit in voxel

%% Use Sphere 1 loc to calculate template Sphere1 loc
dist_Sphere12center = sqrt((Sphere1_center(1)- NIST_sphere_loc.temp_rot_center(1))^2+...
    (Sphere1_center(2)- NIST_sphere_loc.temp_rot_center(2))^2);
NIST_sphere_loc.template_Sphere1_loc = [map_size/2, map_size/2-dist_Sphere12center]; % pixel

%% get coordinates for all spheres
NIST_sphere_loc.template_loc = zeros(NIST_sphere_loc.num_spheres, 2);
NIST_sphere_loc.template_loc(:, 1) = NIST_sphere_loc.template_Sphere1_loc(1, 1) +...
    NIST_sphere_loc.dist_map.*sind(NIST_sphere_loc.ang_map); % x axis
NIST_sphere_loc.template_loc(:, 2) = NIST_sphere_loc.template_Sphere1_loc(1, 2) +...
    NIST_sphere_loc.dist_map.*cosd(NIST_sphere_loc.ang_map); % y axis 

