function NIST_std_par = gen_NIST_std_par(map_type, map_size, fov)
%Input map parameters and return phantom parameters based on map type, map
%size and fov. 
% INPUT
%   map_type  type of map (T1 or T2)
%   map_size  size of map, for example, 128 means the map size is 128x128
%        fov  field of view                                  [mm]
%
% OUTPUT
% phan_std_par  standard phantom parameters.
%
% 7/2019 Enlin Qian
% # Copyright of the Board of Trustees of Columbia University in the City of New York
NIST_std_par.num_spheres = 14;
NIST_std_par.voxel_size = fov/map_size; % unit in mm/voxel
NIST_std_par.sphere_radius = 8; % unit in mm
NIST_std_par.temp_rot_center = [map_size/2, map_size/2];

switch(map_type)
    case 'T1'
        NIST_std_par.true_value = [1.989, 1.454, 0.9841, 0.706, 0.4967,...
            0.3515, 0.24713, 0.1753, 0.1259, 0.089, 0.0627, 0.04453, 0.03084, 0.021719];
        NIST_std_par.plane_radius = 170; % mm
        NIST_std_par.template_Sphere1_loc = [65.23, 27.07];
        NIST_std_par.dist_map = [0, 31.324, 58.893, 81.810, 95.898, 100.317,...
            95.991, 81.077, 58.860, 31.451, 35.380, 35.931, 73.124, 73.157]; % unit in mm
        NIST_std_par.ang_map = [0, 70.0850, 53.3538, 35.3270, 16.9426, -0.6324,...
            -18.7705, -36.9326, -56.4287, -72.7011, 33.6653, -34.0284, 14.9729, -16.4091]; % unit in degree
            
    case 'T2'
        NIST_std_par.true_value = [0.5813, 0.4035, 0.2781, 0.19094, 0.13327,...
            0.09689, 0.06407, 0.04642, 0.03197, 0.02256, 0.015813, 0.011237, 0.007911, 0.005592];
        NIST_std_par.plane_radius = 195; % mm
        NIST_std_par.template_Sphere1_loc = [64.8, 34.5];
        NIST_std_par.dist_map = [0, 31.547, 58.642, 80.697, 95.713, 99.812,...
            96.050, 80.936, 58.824, 30.036, 36.204, 35.882, 73.002, 72.552]; % unit in mm
        NIST_std_par.ang_map = [0, 71.4899, 54.0302, 35.7043, 17.5603, -0.1549,...
            -18.3830, -35.8010, -54.2831, -70.3638, 33.1854, -33.0284, 15.2130, -16.2861]; % unit in degree
end
NIST_std_par.intensity_range = [NIST_std_par.true_value(1)-0.2, NIST_std_par.true_value(1)+0.2];

%% convert mm to voxel based on map_size and fov
NIST_std_par.dist_map = NIST_std_par.dist_map/NIST_std_par.voxel_size; % convert to unit in voxel
NIST_std_par.sphere_radius = NIST_std_par.sphere_radius/NIST_std_par.voxel_size; % convert to unit in voxel
NIST_std_par.plane_radius = NIST_std_par.plane_radius/NIST_std_par.voxel_size; % convert to unit in voxel

%% get coordinates for all spheres
NIST_std_par.template_loc = zeros(NIST_std_par.num_spheres, 2);
NIST_std_par.template_loc(:, 1) = NIST_std_par.template_Sphere1_loc(1, 1) +...
    NIST_std_par.dist_map.*sind(NIST_std_par.ang_map); % x axis
NIST_std_par.template_loc(:, 2) = NIST_std_par.template_Sphere1_loc(1, 2) +...
    NIST_std_par.dist_map.*cosd(NIST_std_par.ang_map); % y axis 

end