function NIST_std_par = gen_NIST_std_par(map_type, map_size, fov)
%Input map parameters and return phantom parameters based on map type, map
%size and fov. 
% INPUT
%   map_type  type of map (T1 or T2)
%   map_size  size of map, for example, 128 means the map size is 128x128
%        fov  field of view                                  [mm]
%
% OUTPUT
%   phan_std_par  standard phantom parameters.
%
% 7/2019 Enlin Qian
% # Copyright of the Board of Trustees of Columbia University in the City of New York
NIST_std_par.num_spheres = 14;
NIST_std_par.voxel_size = fov/map_size; % unit in mm/voxel
NIST_std_par.sphere_radius = 8; % unit in mm
NIST_std_par.temp_rot_center = [map_size/2, map_size/2]; % unit in pixel

switch(map_type)
    case 'T1'
        NIST_std_par.true_value = [1.989, 1.454, 0.9841, 0.706, 0.4967,...
            0.3515, 0.24713, 0.1753, 0.1259, 0.089, 0.0627, 0.04453, 0.03084, 0.021719];
        NIST_std_par.plane_radius = 170; % mm
                 
    case 'T2'
        NIST_std_par.true_value = [0.5813, 0.4035, 0.2781, 0.19094, 0.13327,...
            0.09689, 0.06407, 0.04642, 0.03197, 0.02256, 0.015813, 0.011237, 0.007911, 0.005592];
        NIST_std_par.plane_radius = 195; % mm
end
NIST_std_par.intensity_range = [NIST_std_par.true_value(1)-0.2, NIST_std_par.true_value(1)+0.2];

%% convert mm to voxel based on map_size and fov
NIST_std_par.sphere_radius = NIST_std_par.sphere_radius/NIST_std_par.voxel_size; % convert to unit in voxel
NIST_std_par.plane_radius = NIST_std_par.plane_radius/NIST_std_par.voxel_size; % convert to unit in voxel

end