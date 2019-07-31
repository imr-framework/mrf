%This is the code for DE-CT sinogram restoration
%Modified from spie02.m written by Jeff Fessler

%%%%%%%%%%%%
%Load data %
%%%%%%%%%%%%
Iodine_sel = 1;
geo_sel = 0;
load xtrue512.mat; %512 by 512
if Iodine_sel == 0
    ftab = de_ftab_build({'ps1'}); %DE-Table
elseif Iodine_sel == 1
    ftab = de_ftab_build({'ps1t'}); %TE-Table
    rhoI = 4.930; %Iodine density (from Hubbell&Seltzer)
    xrang = (241:260);yrang = xrang;
    xtrue(xrang,yrang,3) = rhoI;
end %iodine with soft tissue
[n.x n.y] = size(xtrue(:,:,1));
%mask
mask = sum(xtrue, 3) > 0;
mask = imdilate(mask, strel('disk', 5));
mask = logical(mask);
%true density map display
im clf
c.soft = [0 1.2];c.bone = [0 2.2];c.dens = [0 2.2];
im(221, xtrue(:,:,1), 'SoftTissue Density', c.soft), cbar([0 1])
im(222, xtrue(:,:,2), 'Bone Density', c.bone), cbar([0 2])
im(223, sum(xtrue,3), 'Density Map', c.dens), cbar([0 2])
im(224, (sum(xtrue,3)>0.5) + 20*double(mask), 'Reconstruction Support')

%%%%%%%%%%%%%%%%
%System matrix %
%%%%%%%%%%%%%%%%
if geo_sel == 1 %fan-beam geometry
    f.dx = 0.1;%in cm
    n.bf = 888;n.af = 984;ds = 0.1;
    n.bp = 1024;n.ap = 800;dr = ds/2;
    %Image geometry
    ig = image_geom('nx', n.x, 'ny', n.y, 'dx', f.dx, 'mask', mask);
    %Fan-beam geometry
    sgf = sino_geom('fan', 'nb', n.bf, 'na', n.af, 'ds', ds, 'dsd', 949,...
                    'dod', 408, 'offset_s', 1.25); 
    %Parallel-beam geometry
    sgp = sino_geom('par', 'nb', n.bp, 'na', n.ap, 'dr', dr, ...
                    'offset_r', 0, 'orbit', 180); 
    %G = Gtomo2_strip(sg, ig); %small CT
    G = Gtomo2_dscmex(sgf, ig); %full CT (Fan-beam)
    gi = reshape(sum(G'), [n.bf n.af]); %simple projection

elseif geo_sel == 0 %parallel-beam geometry
    f.dx = 0.1; %in cm
    n.bf = 1024;n.af = 800;ds = 0.1;
    n.bp = n.bf;n.ap = n.af;dr = ds/2;
    %Image geometry
    ig = image_geom('nx', n.x, 'ny', n.y, 'dx', f.dx, 'mask', mask);
    %Parallel-beam geometry
    sgp = sino_geom('par', 'nb', n.bp, 'na', n.ap, 'dr', dr, ...
                    'offset_r', 0, 'orbit', 180); 
    G = Gtomo2_dscmex(sgp, ig); %full CT (Para-beam)
    gi = reshape(sum(G'), [n.bf n.af]); %simple projection
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ideal and Noisy Measurement % y_{mi}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Noiseless sinogram and f
%
xtrueL = size(xtrue,3);
tmp = reshape(xtrue, n.x*n.y, xtrueL);
tic; %noiseless forward projection
strue = reshape(G * tmp(mask,:), [n.bf n.af xtrueL]);
printm('forward projection time %0.3g', toc)
file_name = sprintf('test_setup_I%d.mat',Iodine_sel);
if Iodine_sel == 0 %understand this part !!
    ftrue = ftab.fm_fun(ftab, {strue(:,:,1),strue(:,:,2)}); %noiseless f
elseif Iodine_sel == 1
    ftrue = ftab.fm_fun(ftab, {strue(:,:,1),strue(:,:,2),strue(:,:,3)});
end
save(file_name);
