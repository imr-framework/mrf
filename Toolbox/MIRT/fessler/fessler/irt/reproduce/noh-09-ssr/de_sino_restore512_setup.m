%This is the code for DE-CT sinogram restoration
%Modified from spie02.m written by Jeff Fessler

if ~isvar('xtrue'), printm '512 xtrue'
	ddir = path_find_dir([filesep 'data']);
	tmp = fld_read([ddir filesep 'ncat,1024,slice500,0to4.fld']);
	jf1 = 0.2 * (tmp == 1) + 1.0 * (tmp == 2); % lung, body
	jf2 = 2.0 * (tmp == 3) + 2.0 * (tmp == 4); % spine, ribs
	xtrue(:,:,1) = downsample2(jf1, 2); % 512
	xtrue(:,:,2) = downsample2(jf2, 2);
	clear tmp jf ddir

	im pl 1 2
	im(1, xtrue)

	f.down = 4;
	xtrues(:,:,1) = downsample2(xtrue(:,:,1), f.down);
	xtrues(:,:,2) = downsample2(xtrue(:,:,2), f.down);
	im(2, xtrues)
prompt
end

if ~isvar('maskb'), printm 'mask'
	[nx ny] = size(xtrue(:,:,1));
	maskb = sum(xtrue, 3) > 0;
	maskb = imdilate(maskb, strel('disk', 5));
	maskb = logical(maskb);
	im(1, maskb)
prompt
end

if ~isvar('Gbf'), printm 'system models'
	igb = image_geom('nx', nx, 'ny', ny, 'dx', 0.1, 'mask', maskb);
	igs = igb.downsample(f.down);
	masks = igs.mask;
	im(2, masks)

	sgbf = sino_geom('ge1', 'units', 'cm'); % fan
	sgsf = sgbf.downsample(f.down);

	%Parallel-beam geometry
	sgbp = sino_geom('par', 'nb', 1024, 'na', 800, 'dr', sgbf.ds/2, ...
                    'offset_r', 0, 'orbit', 180); 
	sgsp = sgbp.downsample(f.down);

	Gbf = Gtomo2_dscmex(sgbf, igb);
	Gsf = Gtomo2_dscmex(sgsf, igs);
	Gbp = Gtomo2_dscmex(sgbp, igb);
	Gsp = Gtomo2_dscmex(sgsp, igs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f table and true images %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isvar('ftab'), printm 'ftab'
	Iodine_sel = 0;
	if Iodine_sel == 0
%		ftab = de_ftab_build({'ps1'}); %DE-Table

		xrs = xray_read_spectra('ps1');
		sls = de_ftab_sls('max', [40 40], 'n', [120 120]);
		mas = xray_read_mac({'water', 'bone'});
		ftab = de_ftab(xrs, mas, 'sls', sls, 'ctype', 'pre10', ...
			'ftype','exp', 'fit_args', {'kev', 10:5:140, 'mac', []});
	elseif Iodine_sel == 1
    	ftab = de_ftab_build({'ps1t'}); %TE-Table
    	%rhoI = 4.930; %Pure Iodine density (from Hubbell&Seltzer)
    	rhoI1 = 0.08; rhoI2 = 0.01;%0.05 diluted Iodine density
    	ell1 = [-10 0 0.3 0.5 0 rhoI1];
    	ell2 = [0 0 0.5 0.3 0 rhoI1];
    	ell3 = [0 5 0.5 0.5 0 rhoI2];
    	xtrue(:,:,3) = ellipse_im(ig,ell1) + ellipse_im(ig,ell2) + ellipse_im(ig,ell3);
	end %iodine with soft tissue

	%true density map display
	im clf
	c.soft = [0 1.2]; c.bone = [0 2.2]; c.dens = [0 5.2];
	im(231, xtrue(:,:,1), 'SoftTissue Density', c.soft), cbar([0 1])
	im(232, xtrue(:,:,2), 'Bone Density', c.bone), cbar([0 2])
	if Iodine_sel == 1, im(233, xtrue(:,:,3), 'Iodine Density', c.dens), cbar([0 5.2]); end
	im(234, sum(xtrue,3), 'Density Map', c.dens), cbar([0 2])
	im(235, (sum(xtrue,3)>0.5) + 20*double(maskb), 'Reconstruction Support')
end

%For TMI
%figure;im(xtrue(:,:,1),c.soft), cbar([0 1]);
%figure;im(xtrue(:,:,2),c.bone), cbar([0 2]);
%figure;im(xtrue(:,:,1) + xtrue(:,:,3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ideal and Noisy Measurement % y_{mi}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isvar('strue'), printm 'strue'
	%Noiseless sinogram and f
	tic; %noiseless forward projection
	strue = Gbf * xtrue;
	printm('forward projection time %0.3g', toc)
end

if ~isvar('ftrue'), printm 'ftrue'
%	ftrue = de_ftab_fm(strue, ftab.mac.mac, ftab.xray.Ide); %noiseless f
	ftrue = ftab.fit.fmfun(strue);
	im(ftrue)
end

if ~isvar('ymi'), printm 'ymi'
	%Noisy measurement
	ybi = zeros(sgbf.nb, sgbf.na, 2); %bar{y}_{mi}
	ybi(:,:,1) = ftab.xrs.I(1) * exp(-ftrue(:,:,1)); %true mean
	ybi(:,:,2) = ftab.xrs.I(2) * exp(-ftrue(:,:,2));

	%Dose amount
	%f.scale_high = inf; f.scale_low = inf;%for noiseless profile
	f.scale_high = 5*10^5 / ybi(1,1,2); %high dose
	f.scale_low = 2*10^5 / ybi(1,1,2); %low dose

	if isinf(f.scale_high) || isinf(f.scale_low)
		ymi = ybi; se_ymi = ybi;	% noiseless
	else
		ymi = poisson(f.scale_low*ybi, 0) / f.scale_low; % corrupted by poisson %
		se_ymi = poisson(f.scale_high*ybi, 0) / f.scale_high;
	end

	%Borrow neighbor(s) for any log(0) values
	off = 1; %DECT
	while any(ymi(:) == 0)
		ii = find(ymi(:) == 0);
		printm('fixing %d zeros in ymi with %d', length(ii), off)
		ymi(ii) = max(ymi(ii+off), ymi(ii-off));
		off = off + 1;
	end, clear off ii
	off = 1; %SECT
	while any(se_ymi(:) == 0)
		ii = find(se_ymi(:) == 0);
		printm('fixing %d zeros in se_ymi with %d', length(ii), off)
		se_ymi(ii) = max(se_ymi(ii+off), se_ymi(ii-off));
		off = off + 1;
	end, clear off ii
end %measurements look like sinogram


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry Tx & downsampling %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
geo_sel = 1;
if geo_sel == 1 %fan-beam
	ymipd(:,:,1) = downsample2(rebin_fan2par(ymi(:,:,1), sgbf, sgbp),f.down); %fan -> par
	ymipd(:,:,2) = downsample2(rebin_fan2par(ymi(:,:,2), sgbf, sgbp),f.down); %[n.bp n.ap]
	se_ymipd(:,:,1) = downsample2(rebin_fan2par(se_ymi(:,:,1), sgbf, sgbp),f.down); 
	se_ymipd(:,:,2) = downsample2(rebin_fan2par(se_ymi(:,:,2), sgbf, sgbp),f.down);
	for pp=1:size(xtrue,3)
		strue_pd(:,:,pp) = downsample2(rebin_fan2par(strue(:,:,pp), sgbf, sgbp),f.down);
	end
elseif geo_sel == 0 %parallel-beam
	ymipd(:,:,1) = downsample2(ymi(:,:,1),down);  %downsampling
	ymipd(:,:,2) = downsample2(ymi(:,:,2),down); %[n.bp n.ap]
	se_ymipd(:,:,1) = downsample2(se_ymi(:,:,1),down); 
	se_ymipd(:,:,2) = downsample2(se_ymi(:,:,2),down); 
	for pp=1:size(xtrue,3)
		strue_pd(:,:,pp) = downsample2(strue(:,:,pp),down);
	end
end
