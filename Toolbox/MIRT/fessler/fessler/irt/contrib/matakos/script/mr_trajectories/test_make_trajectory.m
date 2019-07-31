% Script for testing EPI and spiral trajectories for the scanner.
fov = 24;
nshot = 1;
etl = 64;
dt = 4e-6;
TE = 0.03; % 30ms echo-time

npix = etl*nshot;
N = [npix npix];

wi_epi = 1/fov^2;


%% Create EPI trajectories

pre = '../../data/ge/trajectories/epi_ss_int/epi_ss';

% Interleaved single shot EPI
[ks1,om1,~,~,ts1,idx1,rlen1] = make_epi(fov,npix,etl,dt,'prefix',pre,...
	'type','int','smpr',0,'out',1,'dirx',-1,'diry',-1);

% Keep points on grid
ks1 = ks1(idx1,:);
om1 = om1(idx1,:);
ts1 = ts1(idx1);

% Save trajectory params
in_traj.trajectory = {ks1 om1 wi_epi};
in_traj.ti = ts1;
in_traj.nufft = {N  [6 6]  2*N  N/2  'table' 2^10 'minmax:kb'};
in_traj.idx = idx1;
in_traj.frsize = rlen1;

save ../../data/ge/trajectories/epi_ss_int/input_traj_epi_ss_int in_traj


pre = '../../data/ge/trajectories/epi_ss_std/epi_ss';

% Standard single shot EPI
[ks2,om2,~,~,ts2,idx2,rlen2] = make_epi(fov,npix,etl,dt,'prefix',pre,...
	'type','std','smpr',0,'out',1,'dirx',-1,'diry',-1);

% Keep points on grid
ks2 = ks2(idx2,:);
om2 = om2(idx2,:);
ts2 = ts2(idx2);

% Save trajectory params
in_traj.trajectory = {ks2 om2 wi_epi};
in_traj.ti = ts2;
in_traj.nufft = {N  [6 6]  2*N  N/2  'table' 2^10 'minmax:kb'};
in_traj.idx = idx2;
in_traj.frsize = rlen2;

save ../../data/ge/trajectories/epi_ss_std/input_traj_epi_ss_std in_traj




% Standard two shot EPI
pre = '../../data/ge/trajectories/epi_s2_std/epi_s2';
[~,~,gr3,sl3,~,idx3,rlen3] = make_epi(fov,npix,etl/2,dt,'prefix',pre,...
	'type','std','smpr',0,'out',1,'offset',0,'dir',-1);

ks3 = [];
om3 = [];
ts3 = [];

parfor ii=1:2

	[ks,om,gr,sl,ts,idx] = make_epi(fov,npix,etl/2,dt,'prefix',pre,...
		'type','std','smpr',0,'out',0,'offset',ii-1,'dir',-1);
	
	ks = ks(idx,:);
	om = om(idx,:);
	ts = ts(idx);
	
% 	figure(3), plot(ks(:,1), ks(:,2), '.')
% 	title 'Standard 2 shot EPI', axis square;
	
	ks3 = [ks3;ks]; 
	om3 = [om3;om]; 
	ts3 = [ts3;ts]; 
end

% Save trajectory params
in_traj.trajectory = {ks3 om3 wi_epi};
in_traj.ti = ts3;
in_traj.nufft = {N  [6 6]  2*N  N/2  'table' 2^10 'minmax:kb'};
in_traj.idx = idx3;
in_traj.frsize = rlen3;

save ../../data/ge/trajectories/epi_s2_std/input_traj_epi_s2_std in_traj



% Standard four shot EPI
pre = '../../data/ge/trajectories/epi_s4_std/epi_s4';
[~,~,gr4,sl4,~,idx4,rlen4] = make_epi(fov,npix,etl/4,dt,'prefix',pre,...
	'type','std','smpr',0,'out',1,'offset',0,'dir',-1);

ks4 = [];
om4 = [];
ts4 = [];

parfor ii=1:4

	[ks,om,gr,sl,ts,idx] = make_epi(fov,npix,etl/4,dt,'prefix',pre,...
		'type','std','smpr',0,'out',0,'offset',ii-1,'dir',-1);
	
	ks = ks(idx,:);
	om = om(idx,:);
	ts = ts(idx);
	
% 	figure(4), plot(ks(:,1), ks(:,2), '.')
% 	title 'Standard 4 shot EPI', axis square;

	ks4 = [ks4;ks]; 
	om4 = [om4;om]; 
	ts4 = [ts4;ts]; 
end

% Save trajectory params
in_traj.trajectory = {ks4 om4 wi_epi};
in_traj.ti = ts4;
in_traj.nufft = {N  [6 6]  2*N  N/2  'table' 2^10 'minmax:kb'};
in_traj.idx = idx4;
in_traj.frsize = rlen4;

save ../../data/ge/trajectories/epi_s4_std/input_traj_epi_s4_std in_traj



% Standard eight shot EPI
pre = '../../data/ge/trajectories/epi_s8_std/epi_s8';
[~,~,gr5,sl5,~,idx5,rlen5] = make_epi(fov,npix,etl/8,dt,'prefix',pre,...
	'type','std','smpr',0,'out',1,'offset',0,'dir',-1);

ks5 = [];
om5 = [];
ts5 = [];

parfor ii=1:8

	[ks,om,gr,sl,ts,idx] = make_epi(fov,npix,etl/8,dt,'prefix',pre,...
		'type','std','smpr',0,'out',0,'offset',ii-1,'dir',-1);
	
	ks = ks(idx,:);
	om = om(idx,:);
	ts = ts(idx);
	
% 	figure(4), plot(ks(:,1), ks(:,2), '.')
% 	title 'Standard 4 shot EPI', axis square;

	ks5 = [ks5;ks]; 
	om5 = [om5;om]; 
	ts5 = [ts5;ts]; 
end

% Save trajectory params
in_traj.trajectory = {ks5 om5 wi_epi};
in_traj.ti = ts5;
in_traj.nufft = {N  [6 6]  2*N  N/2  'table' 2^10 'minmax:kb'};
in_traj.idx = idx5;
in_traj.frsize = rlen5;

save ../../data/ge/trajectories/epi_s8_std/input_traj_epi_s8_std in_traj


% Standard sixteen shot EPI
pre = '../../data/ge/trajectories/epi_s16_std/epi_s16';
[~,~,gr6,sl6,~,idx6,rlen6] = make_epi(fov,npix,etl/16,dt,'prefix',pre,...
	'type','std','smpr',0,'out',1,'offset',0,'dir',-1);

ks6 = [];
om6 = [];
ts6 = [];

parfor ii=1:16

	[ks,om,gr,sl,ts,idx] = make_epi(fov,npix,etl/16,dt,'prefix',pre,...
		'type','std','smpr',0,'out',0,'offset',ii-1,'dir',-1);
	
	ks = ks(idx,:);
	om = om(idx,:);
	ts = ts(idx);
	
% 	figure(4), plot(ks(:,1), ks(:,2), '.')
% 	title 'Standard 4 shot EPI', axis square;

	ks6 = [ks6;ks];
	om6 = [om6;om];
	ts6 = [ts6;ts];
end

% Save trajectory params
in_traj.trajectory = {ks6 om6 wi_epi};
in_traj.ti = ts6;
in_traj.nufft = {N  [6 6]  2*N  N/2  'table' 2^10 'minmax:kb'};
in_traj.idx = idx6;
in_traj.frsize = rlen6;

save ../../data/ge/trajectories/epi_s16_std/input_traj_epi_s16_std in_traj


% Standard sixty-four shot EPI
pre = '../../data/ge/trajectories/epi_s64_std/epi_s64';
[~,~,gr7,sl7,~,idx7,rlen7] = make_epi(fov,npix,etl/64,dt,'prefix',pre,...
	'type','std','smpr',0,'out',1,'offset',0,'dir',-1);

ks7 = [];
om7 = [];
ts7 = [];

parfor ii=1:64

	[ks,om,gr,sl,ts,idx] = make_epi(fov,npix,etl/64,dt,'prefix',pre,...
		'type','std','smpr',0,'out',0,'offset',ii-1,'dir',-1);
	
	ks = ks(idx,:);
	om = om(idx,:);
	ts = ts(idx);
	
% 	figure(4), plot(ks(:,1), ks(:,2), '.')
% 	title 'Standard 4 shot EPI', axis square;

	ks7 = [ks7;ks];
	om7 = [om7;om];
	ts7 = [ts7;ts];
end

% Save trajectory params
in_traj.trajectory = {ks7 om7 wi_epi};
in_traj.ti = ts7;
in_traj.nufft = {N  [6 6]  2*N  N/2  'table' 2^10 'minmax:kb'};
in_traj.idx = idx7;
in_traj.frsize = rlen7;

save ../../data/ge/trajectories/epi_s64_std/input_traj_epi_s64_std in_traj


pre = '../../data/ge/trajectories/epi_ssf_int/epi_ssf1';

% Interleaved "fake" single shot EPI
[ks8,om8,gr8,sl8,ts8,idx8,rlen8] = make_epi(fov,npix,etl/2,dt,'prefix',pre,...
	'type','usmp','smpr',0,'out',1,'offset',0);

% Keep points on grid
ks8 = ks8(idx8,:);
om8 = om8(idx8,:);
ts8 = ts8(idx8);

% Save trajectory params
in_traj_f1.trajectory = {ks8 om8 wi_epi};
in_traj_f1.ti = ts8;
in_traj_f1.nufft = {N  [6 6]  2*N  N/2  'table' 2^10 'minmax:kb'};
in_traj_f1.idx = idx8;
in_traj_f1.frsize = rlen8;

pre = '../../data/ge/trajectories/epi_ssf_int/epi_ssf2';

% Interleaved "fake" single shot EPI
[ks9,om9,gr9,sl9,ts9,idx9,rlen9] = make_epi(fov,npix,etl/2,dt,'prefix',pre,...
	'type','usmp','smpr',0,'out',1,'offset',1);

% Keep points on grid
ks9 = ks9(idx9,:);
om9 = om9(idx9,:);
ts9 = ts9(idx9);

% Save trajectory params
in_traj_f2.trajectory = {ks9 om9 wi_epi};
in_traj_f2.ti = ts9;
in_traj_f2.nufft = {N  [6 6]  2*N  N/2  'table' 2^10 'minmax:kb'};
in_traj_f2.idx = idx9;
in_traj_f2.frsize = rlen9;

save ../../data/ge/trajectories/epi_ssf_int/input_traj_epi_ssf_int ...
	in_traj_f1 in_traj_f2


% Plot kspace samples
figure(1), plot(ks1(:,1), ks1(:,2), '.')
title 'Interleaved single shot EPI', axis square;
figure(2), plot(ks2(:,1), ks2(:,2), '.')
title 'Standard single shot EPI', axis square;
figure(3), plot(ks3(:,1), ks3(:,2), '.')
title 'Standard 2 shot EPI', axis square;
figure(4), plot(ks4(:,1), ks4(:,2), '.')
title 'Standard 4 shot EPI', axis square;
figure(5), plot(ks5(:,1), ks5(:,2), '.')
title 'Standard 8 shot EPI', axis square;
figure(6), plot(ks6(:,1), ks6(:,2), '.')
title 'Standard 16 shot EPI', axis square;
figure(7), plot(ks7(:,1), ks7(:,2), '.')
title 'Standard 64 shot EPI', axis square;
figure(8), plot(ks8(:,1), ks8(:,2), 'b.',ks9(:,1), ks9(:,2), 'r.')
title 'Interleaved fake single shot EPI', axis square;

% Plot gradient waveforms
% figure(5),plot(1:length(gr1),gr1(:,1),1:length(gr1),gr1(:,2))
% figure(6),plot(1:length(gr2),gr2(:,1),1:length(gr2),gr2(:,2))
% figure(7),plot(1:length(gr3),gr3(:,1),1:length(gr3),gr3(:,2))
% figure(8),plot(1:length(gr4),gr4(:,1),1:length(gr4),gr4(:,2))

% Plot slew rate waveforms
% figure(9),plot(1:length(sl1),sl1(:,1),1:length(sl1),sl1(:,2))
% figure(10),plot(1:length(sl2),sl2(:,1),1:length(sl2),sl2(:,2))
% figure(11),plot(1:length(sl3),sl3(:,1),1:length(sl3),sl3(:,2))
% figure(12),plot(1:length(sl4),sl4(:,1),1:length(sl4),sl4(:,2))



%% Create spiral trajectories

pre = '../../data/ge/trajectories/spr_ss_int/spr_ss';

% Interleaved single shot spiral
[ks1,om1,gr1,sl1,ts1,idx1,rlen1] = make_spiral(fov,npix,dt,'prefix',pre,...
	'type','int','out',1);

% Keep points on grid
ks1 = ks1(idx1,:);
om1 = om1(idx1,:);
ts1 = ts1(idx1);
wi1 = mri_density_comp(ks1, 'voronoi');

% Save trajectory params
in_traj.trajectory = {ks1 om1 wi1};
in_traj.ti = ts1;
in_traj.nufft = {N  [6 6]  2*N  N/2  'table' 2^10 'minmax:kb'};
in_traj.idx = idx1;
in_traj.frsize = rlen2;

save ../../data/ge/trajectories/spr_ss_int/input_traj_spr_ss_int in_traj



pre = '../../data/ge/trajectories/spr_ss_std/spr_ss';

% Standard single shot spiral
[ks2,om2,gr2,sl2,ts2,idx2,rlen2] = make_spiral(fov,npix,dt,'prefix',pre,...
	'type','std','out',1);

% Keep points on grid
ks2 = ks2(idx2,:);
om2 = om2(idx2,:);
ts2 = ts2(idx2);
wi2 = mri_density_comp(ks2, 'voronoi');

% Save trajectory params
in_traj.trajectory = {ks2 om2 wi2};
in_traj.ti = ts2;
in_traj.nufft = {N  [6 6]  2*N  N/2  'table' 2^10 'minmax:kb'};
in_traj.idx = idx2;
in_traj.frsize = rlen2;

save ../../data/ge/trajectories/spr_ss_std/input_traj_spr_ss_std in_traj


pre1 = '../../data/ge/trajectories/spr_s2_std/spr_s';

ks3 = [];
om3 = [];
ts3 = [];
idx3 = cell(2,1);

for ii=1:2
	pre = [pre1 num2str(ii)];
	
	[ks,om,gr,sl,ts,idx,rlen3] = make_spiral(fov,npix,dt,'prefix',pre,'shot',2,...
		'type','std','out',1,'rot',2*(ii-1));
	
	ks = ks(idx,:);
	om = om(idx,:);
	ts = ts(idx);

	ks3 = [ks3;ks];
	om3 = [om3;om];
	ts3 = [ts3;ts];
	idx3{ii} = idx;
end

wi3 = mri_density_comp(ks3, 'voronoi');

in_traj.trajectory = {ks3 om3 wi3};
in_traj.ti = ts3;
in_traj.nufft = {N  [6 6]  2*N  N/2  'table' 2^10 'minmax:kb'};
in_traj.idx = idx3;
in_traj.frsize = rlen3;

save ../../data/ge/trajectories/spr_s2_std/input_traj_spr_s2_std in_traj


pre1 = '../../data/ge/trajectories/spr_s4_std/spr_s';

ks4 = [];
om4 = [];
ts4 = [];
idx4 = cell(2,1);

for ii=1:4
	pre = [pre1 num2str(ii)];
	
	[ks,om,gr,sl,ts,idx,rlen4] = make_spiral(fov,npix,dt,'prefix',pre,'shot',4,...
		'type','std','out',1,'rot',ii-1);
	
	ks = ks(idx,:);
	om = om(idx,:);
	ts = ts(idx);

	ks4 = [ks4;ks];
	om4 = [om4;om];
	ts4 = [ts4;ts];
	idx4{ii} = idx;
end

wi4 = mri_density_comp(ks4, 'voronoi');

in_traj.trajectory = {ks4 om4 wi4};
in_traj.ti = ts4;
in_traj.nufft = {N  [6 6]  2*N  N/2  'table' 2^10 'minmax:kb'};
in_traj.idx = idx4;
in_traj.frsize = rlen4;

save ../../data/ge/trajectories/spr_s4_std/input_traj_spr_s4_std in_traj



% Plot kspace samples
figure(9), plot(ks1(:,1), ks1(:,2), '.')
title 'Interleaved single shot spiral', axis square;
figure(10), plot(ks2(:,1), ks2(:,2), '.')
title 'Standard single shot spiral', axis square;
figure(11), plot(ks3(:,1), ks3(:,2), '.')
title 'Standard 2 shot spiral', axis square;
figure(12), plot(ks4(:,1), ks4(:,2), '.')
title 'Standard 4 shot spiral', axis square;



