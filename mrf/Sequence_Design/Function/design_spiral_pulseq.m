function [k,dcf,t,ind,out,grad]=design_spiral_pulseq(fov,npix,arms,ksamp,...
    fname,gmax,smax,nucleus,acq_round2n,do_rot_file,balanced)
% This script designs a spiral with delayed acq for fast CSI
% INPUT
%        fov  field of view                                  [mm]
%       npix  #pixels (Cartesian resolution after gridding)
%       arms  #spatial interleaves                    (1)
%      ksamp  k-space sampling time                   (16)   [us]
%      fname  filename: if given->write_ak_wav 
%             if logical true: generate name          (false)
%       gmax  max gradient amp                        (40)   [mT/m]
%       smax  max slew rate                           (150)  [T/m/s]
%    nucleus  nucleus                                 ('1H')
%acq_round2n  Round up #acq pts to to 2^n             (true)
%do_rot_file  Write single static spiral into gradient 
%             waveform file & rotate via vap_phiXX.fdl (false)
%   balanced  Balancing gradient area                 (true)
%
% OUTPUT
%          k  k-space trajectory  [-0.5..0.5]
%        dcf  density compensation function (calculated with vornoi_area)
%          t  time  (nexc x recon-pts)
%        ind  index (2nd dim) for k-space points on spiral (excl ramps, etc)
%        out  output structure of wrt_wavs
%       grad  gradient waveform [T/m] with dt=10us 
%
% Created 7/2018 Rolf Schulte
% Modified 7/2019 Enlin Qian
if (nargin<1), help(mfilename); return; end


%% fixed parameters
fufa_gmax = 0.99;
fufa_smax = 0.99;
g_offset = 16;               % #grad pts offset for start of trajectory
gsamp = 10e-6;                % [s] gradient update time


%% input parameters
if ~exist('arms','var'),        arms = []; end    % #spatial interleaves
if isempty(arms),               arms = 1; end
if ~exist('ksamp','var'),       ksamp = []; end
if isempty(ksamp),              ksamp = 16; end   % [us] k-space dwell time
if ksamp<2, warning('dcs:ksamp','ksamp(=%g)<2',ksamp); input(''); end
if ~exist('fname','var'),       fname = []; end   % for waveform file
if ~exist('gmax','var'),        gmax = []; end
if isempty(gmax),               gmax = 40; end   % [mT/m]
if ~exist('smax','var'),        smax = []; end
if isempty(smax),               smax = 150; end  % [T/m/s]
if ~exist('nucleus','var'),     nucleus = []; end
if isempty(nucleus),            nucleus = '1H'; end
if ~exist('acq_round2n','var'), acq_round2n = []; end
if isempty(acq_round2n),        acq_round2n = true; end
if ~exist('do_rot_file','var'), do_rot_file = []; end
if isempty(do_rot_file),        do_rot_file = false; end
if ~exist('balanced','var'),    balanced = []; end
if isempty(balanced),           balanced = true; end


%% checks
if gmax<1, warning('dcs:gmax','gmax(=%g)<1',gmax); input(''); end
if fov<1, warning('dcs:fov','fov(=%g)<1',fov); input(''); end
dki = ksamp/gsamp*1d-6;
if dki<0.5
    error('ksamp: dki<0.5');
end
if (abs(dki-round(dki))>1d-10)
    if abs(dki-0.5)>1d-10
        error('ksamp (=%g us) not multiple of gsamp (=%g us)',...
            ksamp*1d6,gsamp*1d6);
    end
else
    dki = round(dki);
end
k_offset = g_offset/dki;


%% generate filename
if (islogical(fname) || isnumeric(fname))
    if fname
        if isnumeric(nucleus), error('Please enter string for nucleus'); end
        fname = sprintf('spiral_%s_fov%g_npix%g_arms%g_ksamp%g_gmax%g_smax%g', ...
            nucleus,round(fov),npix,arms,ksamp,round(gmax),round(smax));
        if do_rot_file, fname = [fname '_rofi']; end
    else
        fname = [];
    end
end


%% convert to standard SI units
fov = fov*1d-3;           % [mm] -> [m]
ksamp = ksamp*1d-6;       % [us] -> [s]
gmax = gmax*1d-3;         % [mT/m] -> [T/m]
res = fov/npix;           % [m] resolution


%% design spiral
rgamma = abs(gyrogamma('1h')/gyrogamma(nucleus));
fov_mns = fov/rgamma;
res_mns = res/rgamma;
fprintf('gamma ratio=%g\n',rgamma);
fprintf('scaled FOV (fov_mns) = %g [m]\n',fov_mns); 
fprintf('scaled res (res_mns) = %g [mm]\n',res_mns*1d3);
gmax_nyquist = 2*pi/(gyrogamma('1h')*ksamp*fov_mns);
fprintf('gmax_nyquist = %g [mT/m]\n',gmax_nyquist*1d3);
if (gmax_nyquist<gmax)
    fprintf('Attention: approaching sampling BW limited regime\n');
    if (ksamp>gsamp)
        fprintf('!!! Undersampling will occur: reduce ksamp !!!\n'); 
        input('press key to continue');
    end
end
pause(1);
[k1,g1,s1,t1,r1,theta1] = vds(fufa_smax*smax*1e2,fufa_gmax*gmax*1e2,gsamp, ...
   arms,[fov_mns*1e2,0],1/(2*res_mns*1e2));


%% calculate single interleave
gspir = 1e-2*g1;                % convert to SI unit [T/m]
kspir = gsamp*cumsum(gspir);    % [s*T/m]


%% rewinders
if balanced
    grewx = gradient_lobe(-real(kspir(1,end)),gsamp,...
        fufa_gmax*gmax/sqrt(2),fufa_smax*smax/sqrt(2),real(gspir(1,end)),false);
    grewy = gradient_lobe(-imag(kspir(1,end)),gsamp,...
        fufa_gmax*gmax/sqrt(2),fufa_smax*smax/sqrt(2),imag(gspir(1,end)),false);
    
    grew = zeros(1,max(length(grewx),length(grewy))+1);
    grew(1:length(grewx)) = 1*grewx;
    grew(1:length(grewy)) = grew(1:length(grewy))+1i*grewy;
else
    nrew = ceil(gmax/smax/gsamp/fufa_gmax);
    grew = linspace(real(gspir(1,end)),0,nrew) + ...
        1i*linspace(imag(gspir(1,end)),0,nrew);
end

% [LinearFill] = CheckSlewRate(gspir, grew, smax, ksamp);
% construct full (not yet rotated) gradient waveform
% gg = [zeros(1,g_offset),gspir,grew].';
ggrealadd1 = (real(grew(1))-real(gspir(end)))/4;
ggimagadd1 = (imag(grew(1))-imag(gspir(end)))/4;
ggrealadd2 = (real(grew(1))-real(gspir(end)))/4*2;
ggimagadd2 = (imag(grew(1))-imag(gspir(end)))/4*2;
ggrealadd3 = (real(grew(1))-real(gspir(end)))/4*3;
ggimagadd3 = (imag(grew(1))-imag(gspir(end)))/4*3;
ggreal = [zeros(1,g_offset),real(gspir),real(gspir(end))+ggrealadd1,real(gspir(end))+ggrealadd2,real(gspir(end))+ggrealadd3, real(grew)].';
ggimag = [zeros(1,g_offset),imag(gspir),imag(gspir(end))+ggimagadd1,imag(gspir(end))+ggimagadd2,imag(gspir(end))+ggimagadd3, imag(grew)].';
gg = complex(ggreal,ggimag);

ng = size(gg,1);
nk = round(length(kspir)/dki);


%% calculate time and index list
if acq_round2n
    acq_pts = 2^ceil(log2(ng/dki));
else
    acq_pts = ceil(length(g)/dki/2)*2;
end
if acq_pts>16384
    warning('dcs:acq','#sampling pts/exc (=%g)>16384: exceeds fidall limit',nk);
end
ind = false(arms,acq_pts);             % index for acquired data 
ind(:,(k_offset+(1:nk))) = true;
t = repmat((0:nk-1)*ksamp,[arms 1]);   % time of acq pts


%% rotate trajectory
phi = [];
if arms>1
    if do_rot_file
        phi = mod(-(0:arms-1)/arms*360,360);
    else
        gg = bsxfun(@times,gg,exp(1i*(0:arms-1)/arms*2*pi));
    end
end

%% k-space trajectory
if dki>1, k = kspir(1:dki:length(kspir))/2/pi*res_mns*gyrogamma('1H'); 
else,     k = kspir/2/pi*res_mns*gyrogamma('1H');
end
if abs(dki-0.5)<1d-10
    k = interp1((1:length(k)),k,(1:0.5:length(k)+0.5),'spline');
end
if arms>1
    % tmp = [];
    % for l=1:arms, tmp = [tmp , k*exp(1i*(l-1)/arms*2*pi)]; end
    % k = tmp;
    k = bsxfun(@times,k,exp(1i*(0:arms-1)/arms*2*pi).').';
    %k(arms)
    k = k(:).';
end
dcf = voronoi_area(k*npix);       % density compensation function


%% print info about waveform
gmax_act = max(max(abs(gg)));
smax_act = max(max(abs(diff(gg,1))))/gsamp;
fprintf('actual gmax  = %g [mT/m]\n',gmax_act*1d3);
if (gmax_act>gmax), warning('gmax exceeded'); end
fprintf('actual smax  = %g [T/m/s]\n',smax_act);
if (smax_act>smax), warning('smax exceeded'); end

desc1 = sprintf('Sequence details\n');
desc2 = sprintf('Acq BW = %g [kHz] (full)\n',1d-3/ksamp);
desc3 = sprintf('gsamp = %g [us]; ksamp = %g [us]\n',gsamp*1d6,ksamp*1d6);
desc4 = sprintf('g_pts = %gx%g; k_pts = %g; acq_pts/exc = %g\n',...
    size(gg,1),size(gg,2),size(k,2),acq_pts);
t_arm = t1(1,end); t_rew = (length(grew)+1)*gsamp;
desc5 = sprintf('t_arm = %g [ms]; t_rew = %g [ms]\n',t_arm*1d3,t_rew*1d3);
desc6 = sprintf('t_seq = %g [ms]\n',size(gg,1)*gsamp*1d3);
desc = [desc1 desc2 desc3 desc4 desc5 desc6];
fprintf('\n%s\n',desc);


%% checks
if any(size(dcf)~=size(k))
    warning('design_spiral:size','size(dcf)~=size(k): interpolating dcf'); 
    dcf = interp1(linspace(0,1,size(dcf,2)),dcf,linspace(0,1,size(k,2)));
end
if size(k)~=sum(ind(:)),     error('size(k)~=sum(ind(:))'); end
if prod(size(t))~=size(k,2), error('prod(size(t))~=size(k,2)'); end


%% export waveforms
if ~isempty(fname)
    fprintf('\nWriting out gradient waveforms + .mat\n');
    fprintf('fname = %s\n',fname);
    out = write_ak_wav([fname '.wav'],gg,1/ksamp,fov_mns,desc);
    save(fname,'out','k','dcf','t','ind','fov',...
        'npix','gsamp','ksamp','gmax','smax','nucleus',...
        't_rew','t_arm','g_offset','gmax_nyquist','rgamma','fov_mns','phi');
    if do_rot_file, write_fdl([fname '_phi.fdl'],phi,'phi'); end
else
    out.gmax = gmax;
    out.smax = smax;
    out.gdt = gsamp;
    out.kdt = ksamp;
    out.grad = gg;
    out.bw = 1/ksamp;
    out.fov = fov;
end


%% reshape + output grad
grad = zeros(2,size(gg,1),size(gg,2));
grad(1,:,:) = real(reshape(gg,[1 size(gg,1) size(gg,2)]));
grad(2,:,:) = imag(reshape(gg,[1 size(gg,1) size(gg,2)]));
