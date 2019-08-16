function out = write_ak_wav(fname,grad,bw,fov,desc,n_kpts)
%WRITE_AK_WAV  Saves waveforms to external file using Stanford format
% out = write_ak_wav(fname,grad,bw,fov,desc,n_kpts,rfs_fastfid,ak_dynwf_mode)
% INPUT (SI Units)
%        fname  File name of output file                      [string]
%         grad  Gradient waveforms with 4us time resolution   [T/m]
%               [#pts/interleave,#interleaves,#groups]
%               with: #groups = 2 for 2d-imaging, =3 for 3d-imaging
%               if #groups=1 and complex grad -> 2D
%           bw  Full sampling bandwidth; scalar               [Hz]
%          fov  Field-of-view relative to 1H; scalar          [m]
%          des  Description string (opt; up to 254 chars)
%       n_kpts  Number of acq points (opuser1) (opt)
%               default = ceil(N.gpts*gdt/kdt)
%
% OUTPUT (GE Scanner Units):
% out               Output structure
% out.descr         Description
% out.N.gpts        # input gradient pts/interleave
% out.N.kpts        # readout pts
% out.N.groups      # groups (=2; real and imaginary)
% out.N.intl        # interleaves
% out.N.params      # parameters (=14)
% out.parms=params  Header file parameters
% out.wave=wave     Output waveform
%
% Difference to previous versions: save indices to worst (3D) trajectories 
% with >1 interleaves into header for heating calculations (pgen on host)
%
% 2007       Adam@mrsrl.stanford.edu 
% 7/2018     Rolf Schulte
if (nargin<1), help(mfilename); return; end

if ~exist('n_kpts','var'), n_kpts = []; end


%% check parameters
if length(bw)~=1,         warning('length(bw)~=1'); end
if length(fov)~=1,        warning('length(fov)~=1'); end
if ~isempty(n_kpts)
    if length(n_kpts)~=1, warning('length(n_kpts)~=1'); end
end
if ~ischar(desc),         warning('~ischar(desc)'); end


%% FIXED PARAMETERS
gdt = 4;                  % [us] Gradient Sampling time; scalar
grad_type = 0;


%% check trajectory
switch length(size(grad))
    case 2, 
        if ~isreal(grad),
            tmp(:,:,1) = real(grad);
            tmp(:,:,2) = imag(grad);
            grad = tmp; clear tmp;
        end
    case 3,
    otherwise, error('length(grad)');
end
if ~isreal(grad), error('~isreal(grad)'); end
% if (size(grad,1)<size(grad,2)), 
%     warning('write_ak_wav:size','(size(grad,1)<size(grad,2))'); 
% end
if (size(grad,1)<size(grad,3)),
    warning('write_ak_wav:size','size(grad,1)<size(grad,3)'); 
end
if any(any(abs(grad(1,:,:))>1d-10)), error('grad must start with 0'); end
if any(any(abs(grad(end,:,:))>1d-10)), error('grad must end with 0'); end
if isodd(size(grad,1)),
    warning('size(grad,1)(=%g) is odd; adding zero to end',size(grad,1));
    grad = [grad ; zeros(1,size(grad,2),size(grad,3))];
end


%% DERIVED PARAMETERS: Plus conversion to GE scanner units
grad = grad*1d2;                              % [T/m] -> [G/cm]
gmax = max(abs(grad(:)));
smax = max(max(max(abs(diff(grad,[],1)))))/(gdt*1d-3);
fov = fov*1d2;                                % [m] -> [cm]
kdt = 1d6/bw;                                 % Sampling time; scalar [us]
if abs(ceil(kdt)-kdt)>1d-10, 
    error('kdt (=%g[us]) must be multiple of 1[us]',kdt); 
end
[N.gpts,N.intl,N.groups]=size(grad);          % [#pts/intl,#intl,#groups]
if isempty(n_kpts),
    n_kpts = ceil(N.gpts*gdt/kdt);            % # readout pts
end
if n_kpts<64, 
    warning('write_ak_wav:kpts',...
        'n_kpts (=%g) < 64; setting n_kpts=64 (fidall minimum)',n_kpts); 
    n_kpts = 64;
end
N.kpts = n_kpts;                              % # readout pts
if N.gpts > 32766,
    warning('write_ak_wav:gpts',...
        'n_gpts (=%g) > 32766; exceeding maximum of fidall',N.gpts); 
end
if N.intl > 16382,
    warning('write_ak_wav:intl',...
        'n_intl (=%g) > 16382; exceeding maximum of fidall',N.intl); 
end

fprintf('gmax = %g [mT/m] = %g [G/cm]\n',gmax*10,gmax);
fprintf('smax = %g [T/m/s] = %g [G/cm/ms]\n',smax*10,smax);


%% CREATE HEADER INFORMATION
if ~exist('desc','var'), desc = ' '; end
if isempty(desc), desc = ' '; end
if ~ischar(desc), error('~ischar(desc)'); end
% make length(des) = 256
if length(desc)>255,
    desc = desc(1:256);
    warning('write_ak_wav:desc','length(desc)>256 -> truncating');
else
    desc(length(desc):256) = ' ';
end
desc(254:256) = sprintf('\n\f\n');


%% Determine worst (sum-of-squares) sub-waveform for pulsegen on host
wf4pg = zeros(1,3);   % index to waveform [Gx,Gy,Gz]
% following c-convention 0=first, etc
if N.intl>1,
    tmp = sum(grad.^2,1);
    for l1=1:N.groups,
        sos_val = 0;
        for l2=1:N.intl,
            if tmp(1,l2,l1) > sos_val*1.001
                % slightly favour old solutions to catch errors in case of
                % identical waveforms (same in rth_pgen.c)
                sos_val = tmp(1,l2,l1);
                wf4pg(1,l1) = l2-1;
            end
        end
    end
    fprintf('Worst trajectories: ak_wf4pg_x=%g, y=%g, z=%g\n', wf4pg);
end


%% header parameters
params = [grad_type fov N.intl str2num(sprintf('%4.6f',gmax)) N.gpts ...
    gdt N.kpts kdt 0 0 0 wf4pg];
% hw_limits removed
N.params = length(params);            % # parameters


%% WRITE HEADER + WAVEFORMS
% Int16 rounding + check that waveform is even
% wave=2*round((2^14-1)/gmax*grad);
wave=2*int16((2^14-1)/gmax*grad);

% Stop criterion: Last pt. of waveform must be odd (e.g.=1)
wave(end,:,:)=1;
% File IO: Open + Write + Close
fid = fopen(fname,'w','b');
if (fid==-1)
  error(['Could not write file: ',fname]);
end
fwrite(fid, desc, 'char');
fwrite(fid, N.gpts, 'uint16');
fwrite(fid, N.groups, 'uint16');
fwrite(fid, N.intl*ones(1,N.groups), 'uint16');
fwrite(fid, N.params, 'uint16');
fwrite(fid, params, 'float64');
fwrite(fid, wave, 'int16');
fclose(fid);


%% Create come output/feedback for control purposes
out.N = N; 
out.parms = params; 
out.wave = wave;
out.gmax = gmax;
out.smax = smax;
out.gdt = gdt;
out.kdt = kdt;
out.grad = grad;
out.bw = bw;
out.fov = fov;
out.desc = desc;

