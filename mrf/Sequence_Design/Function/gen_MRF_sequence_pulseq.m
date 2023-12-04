function [kshot, dcf, ind, TR_all, FA_all, TE_all] = gen_MRF_sequence_pulseq()
% This script uses pulseq to generate a sequence proposed in Jiang's paper[1]
% INPUT
%       None  User does not need to input any parameter. All parameters for
%             the sueqnece are pre set. If user wants to customize their own
%             MRF sequence, please modify the double for loops starting at
%             line 98 to add different blocks.
%
% OUTPUT
%      kshot  Spiral trajector in k-space
%        dcf  Density compensation function
%        ind  Index (2nd dim) for k-space points on spiral
%     TR_all  All TRs used for the sequence.
%     FA_all  All FAs used for the sequence.
%     TE_all  All TEs used for the sequence.

% 7/2019 Enlin Qian
% # Copyright of the Board of Trustees of Columbia University in the City of New York
% [1] Jiang Y, Ma D, Seiberlich N, Gulani V, Griswold M. MR fingerprinting using
% fast imaging with steady state precession (FISP) with spiral readout. Magn Reson Med.
% 2015;74(6):spcone-spcone. doi:10.1002/mrm.26048

%% Set system limits
gamma = 42576000; % in Hz %Determined from Pulseq - do not change,%Hz/T
Gmax = 32; % mT/m
SRmax = 130;% T/m/s
system = mr.opts('MaxGrad',Gmax,'GradUnit','mT/m',...
    'MaxSlew',SRmax,'SlewUnit','T/m/s', 'gradRasterTime', 10e-6, 'rfRingdownTime',10e-6);
seq=mr.Sequence(system);
total_acq = 1000; % total 1000 time points

%% Load data
fov=225e-3;
Nx=256;% Ny=256;
sliceThickness=5e-3;
load('method_orig_147.mat')
TR_all = method.VariableTR*1e-3; % converted to s (make min TR 14.7)
TE_all = method.VariableTE*1e-3; % converted to s
FA_all = method.VariableFlip+0.01; % avoid 0 flip angles
SRmax = SRmax.*0.8;

%% setting up trajectory
Nshots = 48; %# of spirals
[ktraj,dcf,t,ind,out,G]=design_spiral_pulseq(fov*1000,Nx,Nshots,10,...
    '',Gmax,SRmax,'1H',true,true,true);%fov is in mm
G = G.*gamma; %converting T/m to Hz/m

%% Rotating G, the G from design_spiral_pulseq function has rewinders but no rotation
% ktraj has rotations but no rewinders
G = complex(G(1,:),G(2,:));
G = [0,G]; %add a zero because Siemens does not like odd gradient length
Gn = zeros(size(G,2),Nshots);
for ks = 1:Nshots
    Gn(:,ks) = (G.*exp(1i*(ks-1)/Nshots*2*pi))';
end

%% Break ktraj into Nshots
kshot = zeros(length(ktraj)/Nshots,Nshots);
OneSpiralPoints = length(ktraj)/Nshots; %# of points in one spiral
for n1 = 1:Nshots
    kshot(:,n1) = ktraj((n1-1)*OneSpiralPoints+1:n1*OneSpiralPoints).*1e3;% per m
    plot (kshot(:,n1));hold on;grid on;xlabel('kx (m-1)'); ylabel('ky(m-1)');
end

%% Get gradient based on k data
ktrajs = zeros(size(ktraj,1), size(ktraj,2), 2);
ktrajs(:,:,1) = real(ktraj);
ktrajs(:,:,2) = imag(ktraj);

%% Construct flip angles and construct RF angles (RF shapes cant be more than 128 (Siemens limit))
adc = mr.makeAdc(length(Gn),'Dwell',system.gradRasterTime);
for itr = 1:total_acq
    [rf, gz] = mr.makeSincPulse(deg2rad(FA_all(itr)),system,'Duration', 1e-3,...
        'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4);
    rf_all(itr) = rf;
    gz_all(itr) = gz;
end

%% Need to determine multi-slice frequency offsets
Nslices=1;
deltaz=Nslices*sliceThickness;
z=-(deltaz/2):sliceThickness:(deltaz/2);
gzReph = mr.makeTrapezoid('z',system,'Area',-gz.area/2,'Duration',1e-3);
gx = mr.makeArbitraryGrad('x', squeeze(real(Gn(:,1))),system);
% [rf2] = mr.makeBlockPulse(pi,system,'Duration',500e-6);
[rf2, gz2] = mr.makeSincPulse(pi,system,'Duration', 2e-3,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4);
rf2.deadTime = 100e-6;rf2.ringdownTime=30e-6;%changing properties only for the rf block
gzSpoil = mr.makeTrapezoid('z',system,'Area',gz.area*2,'Duration',3*8e-4);

%% Define sequence block
numfullim = ceil(total_acq/Nshots);
TI = method.InversionTime*1e-3; % unit in s
% delay18 = TI - mr.calcDuration(gzSpoil)-mr.calcDuration(rf2)/2;
delay18 = TI - mr.calcDuration(rf2)/2; %changed due to TI calibration to get TI = 18 ms
delaySpoil = mr.makeDelay(delay18);
% seq.addBlock(gzSpoil);
seq.addBlock(rf2, gz2);
% seq.addBlock(gzSpoil);
seq.addBlock(delaySpoil);
time = 0;
for ds = 1:numfullim
    for ns=1:Nshots
        n = (ds-1)*Nshots+ns;
        if n <= total_acq %GE limitation, only 1000 acq point
            delayTE = mr.calcDuration(gzReph) + (mr.calcDuration(rf_all(n))/2);
            %% Calculate timing
            delayTR= TR_all(n) - TE_all(n) - mr.calcDuration(gx)/2;
            delayTR = delayTR - mod(delayTR, 1e-5);
            
            delay1 = mr.makeDelay(delayTE);
            delay2 = mr.makeDelay(delayTR);
            
            seq.addBlock(rf_all(n),gz_all(n));
            seq.addBlock(gzReph);
            gx = mr.makeArbitraryGrad('x', squeeze(real(Gn(:,ns))),system);
            gy = mr.makeArbitraryGrad('y', squeeze(imag(Gn(:,ns))),system);
            seq.addBlock(delay1);
            seq.addBlock(gx,gy,adc);
            seq.addBlock(delay2);
            time = time + TR_all(n);
        end
    end
end
disp(time/60);
seq.plot('TimeRange',[0 10*0.028]);
fname = ['Spiral_2D_vari_1000_single_slice_Orig_Par', num2str(Nshots),'_',num2str(TE_all(1))];
seq.write([fname,'.seq']);
end

