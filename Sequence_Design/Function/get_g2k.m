function  [kx,ky,kz] = get_g2k(Gxp,Gyp,Gzp)
% gamma = 42576000; %From Pulseq 

dt = diff(Gxp.t);
if(nargin < 3)
kx = cumsum(Gxp.waveform)*dt(1);
ky = cumsum(Gyp.waveform)*dt(1);
kz =0;
elseif(nargin <4)
kx = cumsum(Gxp.waveform)*dt(1);
ky = cumsum(Gyp.waveform)*dt(1);
kz = cumsum(Gzp.waveform)*dt(1);
end
% for l = 1:length(Gxp.waveform)
%         kx(l) = trapz(Gxp.waveform(1:l))*dt(1);
%         ky(l) = trapz(Gyp.waveform(1:l))*dt(1);
%         kz(l) = trapz(Gzp.waveform(1:l))*dt(1);
% end