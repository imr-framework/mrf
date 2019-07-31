 function sound(x, S)
%function sound(x, S)
%|
%| sound command for octave on a mac
%| requires "lame" and "mpg123" packages to be installed

if nargin < 1, help(mfilename), error(mfilename), end

if nargin < 2
	S = 8192; % default sampling rate
end

tmp_wav = '/tmp/sound1.wav';
tmp_mp3 = '/tmp/sound1.mp3';

com1 = ['lame --quiet --preset standard ' tmp_wav ' ' tmp_mp3];
com2 = ['mpg123 --quiet ' tmp_mp3];


wavwrite(x, S, 8, tmp_wav)

if ~exist(tmp_wav, 'file')
	ir_fail('cannot find file "%s"', tmp_wav)
end

os_run(com1)

if ~exist(tmp_mp3, 'file')
	ir_fail('cannot find file "%s"', tmp_mp3)
end

os_run(com2)
