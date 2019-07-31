function rng(seed)
% 2012-09-22 octave 3.6.3 does not provide rng()

rand('seed', seed)
randn('seed', seed)
