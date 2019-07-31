function y = flip(varargin)
% function y = flip(varargin)
% Matlab 2015 or so says "flipdim" is obsolete and replace with "flip"
y = flipdim(varargin{:});
