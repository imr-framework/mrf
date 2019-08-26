function [CircleMap, PointsIndex] = Circle_Analysis(Center, Radius, ImXSize, ImYSize)
% Input center, radius, and map size and return index of pixels inside the
% circle.
% INPUT
%     Center  Center of circle [x, y]
%     Radius  Radius of circle (in pixels)
%    ImXSize  Image size along x direction
%    ImYSize  Image size along y direction
%
% OUTPUT
%  CircleMap  Binary map for each circle
% PointsIndex  Pixel index for all pixels inside input circle
%
% 7/2019 Enlin Qian
% # Copyright of the Board of Trustees of Columbia University in the City of New York
% (x,y) = (colum, row)
[X,Y] = meshgrid(1:ImXSize,1:ImYSize);
DistanceMap = sqrt((X-Center(1)).^2+(Y-Center(2)).^2);
CircleMap = DistanceMap<=Radius;
PointsIndex = find(CircleMap);
end