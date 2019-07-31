% ncat_segment
% the ncat phantom is not binary and is 140kev 

%file = 'ncat,256,slice,140.fld';
clear all;
file = 'ncat,1024,slice.fld';
x = 1000 * double(fld_read(file));
x = round(x);
im(x)
hist(x(x~=0),100);
un = unique(x(:))'

%	0	background
%	3	lungs   
%	9	body    
%	11	spine   
%	15	ribs    

map = zeros(1+max(un),1);
%map(un+1) = [0 0.2 1 1.9 2];
map(un+1) = [0:length(un)-1];
y = map(x+1);
im(y)
clear file map un x;

ribs = zeros(size(y));
ribs(y == 4) = 2; %follow Jeff's 
spine = zeros(size(y));
spine(y == 3) = 1.9;
body = zeros(size(y));
body(y == 2) = 1;
lungs = zeros(size(y));
lungs(y == 1) = 0.2;
xtrue(:,:,1) = body+lungs;
xtrue(:,:,2) = ribs+spine;

clear y
save 'xtrue1024.mat';

%if 1
%	y = single(y);
%	fld_write('ncat,256,slice,140,ct,dens.fld', y)
%	y = uint8(y);
%	fld_write('ncat,1024,slice500,0to4.fld', y, 'datatype', 'byte')
%end
