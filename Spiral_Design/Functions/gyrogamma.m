function g = gyrogamma(nucleus)
%GYROGAMMA  Gyromagnetic ratio of various spins (in radians)
% g = gyrogamma(nucleus)
%        g [rad/(s T)]
%  nucleus = 'h1';'h2';'he3';'c13';'n15';'f19';'na23';'p31';'xe129';...
%            '13c1h': ratio gyrogamma('13C')/gyrogamma('1h')
%            'h1' is default
%
% Literature:
%        Levitt; "Spin Dynamics"; Wiley
%        DeGraaf; "In Vivo NMR Spectroscopy"; Wiley
%        Ratio: jbnmr6_135_1995_wishart,pac70_117_1998_markley
%        electron: wikipedia
%
% 11/2009 Rolf Schulte
if ~exist('nucleus','var'), nucleus = 'h1'; end
switch lower(nucleus)
    case {'h1','1h',1},      g = 267.522128d6;
    case {'h2','2h',2},      g =  41.066d6; fprintf('Spin 1\n');
    case {'h3','3h'},        g = 285.349d6;
    case {'he3','3he',3},    g =-203.789d6;
    case {'li7','7li',7},    g = 103.962d6; % source wikipedia
    case {'b10','10b',10},   g =  28.747d6; fprintf('Spin 3\n');
    case {'b11','11b',11},   g =  85.847d6; fprintf('Spin 3/2\n');
    case {'c13','13c',13},   g =  67.283d6;
    case {'n14','14n',14},   g =  19.338d6; fprintf('Spin 1\n');
    case {'n15','15n',15},   g = -27.126d6;
    case {'o17','17o',17},   g = -36.281d6; fprintf('Spin 5/2\n');
    case {'f19','19f',19},   g = 251.815d6;
    case {'ne21','21ne',21}; g = -2.113d7;  fprintf('Spin 3/2\n');
    case {'na23','23na',23}, g =  70.808d6; fprintf('Spin 3/2\n');
    case {'al27','27al',27}, g =  69.763d6; fprintf('Spin 5/2\n');
    case {'si29','29si',29}, g = -53.190d6;
    case {'p31','31p',31},   g = 108.394d6;
    case {'cl35','35cl',35}, g =  10.610d6; fprintf('Spin 3/2\n');
    case {'cl37','37cl',37}, g =   8.832d6; fprintf('Spin 3/2\n');
    case {'k39','39k',39},   g =  12.5d6;   fprintf('Spin 3/2\n');
    case {'mn55','55mn',55}, g = 6.6317e+07; fprintf('Spin 5/2\n');
    case {'cu63','63cu',63}, g =  71.118d6; fprintf('Spin 3/2\n');
    case {'cu65','65cu',65}, g =  76.044d6; fprintf('Spin 3/2\n');
    case {'kr83','83kr',83}, g = -1.033d7;  fprintf('Spin 9/2\n');
    case {'ag107','107ag',107}, g = -10.889d6;
    case {'ag109','109ag',109}, g = -12.518d6;
    case {'xe129','129xe',129}, g = -74.521d6;
    case {'xe131','131xe',131}, g = 2.209d7;   fprintf('Spin 3/2\n');
    case {'cs133','133cs',133}, g = 35.0878d6; fprintf('Spin 7/2\n');
    case {'pb207','207pb',207}, g =  55.805d6;
    case {'13c1h'},         g = 0.251449530; fprintf('Ratio 13C/1H\n');
    case {'e'},             g = 1.760859770d11;
    otherwise, error('Input undefined');
end
