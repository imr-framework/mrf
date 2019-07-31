% newarr = hflip(arr [,cp]);
%
% Flips the array about the horizontal axis.
%
% If cp is given, the array is flipped about the point 
% cp, where 2*cp must be an integer.
%
% $Id: hflip.m,v 1.1 2010/03/15 02:18:35 jfnielse Exp $


% =============== CVS Log Messages ==========================
%	This file is maintained in CVS version control.
%
%	$Log: hflip.m,v $
%	Revision 1.1  2010/03/15 02:18:35  jfnielse
%	
%	: Committing in .
%	:
%	: Added Files:
%	: 	README hflip.m make2dft.m makeEPI.m makecrusher.m makegz.m
%	: 	readjfnwav.m trapwave.m writejfnwav.m
%	
%	Revision 1.1  2008/09/21 14:37:51  jfnielse
%	
%	: Modified Files:
%	: 	makeEPI.m
%	: Added Files:
%	: 	hflip.m trapwave.m
%	
%	Revision 1.1  2004/04/05 16:20:35  knayak
%	initial Matlab Directory
%	
%	Revision 1.1  2002/03/28 01:35:30  bah
%	Added to CVS
%	
%
% ===========================================================


function newarr = hflip(arr, cp)

s = size(arr);
s2 = s(2);

if (nargin < 2)
	cp = (s2+1)/2;
end;

cp = cp+s2;
temp = [arr, arr, arr];

newarr = 0*arr;

for k = 1:s2
  newarr(:,k) = temp(:,2*cp-k-s2);
end;






