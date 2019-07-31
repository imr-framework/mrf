%
%	function [waveform,t] = trapwave(area,T,gmax,smax)
%
%	Function returns the minimum-duration waveform with
%	the given area.
%
%	INPUT:
%		area = area in G*s/cm.
%		T = gradient sampling period (s)
%		gmax = maximum gradient amplitude (G/cm)
%		smax = maximum gradient slew rate (G/cm/s)
%
%	OUTPUT:
%		waveform = gradient waveform (G/cm).
%
%	Brian Hargreaves, July 2002.
%
% $Id: trapwave.m,v 1.1 2010/03/15 02:18:35 jfnielse Exp $

% ---------------------------------------------------------------
%  CVS Log Messages
%
%  Log: trapwave.m,v 
%  Revision 1.1  2004/04/05 16:20:36  knayak
%  initial Matlab Directory
%
%  Revision 1.3  2003/10/20 22:28:16  brian
%  minor edits
%
%  Revision 1.2  2002/09/04 00:34:55  bah
%  Added warning note if passed area appears to be in wrong units.
%
%  Revision 1.1  2002/07/26 00:03:43  bah
%  New moment rewinding functions added
%
% ---------------------------------------------------------------


function [waveform,t] = trapwave(area,T,gmax,smax,gstart)


% ============	Defaults for gmax and smax.   =============

if nargin < 2
	T = 4e-6;	% s.
end
if nargin < 3
	gmax = 3.9;	% G/cm.
end
if nargin < 4
	smax = 14500;	% G/cm/s.
end
if nargin < 5
	gstart = 0;
end

dg = smax*T;

% Bring gradient to 0 using max slew rate and recalculate area
if gstart == 0
	sramp = [];
	sarea = 0;
else
	nr = floor(gstart/dg);
	sramp = (nr:-1:0)*dg*sign(gstart);
	sarea = sum(sramp)*T;
end

sgn = sign(area - sarea);
garea = abs(area - sarea);	% Same for positive and negative areas.

if garea > .1
	disp('Note:  trapwave expects area in G*s/cm, not cm^(-1)');
end

% ============  Calculate Rise Time to full scale ============

tr = gmax/smax;		% Time for ramp to full scale.
Acrit = smax*tr^2;	% Area of maximum width triangle.

if (garea < Acrit) 	% =========== Case 1, time < tr ==============
	rtime = sqrt(garea/smax);
	nr = ceil(rtime/T);
	ramp = (0:nr)*dg;
	waveform = [ramp flipdim(ramp, 2)];

else			 % =========== Case 1, time > tr ==============

	nr = ceil(tr/T);
	ramp = (0:nr)*dg;
	np = ceil((garea-Acrit)/gmax/T);
	plat = gmax*ones(1,np);
	waveform = [ramp plat flipdim(ramp, 2)];
end;

% =========== Shift waveform down until area is correct. =========

actarea = sum(waveform)*T;
%tt = sprintf('Actual |area| = %f, desired |area| = %f ',actarea, garea);
%disp(tt);

while (actarea > garea+4/32767*T)
	waveform = waveform - (actarea-garea)/(length(waveform)*T);
	waveform = max(waveform,0);
	while (waveform(1)==0)
		waveform=waveform(2:end);
	end;
	while (waveform(end)==0)
		waveform=waveform(1:end-1);
	end;
	actarea = sum(waveform)*T;
	%tt = sprintf('Actual |area| = %f, desired |area| = %f ',actarea, garea);
	%disp(tt);
end;


waveform = waveform * sgn;
waveform = [sramp waveform];
t = (1:length(waveform))*T;



% actarea = sum(waveform)*T;
%tt = sprintf('Actual area = %f, desired area = %f ',actarea, area);
%disp(tt);


% Test for slew violation.

stest = [waveform 0] - [0 waveform];
[y,I] = max(abs(stest)/T);
if (y>1.01*smax)
        tt=sprintf(' --- Warning (trapwave):  Slew Violation, sample %d ------',I);
        disp(tt);
end;

