function [y,Fs,Format]=aiffread(aiffile,start,stop)
%AIFFREAD  Load AIFF  format sound files.
%   Note: The returned samples are in the range [-32768:32767]
%         and NOT scaled in any fashion 
%   [y]=AIFFREAD(aiffile,start,stop) loads a AIFF format file specified 
%       by "aiffile" from  frame no. start (Def=1) up to frame no. stop
%       (Def=length(file)) and returns the sampled data in variable "y". 
%
%   [y,Fs]=AIFFREAD(aiffile,start,stop) loads a AIFF format file specified 
%       by"aiffile" from frame no. start (Def=1) up to frame no. stop
%       (Def=length(file)) and  returns the sampled data in variable "y" 
%       and the sample rate in variable "Fs".
%   
%   [y,Fs,Format]=AIFFREAD(aiffile,start,stop) loads a AIFF format file 
%        specified by "aiffile" from  frame no. start (Def=1) up to  
%        frame no. stop (Def=length(file)) and  returns the sampled 
%        data in variable "y", the sample rate in variable "Fs", and the 
%        AIFF file format information in variable "Format". 
%        The format information is returned as a 2 element
%        vector with the following order:
%
%               Format(1)       Number of channels
%               Format(2)       Bits per sample
%
% JF: I found this on the web somewhere but it did not work for me.


if nargin<1 || nargin> 3
  error(sprintf('AIFFREAD takes up to three arguments, which are\n\tthe name of the AIFF file (required)\n\tthe first frame to read\n\tthe last frame to read'));
end

notend = 1;
comread = 0;
ssndread = 0;
if nargin==1
	start=1;
	stop=-1;
end

if nargin ==2 
	stop=-1;
end

if( stop> 0 && stop <= start)
	error('aiffread: No samles selected\n');
	return
end 

fid=fopen(aiffile,'rb','b');
if fid == -1
	error(['Can''t open AIFF file >',aiffile,'< for input!']);
end

% read aiff chunk
header=fscanf(fid,'%4s',1);
if  ~strcmp(header, 'FORM')
	error('No AIFF File!')
end

fsize = fread(fid, 1, 'long')
        
if 1
	header=fscanf(fid,'%4s',1);
	if ~strcmp(header, 'AIFF')
		disp(['FORM =', header])
		error('No AIFF File!')
	end
end

return
'here'


% read format sub-chunk
while notend == 1
	header=fscanf(fid,'%4c',1);
	len=fread(fid,1,'long');
	if strcmp(header,'COMM')
		if comread == 1
       			error('Error reading AIFF chunks: 2 COMM-Chunks detected!')
		end
		comread = 1;
       		Format(1)=fread(fid,1,'short');           % Channel
		disp(sprintf('Channels   : %d',Format(1)))
       		if Format(1) > 4
               		error('File contains more then 4 channels!')
       		end
		frames = fread(fid,1,'ulong');            % Frames
		Format(2)=fread(fid,1,'short');            % bits per sample
       		disp(sprintf('Bits/Sample: %d', Format(2)))
		expon=fread(fid,1,'short');
		tt   = fread (fid,8,'unsigned char');

		himant=((((tt(1)*256+tt(2))*256)+tt(3))*256+tt(4));
		himant=(((((himant *256 +tt(5))*256+tt(6))*256)+tt(7))*256+tt(8))/2^63;
		if( expon == 0 && himant == 0)
	      		Fs = 0;
		else
	      		Fs=(himant) * 2^ (abs(expon) -16383);
		end
       		disp(sprintf('Rate       : %.0f', Fs))

	elseif strcmp(header, 'SSND')
		if ssndread == 1
			error('Error reading AIFF chunks: 2 SSND-Chunks detected!')
		end
		ssndread = 1;
	  
		offset=fread(fid,1,'ulong');
		dummy=fread(fid,1,'ulong');
		if dummy ~= 0 
			disp('AIFF-SSND Chunk specify non zero blocksize ??')
		end
		fseek(fid,offset,0);


		if     Format(2) == 8
			formstr = 'schar';
		elseif Format(2) == 16
				formstr = 'short';
		elseif Format(2) == 32
			formstr = 'slong';
		end

		if start>1
			yd=fread(fid,Format(1)*(start-1),formstr);
		end	      
		if stop >start && stop < frames
			frames =stop;
		end	      
		yd=fread(fid,Format(1)*(frames-start+1),formstr);
		notend = 0;
	else
		disp(sprintf(' Skip over %s',header));
		fseek(fid,len,0);
	end
end
fclose(fid);

if	Format(1) == 1
	y = yd;
elseif	Format(1) == 2
	ss=size(yd);
	y = [yd(1:2:ss) yd(2:2:ss)];
elseif	Format(1) == 3
	ss=size(yd);
	y = [yd(1:3:ss) yd(2:3:ss) yd(3:3:ss)];
elseif	Format(1) == 4
	ss=size(yd);
	y = [yd(1:4:ss) yd(2:4:ss) yd(3:4:ss) yd(4:4:ss)];
end
