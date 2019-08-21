function out = MRF_Dict_Matching(dict, data)
%% This function is a re-implementation based on P. Gomez's code
%This script uses vector dot product to find the most similar entry inside
%a dictionary and retrieve the parameters from that dictionary.
% Input
%  dict  Dictionary file
%  data  Reconstructed data file
% OUTPUT
%   out  An output structure that contains T1 and T2 maps 

%% Reshape data
iz = 1;
[ix,iy,T] = size(data);
datadims = [ix,iy,T];
N = ix*iy*iz;
x = reshape(data,[N,T]);
mask = ones(datadims(1:end-1)) > 0;
x = single(x(mask>0,:));

%% Reading in dictionary and checking fields in dictionary
D = dict.D;
if isfield(dict,'normD'); normD = dict.normD; end
if ~isfield(dict,'lut')
    error('dict.lut required for parameter estimation');
end
if ~exist('normD','var')
    warning('normD not found, normalizing dictionary \n');
    normD = zeros(1,size(D,1));
    for l = 1:size(D,1)
        normD(l)=norm(D(l,:));
        D(l,:)=D(l,:)/normD(l);
    end
    D(isnan(D))=0;
end

for n1 = 1:size(D,1)
    temp1 = D(n1,:);
    temp2 = x(n1,:);
    temp1 = temp1./max(temp1(:));
    temp2 = temp2./max(temp2(:));
    D(n1,:) = temp1;
    x(n1,:) = temp2;
end

%% Setting up parameters for dic matching
S = whos('x');
dataBytes = S.bytes;
mem = memory;
percentageOfUse = 0.6;
memPercentage = mem.MemAvailableAllArrays * percentageOfUse;
blockSize = memPercentage/dataBytes;

%% Dictionary Matching
blockSize = min(max(floor(blockSize/size(D,1)),1),size(x,1));
iter = ceil(size(x,1)/blockSize);
mt = zeros(size(x,1),1);
dm = zeros(size(x,1),1);
pd = zeros(size(x,1),1);
X = zeros(size(x));
fprintf('Matching data \n');
maxDict =  max(dict.D,[],2);
dtemp = abs(dict.D./maxDict);
for i = 1: iter
    if(mod(i,1000) ==0)
        disp(i)
    end
    if i<iter
        cind = (i-1)*blockSize+1:i*blockSize;
    else
        cind = (i-1)*blockSize+1:size(x,1);
    end
    xtemp = x(cind,:)./max(abs(x(cind,:)));
    ip = zeros(1,size(D,1));
    for m = 1:size(D,1)
        ip(m) = sum(abs(dtemp(m, :).^2 - xtemp.^2));
    end
    [mt(cind),dm(cind)] = min(abs(ip.'),[],1);
    if(mod(i, 100)==0)
        disp(i);
        disp(dm(cind));
    end
end
clear match;
clear dictentry;

%% Generate output
Q=6;
dm(dm==0) =1;
dict.lut(1,:) =0;

qmap = zeros(N,Q,'single');
qmap(mask>0,:) = dict.lut(dm,:);
qmap(isnan(qmap)) = 0; %remove possible NaN from infeasible search window
out.qmap = reshape(qmap,[datadims(1:end-1),Q]);
out.mask = mask;
figure; imagesc(abs(squeeze(out.qmap(:,:,1)))); axis equal tight; colormap hot;
figure; imagesc(abs(squeeze(out.qmap(:,:,2)))); axis equal tight; colormap hot;
end
