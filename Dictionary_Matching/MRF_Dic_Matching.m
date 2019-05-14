% function out = MRF_Dic_Matching(dict, data)
%% This function is a re-implementation based on P. Gomez's code

%% Parameters for par because generation of par is not avaiable so need to assign values to par
par.fp.percentageOfUse = 0.6;
par.f.verbose = 1;
par.f.verbose_iter = 0;
par.f.Xout = 1;
par.f.dmout = 1;
par.f.mtout = 1;
par.f.qout = 1;
par.ind.Q = 6;

%% Reshape data
iz = 1;
% % data.X = image_data_final;
% dict.D = dict.D(:,1:961);
%%
[ix,iy,T] = size(data.X);
par.ind.datadims = [ix,iy,T];
Q = par.ind.Q;
datadims = par.ind.datadims;
N = ix*iy*iz;
x = reshape(data.X,[N,T]);
x = single(x(data.mask>0,:));
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
memPercentage = mem.MemAvailableAllArrays * par.fp.percentageOfUse;
par.fp.blockSize = memPercentage/dataBytes;

%% Dictionary Matching 
blockSize = min(max(floor(par.fp.blockSize/size(D,1)),1),size(x,1));
iter = ceil(size(x,1)/blockSize);
mt = zeros(size(x,1),1);
dm = zeros(size(x,1),1);
pd = zeros(size(x,1),1);
X = zeros(size(x));
if par.f.verbose; fprintf('Matching data \n'); end
maxDict =  max(dict.D,[],2);
for i = 1:iter
    if(mod(i,1000) ==0)
        disp(i)
    end
    if par.f.verbose_iter; fprintf(['Matching block ', num2str(i) '/' num2str(iter), '\n']); end
    if i<iter
        cind = (i-1)*blockSize+1:i*blockSize;    
    else
        cind = (i-1)*blockSize+1:size(x,1); 
    end

%%
xtemp = x(cind,:)./max(abs(x(cind,:)));
ip = zeros(1,size(D,1));
if(sum(abs(x(cind,:))) > 16 )

        for m = 1:size(D,1)
%           ip=dict.D*ctranspose(x(cind,:)); 
            dtemp = abs(dict.D(m,:)./maxDict(m));%This line can be done apriori
            
            ip(m) = sum(abs(dtemp.^2 - xtemp.^2));
%             ip(m) = max(abs(ipt));
%         for ic = 1:length(cind)
%             pd(cind(ic)) = ip(dm(cind(ic),1),ic);
%             X(cind(ic),:) = pd(cind(ic)).*D(dm(cind(ic),1),:); %re-scale X as in Davies et al.
%             pd(cind(ic)) = pd(cind(ic))./normD(dm(cind(ic),1)); %get PD after signal scaling. 
%         end       
        end
           [mt(cind),dm(cind)] = min(abs(ip.'),[],1);
          
           if(mod(i, 100)==0)
               disp(i);
                disp(dm(cind));
           end
end
end
clear match;
clear dictentry;

%% Generate output
% if par.f.Xout
% %     Xout = zeros(N,T,'single');
%     Xout = single(X);
%     out.Xfit = reshape(Xout,[datadims(1:end-1),T]);
%     out.X = data.X;
% end
% 
% if par.f.mtout
% %     mfit = single(zeros(N,1));
%     mfit = single(mt);
%     out.mt =  reshape(mfit,[datadims(1:end-1),1]);
% end
% 
% 
% if par.f.dmout
% %     dmfit = single(zeros(N,1));
%     dmfit = single(dm);
%     out.dm =  reshape(dmfit,[datadims(1:end-1),1]);
% end
dm(dm==0) =1;
dict.lut(1,:) =0;
if par.f.qout
    qmap = zeros(N,Q,'single');
    qmap(data.mask>0,:) = dict.lut(dm,:);
    qmap(isnan(qmap)) = 0; %remove possible NaN from infeasible search window
    out.qmap = reshape(qmap,[datadims(1:end-1),Q]);
    out.mask = data.mask;
end
figure; imagesc(abs(squeeze(out.qmap(:,:,1)))); axis equal tight; colormap hot;
figure; imagesc(abs(squeeze(out.qmap(:,:,2)))); axis equal tight; colormap hot;
% 
% end