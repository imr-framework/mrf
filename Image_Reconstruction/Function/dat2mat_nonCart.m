function [kspace_data]= dat2mat_nonCart()


%% Author - Sairam Geethanath
% Date 20130824
% Input - *.dat file
% Output - Kspace, Image space
% We may have to have a relook at this code after ISMRM 2014 due to its
% rigidity for N x N acquisitions


%% Modification history
% to make it generic - remove hardcoding
% 06/09/2014
%% Obtain the file
[Filename,Pathname] = uigetfile('*.dat','Pick the raw data file');


%% Read data using mapVBVD
% image_obj = mapVBVD(fullfile(Pathname,Filename));
twix_obj = mapVBVDVE(fullfile(Pathname,Filename));

% image_obj = twix_obj.image; %Body coil
if(~iscell(twix_obj))
% if ( size(twix_obj.image) ==1)
    image_obj = twix_obj.image;
else
    image_obj = twix_obj{2}.image;
end

sizeData = image_obj.sqzSize; %kx, nch, ky, slices, partitions
dimsData = image_obj.sqzDims;
dimsallData = image_obj.dataDims;

%% Concatenate slices and partitions together and has a N x N acquisition set
kspace_data = zeros(sizeData);

sl = 1;

partitions=1;

%% Determine k-space shift to the center
% temp = squeeze(image_obj(:,1,:,round(sizeData(4)/2),round(sizeData(5)/2)));
% [~,idx] = max(abs(temp(:)));
% [idx,idy] = ind2sub(size(temp),idx);

%%


tic;
% for nch=1:size(kspace_data,2) %number of channels

%        for narms=1:size(kspace_data,1) %number of arms
         kspace_data = squeeze(image_obj());
         
%         temp = squeeze(image_obj(:,nch,:,sl,partitions));
%         temp = circshift(temp, [round(size(temp,1)/2) - idx,round(size(temp,2)/2) - idy ]);
%         trunc_part =round(0.5 .*(size(temp,1) - size(temp,2)));
%         temp = squeeze(temp(trunc_part+1:end-trunc_part,:));
%         
%         kspace_data(:,:,nslices,nch) = temp;

      


% end
toc;

kspace_data = permute(kspace_data,[1  3   2   4]); %kx, ky, slices, partitions, nch
kspace_data = squeeze(kspace_data);











