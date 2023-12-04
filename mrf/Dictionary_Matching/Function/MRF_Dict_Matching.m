function Par_map = MRF_Dict_Matching(MRF_dict, data)

% This script uses vector dot product to find the most similar entry inside
% a dictionary and retrieve the parameters from that dictionary.
% INPUT
%  MRF_dict  Dictionary file
%      data  Reconstructed data file
%
% OUTPUT
%   Par_map  An output structure that contains T1 and T2 maps
%
% 7/2019 Enlin Qian
% # Copyright of the Board of Trustees of Columbia University in the City of New York

%% Flatten image data
data = abs(data);
image_size = size(data);
num_pixel = image_size(1)*image_size(2);
num_timepts = image_size(3);
data_flat = reshape(data,[num_pixel,num_timepts]);
mask = ones(image_size(1:end-1)) > 0;
data_flat = data_flat(mask>0,:);
data_size = size(data_flat);

%% Reading in dictionary
dict = abs(MRF_dict.dict_SW_norm);
dict_size = size(dict);

%% Dictionary Matching Data Preparation
dict_max = max(dict,[],2);
dict_temp = abs(dict./dict_max);
data_max = max(data_flat,[],2);
data_temp = abs(data_flat./data_max);
ip = zeros(dict_size(1),1);
match_point = zeros(data_size(1),1);
dict_match = zeros(data_size(1),1);
dict_progress = ceil((1:1:data_size(1))./data_size(1)*100);

%% Dictionary Matching
fprintf('Dictionary match starts \n');
for n1 = 1:data_size(1)
    if mod(n1, floor((data_size(1))/10))==0
        fprintf('%d percent has been completed. \n' ,dict_progress(n1));
    end
    ip = sum(abs(dict_temp.^2 - data_temp(n1, :).^2), 2);
    [match_point(n1),dict_match(n1)] = min(abs(ip),[],1);
end

%% Generate Output
Par_map.T1_map = reshape(MRF_dict.lut(dict_match,1), [image_size(1), image_size(2)]);
Par_map.T2_map = reshape(MRF_dict.lut(dict_match,2), [image_size(1), image_size(2)]);
figure; imagesc(abs(squeeze(Par_map.T1_map))); axis equal tight; colormap hot; title('T1 map');
figure; imagesc(abs(squeeze(Par_map.T2_map))); axis equal tight; colormap hot; title('T2 map');
end
