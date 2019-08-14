%%
ktraj = spiral_data(:,1:48,:);
ktraj = reshape(ktraj,[1092*48,2]).*1e3;
plot(squeeze(ktraj));

%% 
om = squeeze(ktraj.*2.*pi);
st = nufft_init(om, Nd, Jd, Kd, Nd/2, 'minmax:kb');
kspace_data_1pt = squeeze(kspace_data_ch(:,1:48,:));
for ch =1:size(kspace_data_1pt,3)
    kdata_temp = reshape(squeeze(kspace_data_1pt(:,:,ch)),[1 1092*48]).';
    image_data_ch(:,:,ch) = nufft_adj(kdata_temp,st);   
end

im_sos = sum(abs(image_data_ch),3);
figure; imagesc(im_sos);