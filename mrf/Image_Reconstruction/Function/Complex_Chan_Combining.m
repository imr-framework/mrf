function coil_filter = Complex_Chan_Combining(image1)
% This script correct phase differences in each coil using method in [1].
% IUPUT
%     image1  reconstructed image 
%
% OUTPUT
% coil_filter  filter used to correct phase from each coil
% 
% [1] Roemer P, Edelstein W, Hayes C, Souza S, Mueller O. The NMR phased
% array. Magn Reson Med. 1990;16(2):192-225. doi:10.1002/mrm.1910160203
image1 = permute(image1,[1 2 4 3]);
[Nx, Ny, Nt, C] = size(image1);
filtsize = 5;
Rs = zeros(Nx, Ny, C, C);

% Correlation matrix
for n1 = 1:C
    for n2 = 1:C
        for n3 = 1:Nt
            Rs(:,:,n1,n2) = Rs(:,:,n1,n2) + filter2(ones(filtsize),...
                image1(:,:,n3,n1).*conj(image1(:,:,n3,n2)),'same');
        end
    end
end

% Get filter
Rs = permute(Rs,[3,4,1,2]);
Rs = Rs(:,:,:);
coil_filter = zeros(size(Rs,3),C);
parfor n4 = 1:size(Rs,3)
    [U,~] = svd(squeeze(Rs(:,:,n4)));
    coil_filter(n4,:) = U(:,1);
end

coil_filter = reshape(coil_filter, Nx, Ny, C);

end