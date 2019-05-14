function [acquired_raw_data_con2, MaxValue, MaxLoc] = Optimized_Shift(acquired_raw_data_con)
% This function takes in an FID signal raw data file and finds the peak.
% Then remove the points before the peak and add same number of zeros in
% the end. It is because the density of the spiral is larger in the middle
% and smaller when it goes out, thus DCF is used to balanced the density
% differences so that Fourier transform can be done correctly. However, FID
% signal takes time to reach its maximum and DCF is 0 at the beginning.
% When FID reaches maximum, DCF is already 1 and does not have the ability
% to balance the density differences anymore. So we finds the peak and
% match them with the DCF manully here.

%% Check dimension of input
dim1 = size(acquired_raw_data_con);
acquired_raw_data_con2 = zeros(dim1);
MaxValue = zeros(dim1(2:3));
MaxLoc = zeros(dim1(2:3));
for n1 = 1:dim1(2)
    for n2 = 1:dim1(3)
        [MaxValueTemp, MaxLocTemp] = max(acquired_raw_data_con(:,n1,n2));
        ZeroTail = zeros(MaxLocTemp-1,1,1);
        DataTemp = acquired_raw_data_con(MaxLocTemp:end,n1,n2);
        acquired_raw_data_con2(:,n1,n2) = cat(1,DataTemp,ZeroTail);
        %         acquired_raw_data_con2(:,n1,n2) = [DataTemp,ZeroTail];
        
        MaxValue(n1,n2) = MaxValueTemp;
        MaxLoc(n1,n2) = MaxLocTemp;
    end
end