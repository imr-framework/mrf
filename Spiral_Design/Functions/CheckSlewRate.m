function [LinearFill] = CheckSlewRate(Array1, Array2, Array3, SlewRate, gsamp)
% This function takes in three arrays and checks if the slope is bigger than the SlewRate. If
% no, concatenate these three arrays. If yes, 

A = [Array1, Array2, Array3];
% DiffA = diff(real(A))./gradRasterTime;
smax_act = max(max(abs(diff(A,1))))/gsamp;
n = 1;

while smax_act>SlewRate
    smax_act = max(max(abs(diff(A,1))))/gsamp;
    LinearFill = [length(Array1)+1,(Array1(end)+Array2(1))/2];
%     LinearFill = [Array1,point(2),Array2]
else
    LinearFill = [];
end

end