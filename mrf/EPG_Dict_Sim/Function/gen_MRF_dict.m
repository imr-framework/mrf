function MRF_dict = gen_MRF_dict(TR, TE, FA, T1, T2, TI, kshots)
% This script uses EPG to generate a dictionary based on the input of TR, 
% flip angle, TE, and T1/T2 range. The dictionary is processed using
% sliding window to match with the image reconstruction method. The
% dictionary is then normalized based on the norm.
% IUPUT
%       TR  Repetition time (s)
%       TE  Echo time (s)
%       FA  Flip angles (degree)
%       T1  T1 values for dictionary (s)
%       T2  T2 values for dictionary (s)
%       TI  Inversion time, if TI=0, no inversion pulse is added. If TI~=0,
%           a 180 pulse with input TI value is added. (s)
%   kshots  number of shots of spiral
%
% OUTPUT
% MRF_dict  A structure that contains a look up table for T1 and T2 values
%           used for each entry and a EPG simulated dictionary before
%           sliding window, after sliding window, and after being
%           normalized by the norm. 
% 
% 7/2019 Enlin Qian
% # Copyright of the Board of Trustees of Columbia University in the City of New York

%% Check if 180 pulse needs to be added
if TI ~= 0  % if TI is not 0, add 180 pulse to Flip angles
    alphas = [180, FA'];
    TRs = [TI, TR'.*1000];
    TEs = [0, TE'.*1000];
else
    alphas = FA';
    TRs = TR'.*1000;
    TEs = TE'.*1000;
end

%% Prepare variables for dictionary simulation
T1 = T1.*1000;
T2 = T2.*1000;
num_TRs = size(TRs, 2);
num_dict_entry = length(T1)*length(T2);
dict_size = [num_dict_entry, num_TRs];
phis = zeros(1, num_TRs);
Echo_Final = zeros(dict_size);
dict_lut = zeros(num_dict_entry, 2);
dict_SW = zeros(dict_size);
dict_SW_norm = zeros(dict_size);

%% Generate EPG dictionary
n3 = 1;
dict_progress = ceil(((1:1:num_dict_entry)./num_dict_entry)*100);
for n1 = 1:length(T2)
    for n2 = 1:length(T1)    
        if mod(n3, floor(num_dict_entry/10))==0
            fprintf('%d percent has been completed. \n' ,dict_progress(n3));
        end
        [om_store,echoes,seq] = EPGsim_MRF(phis,alphas,TRs,[T1(n2), T2(n1)]);
        Echo_Final(n3, :) = Var_TE(om_store, TEs, T2(n1));
        dict_lut(n3, 1) = T1(n2);
        dict_lut(n3, 2) = T2(n1);
        n3=n3+1;        
    end
end
MRF_dict.lut = dict_lut;
MRF_dict.dict = Echo_Final;
clear dict_lut dict_progress

%% Sliding window
kshot_win = floor(kshots.*0.5);
for n4 = 1: num_TRs
    SelectedRange = (n4-kshot_win-1):(n4+kshot_win-2);
    SelectedRange = SelectedRange(SelectedRange > 0 & SelectedRange <= num_TRs);
    dict_SW(:, n4) = sum(Echo_Final(:,SelectedRange),2);
end
MRF_dict.dict_SW = dict_SW;
clear Echo_Final

%% Normalized by the norm
dict_norm = zeros(num_dict_entry, 1);
for n5 = 1: num_dict_entry
    dict_norm(n5) = norm(dict_SW(n5,:));
    dict_SW_norm(n5, :) = dict_SW(n5, :) / dict_norm(n5);
end
dict_SW_norm = dict_SW_norm(:, 2:1001);
MRF_dict.dict_SW_norm = dict_SW_norm;
clear dict_SW_norm dict_SW

end