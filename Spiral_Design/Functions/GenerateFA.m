function FA_shots_final = GenerateFA

Scales = [4,8,11,6];
NumOfFA = 250;
FA_shots = zeros(NumOfFA, length(Scales));
for n1 = 1:length(Scales)
    FA_read = Scales(n1).*hanning(NumOfFA);%This is to limit # of flip angles so it does not go over 128 limit for Simmens
    for n2 = 1: NumOfFA/10
        temp = 10*(n2-1)+1;
        FA_shots(temp:temp+9,n1) = FA_read(temp);
    end
end
FA_shots_final = cat(1,FA_shots(:,1),FA_shots(:,2),FA_shots(:,3),FA_shots(:,4));
end