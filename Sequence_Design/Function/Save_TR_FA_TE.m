function Save_TR_FA_TE(TR_all, FA_all, TE_all)
    dim1 = size(TR_all);
    
    method.Method = 'vFISP';
    method.ExcPulse = {2};
    method.InvPulse = {10,'adiabatic'};
    method.FISP_inversion_enable = 'Yes';
    method.VariableFlip = FA_all;
    method.VariableTR = TR_all*1e3;
    method.TE = 1.9080;
    method.InversionTime = 18;
    method.VariablePhase = zeros(dim1(1),1);
    method.VariableTE = TE_all*1e3;
    
    save ('Save_TR_FA_TE.mat','method')

end