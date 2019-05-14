phis = zeros(1,1001);
alphas = [180, FA_all'];
TRs = [18, TR_all'.*1000];
% dt = 1e-3;
rlxWM = [848.5, 62.74]; %GM: T1 950 T2 100 WM: T1 600 T2 80 CSF: T1 4500 T2 2200
%https://onlinelibrary.wiley.com/doi/epdf/10.1002/%28SICI%291522-2586%28199904%299%3A4%3C531%3A%3AAID-JMRI4%3E3.0.CO%3B2-L
rlxGM = [1687.02, 106.8];
rlxCSF = [4000, 2000];

[om_store_WM,echoes,seq] = EPGsim_MRF(phis,alphas,TRs,rlxWM);
Echo_Final_WM = Var_TE(om_store_WM, TE_all, rlxWM(2));
[om_store_GM,echoes,seq] = EPGsim_MRF(phis,alphas,TRs,rlxGM);
Echo_Final_GM = Var_TE(om_store_GM, TE_all, rlxGM(2));
[om_store_CSF,echoes,seq] = EPGsim_MRF(phis,alphas,TRs,rlxCSF);
Echo_Final_CSF = Var_TE(om_store_CSF, TE_all, rlxCSF(2));

plot(Echo_Final_WM)
hold on
plot(Echo_Final_GM)
plot(Echo_Final_CSF)
hold off
legend ('WM','GM','CSF')