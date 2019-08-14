function Echo_Final = Var_TE(om_store, TE_all, T2)
F_plus = cellfun(@(v)v(1),om_store);
F_plus = F_plus(1:3:end);
F_plus = abs(F_plus);
F_plus = F_plus(2:1001);
Echo_Final = F_plus.*exp(-TE_all'/T2);
end