function Echo_Final = Var_TE(om_store, TEs, T2)
F_plus = cellfun(@(v)v(1), om_store);
F_plus = F_plus(1:3:end);
Echo_Final = F_plus.*exp(-TEs./T2);
end