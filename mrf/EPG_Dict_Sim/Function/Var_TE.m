function Echo_Final = Var_TE(echoes, TE_all, T2)
Echo_Final = echoes(:,2).*exp(-TE_all'/T2);
end