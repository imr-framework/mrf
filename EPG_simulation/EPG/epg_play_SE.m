seq.rf = [0 90
     90 180];
seq.grad = [1 1];
seq.events = {'rf','grad','relax','rf','grad','relax'};
seq.time = [0, 10, 10, 10, 20, 20];
seq.T1 = 1000;
seq.T2 = 0;
seq.name = "SE";
[om_store,echoes] = EPG_custom(seq);
display_epg(om_store,seq,1)