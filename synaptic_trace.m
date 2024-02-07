function F=synaptic_trace(Tau,dt,N_e,N_i)

    F=struct();
    tau_sE=Tau.tausyn_e; % exc synaptic time (fall time)
    tau_sI=Tau.tausyn_i; % inh synaptic time (fall time)
    fexp=[repmat(exp(-dt/tau_sE),N_e,1); repmat(exp(-dt/tau_sI),N_i,1)]; % Multiplicative step (fp)
    fSpike=[repmat((1/tau_sE),N_e,1); repmat((1/tau_sI),N_i,1)]; % add to fp with a spike
    f=zeros(N_e+N_i,1);
    F.fexp=fexp;
    F.fSpike=fSpike;
    F.f=f;

end