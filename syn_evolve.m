function F=syn_evolve(F,fired)


% update synaptic filter
    F.f=F.fexp.*F.f;
    if ~isempty(fired)      
    % update of synaptic filter with spikes
        F.f(fired)=F.f(fired)+F.fSpike(fired);
    end
end

