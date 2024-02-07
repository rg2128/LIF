% SIM of one trials given parameters
%
% Luca Mazzucato March 2014

% SET OPTIONS
% ParamsRun = structure containing parameters for simulation


function [all_firings, PlotData, J]=fun_LIF_SIM(ParamsRun,Inh,i,CondIter)

Odor = i;
Theta=ParamsRun.Theta; Sim=ParamsRun.Sim; Stimulus=ParamsRun.Stimulus;
Ext=ParamsRun.Ext; J=ParamsRun.J; N_e=ParamsRun.N_e; N_i=ParamsRun.N_i; p=ParamsRun.p; He=ParamsRun.He;
Hi=ParamsRun.Hi; tau_e=ParamsRun.tau_e; tau_i=ParamsRun.tau_i; tausyn_e=ParamsRun.tausyn_e;
tausyn_i=ParamsRun.tausyn_i; tau_arp=ParamsRun.tau_arp;
%
all_firings=[];
dt=Sim.dt_step;            % time step (s)
Tseq=Sim.t_Start:dt:Sim.t_End;

%--------------------
% PARAMETERS
%--------------------
% CELL
VEreset=He*Theta(1); % Since 1:14 are all excitatory
VIreset=Hi*Theta(end);
%
%----------
% STIMULUS
%----------
% build external current
% baseline current for all neurons
% add stimuli on top of baseline: for each stimulus provide
%              - profile (perc. increase on baseline current)
%              - index of selective neurons
%
% BASELINE EXTERNAL CURRENT
mu=Ext.Mu; % mu is an (N_e+N_i)-th array
if strcmp(Stimulus.option,'on')
    stim=Ext.stim;
    nstim=numel(stim); % number of stimuli in current trials
end
%----------------
% SYNAPTIC FILTER
%----------------

Tau.tausyn_e=tausyn_e; % exc synaptic time (fall time)
Tau.tausyn_i=tausyn_i; % exc synaptic time (fall time)
F=synaptic_trace(Tau,dt,N_e,N_i); % traces for recurrent connections

%--------------------------
% SIMULATION
%--------------------------
% preallocate memory for stored variable firings_tmp
% INITIAL CONDITIONS: random
v=[(Theta(1)-VEreset)/2*ones(N_e,1)+(Theta(1)-VEreset)/2*(2*rand(N_e,1)-1);...
    (Theta(end)-VIreset)/2*ones(N_i,1)+(Theta(end)-VIreset)/2*(2*rand(N_i,1)-1)];  % Initial values of v
% aa = 'vFix.mat';
% save(aa,'v');
% v = load('vFix.mat');
% v = v.v;
% THRESHOLD VECTOR

% v =v/3;
VTh=[Theta(1)*ones(N_e,1); Theta(end)*ones(N_i,1)];
c=[VEreset*ones(N_e,1);  VIreset*ones(N_i,1)]; % reset potentials
% fprintf('\nVEThres=%g --- VIThres=%g',VTh(1),VTh(end));
% fprintf('\nVEreset=%g --- VIreset=%g \n',c(1),c(end));
% Excitatory neurons        Inhibitory neurons
tau=[tau_e*ones(N_e,1);       tau_i*ones(N_i,1)];
%
firings=zeros(10*numel(Tseq),2);
firings_cnt=0;
% tic
%--------------------
% PLOT
%--------------------
PlotData=[];
PlotData.Ne_plot=N_e; % number of exc neuron to plot
PlotData.Ni_plot=N_i; % number of inh neurons to plot
ind_plot=[5; N_e+5]; % indices of neurons to plot
indcue=find(cellfun(@(x)~isempty(x),strfind({stim(:).name},'CS')));
if ~isempty(indcue)
    ind_plot(1)=stim(indcue).ind(1);
end
nplot=numel(ind_plot); % number of neurons to plot (membrane potential plot)
vi=0; % running index for vplot
PlotData.vplot = zeros(nplot,round(Sim.plot_length/dt)); % store membrane potential for plots; rows=neurons, cols=time steps;
PlotData.iEplot = zeros(2,round(Sim.plot_length/dt)); % store EPSC for plots; rows=neurons, cols=time steps;
PlotData.iExtplot = zeros(2,round(Sim.plot_length/dt)); % store IPSC for plots; rows=neurons, cols=time steps;
PlotData.iIplot = zeros(2,round(Sim.plot_length/dt)); % store IPSC for plots; rows=neurons, cols=time steps;
PlotData.p=p;
PlotData.VTh=VTh;
PlotData.tau=tau;
PlotData.ind_plot=ind_plot;
%----------------------------
% RUN
%----------------------------
NcUnits = 50;
popsize=repmat(NcUnits,10,1);
cusumNcE=[0 cumsum(popsize)'];

syncount = 0;
count3=1;
count4=1;
GstoreDec = cell(40,1);
Gstore = cell(40,1);
TimeStoreDec = [];
TimeStore = [];
Dstore = cell(40,1);
DstoreDec = cell(40,1);

JEI = J(1:500,501:1000);
JIE = J(501:1000,1:500);

% switch CondIter
%     case 1
%         jinc1 = 0;
%         jdec1 = 1;
%         jdec2 = 1;
%     case 2
%         jinc1 = 1;
%         jdec1 = 0;
%         jdec2 = 1;
%     case 3
%         jinc1 = 1;
%         jdec1 = 1;
%         jdec2 = 0;
%     case 4
%         jinc1 = 0.1;
%         jdec1 = 0.1;
%         jdec2 = 0.1;
%     case 5
%         jinc1 = 0;
%         jdec1 = 0.1;
%         jdec2 = 0.1;
%     case 6
%         jinc1 = 0.1;
%         jdec1 = 0;
%         jdec2 = 0.1;
%     case 7
%         jinc1 = 0.1;
%         jdec1 = 0.1;
%         jdec2 = 0;
%     case 8
%         jinc1 = 0.01;
%         jdec1 = 0.01;
%         jdec2 = 0.01;
%     case 9
%         jinc1 = 0;
%         jdec1 = 0.01;
%         jdec2 = 0.01;
%     case 10
%         jinc1 = 0.01;
%         jdec1 = 0;
%         jdec2 = 0.01;
%     case 11
%         jinc1 = 0.01;
%         jdec1 = 0.01;
%         jdec2 = 0;
% end

% switch CondIter
%     case 1
%         jinc1 = 0.1;
%         jdec1 = 0;
%         jdec2 = 0.1;
%    case 2
%         jinc1 = 0.1;
%         jdec1 = 0;
%         jdec2 = 0.05;
%    case 3
%         jinc1 = 0.05;
%         jdec1 = 0;
%         jdec2 = 0.1;
%    case 4
%         jinc1 = 0.1;
%         jdec1 = 0;
%         jdec2 = 0.01;
%    case 5
        jinc1 = 0.01;
        jdec1 = 0;
        jdec2 = 0.1;
%     case 6
%         jinc1 = 0.01;
%         jdec1 = 0;
%         jdec2 = 0.01;
%  end
%  switch CondIter
%     case 1
%         if rem(Odor,2) == 1
%             jinc1 = 0.01;
%             jdec1 = 0;
%             jdec2 = 0.1;
%         else
%             jinc1 = 0.01;
%             jdec1 = 0;
%             jdec2 = 0.01;
%         end
%         
%     case 2
%         if rem(Odor,2) == 1
%             jinc1 = 0.01;
%             jdec1 = 0;
%             jdec2 = 0.05;
%         else
%             jinc1 = 0.01;
%             jdec1 = 0;
%             jdec2 = 0.01;
%         end
%      case 1
%         if rem(Odor,2) == 1
%             jinc1 = 0.01;
%             jdec1 = 0;
%             jdec2 = 0.1;
%         else
%             jinc1 = 0.01;
%             jdec1 = 0;
%             jdec2 = 0.05;
%         end
%     case 4
%         if rem(Odor,2) == 1
%             jinc1 = 0.01;
%             jdec1 = 0;
%             jdec2 = 0.2;
%         else
%             jinc1 = 0.01;
%             jdec1 = 0;
%             jdec2 = 0.01;
%         end
%     case 5
%         if rem(Odor,2) == 1
%             jinc1 = 0.01;
%             jdec1 = 0;
%             jdec2 = 0.05;
%         else
%             jinc1 = 0.01;
%             jdec1 = 0;
%             jdec2 = 0.0;
%         end
%     case 6
%         if rem(Odor,2) == 1
%             jinc1 = 0.01;
%             jdec1 = 0;
%             jdec2 = 0.1;
%         else
%             jinc1 = 0.01;
%             jdec1 = 0;
%             jdec2 = 0.0;
%         end
%             case 7
%         jinc1 = 0.1;
%         jdec1 = 0.1;
%         jdec2 = 0.1;
%     case 8
%         jinc1 = 0;
%         jdec1 = 0.1;
%         jdec2 = 0.1;
%     case 2
%         jinc1 = 0.1;
%         jdec1 = 0;
%         jdec2 = 0.1;
%     case 10
%         jinc1 = 0.1;
%         jdec1 = 0.1;
%         jdec2 = 0;
%     case 11
%         jinc1 = 0.01;
%         jdec1 = 0.01;
%         jdec2 = 0.01;
%     case 12
%         jinc1 = 0;
%         jdec1 = 0.01;
%         jdec2 = 0.01;
%     case 13
%         jinc1 = 0.01;
%         jdec1 = 0;
%         jdec2 = 0.01;
%     case 14
%         jinc1 = 0.01;
%         jdec1 = 0.01;
%         jdec2 = 0;
%            case 15
%         jinc1 = 0.1;
%         jdec1 = 0;
%         jdec2 = 0.01;
%    case 16
%         jinc1 = 0.01;
%         jdec1 = 0;
%         jdec2 = 0.1;
%  end

threshclus = 8;

refr=zeros(size(mu,1),1);       % neurons in refractory state
for t=1:1000%numel(Tseq)         % siMulation of 1000 ms
    fired=find(v>VTh); % indices of spikes
    Isyn=zeros(N_e+N_i,1);
    % spikes
    if ~isempty(fired)
        v(fired)=c(fired);
        refr(fired)=tau_arp;
    end
    SeriousD = [];
    SeriousG = [];
    DClustAct = [];
    GClustAct = [];
    aa = [];
    yy = [];
    low = [];
    high = [];
    
    
    if t> 500
        syncount = syncount+1;
        high = find(nonzeros(firings(:,1))< Tseq(t));
        low=find(nonzeros(firings(:,1)) > Tseq(t-499));
        yy = intersect(low,high);
        aa = firings(yy,2);
        count1 =1;
        count2= 1;
        for kk =1:length(yy)
            if aa(kk)>500
                DClustAct(count1) = ceil(aa(kk)/50)-10;
                count1 = count1+1;
            else
                GClustAct(count2) = ceil(aa(kk)/50);
                count2 = count2+1;
            end
        end
        count1 =1;
        count2= 1;
        if ~isempty(DClustAct)
            %                 [a1,b1] = groupcounts(DClustAct);
            tb = tabulate(DClustAct);
            for kl =1:length(tb(:,1))
                if tb(kl,2)>threshclus
                    SeriousD(count1) = tb(kl);
                    count1= count1+1;
                end
            end
        end
        if ~isempty(GClustAct)
            %                 [a2,b2] = groupcounts(GClustAct);
            tb2 = tabulate(GClustAct);
            for kl =1:length(tb2(:,1))
                if tb2(kl,2)>threshclus
                    SeriousG(count2) = tb2(kl);
                    count2= count2+1;
                end
            end
        end
    end
    
    for km =1:length(SeriousD)
        tempG =  Inh(SeriousD(km),:);
        ID = intersect(SeriousG,tempG);
        ID2 = setdiff(tempG,SeriousG);
        ID2 = nonzeros(ID2)';
        ID3 = setdiff(SeriousG,tempG);
        
        %% Weight increase for overlap
        if ~isempty(ID)
            for kc = 1: length(ID)
                JEI(cusumNcE(ID(kc))+1:cusumNcE(ID(kc)+1),cusumNcE(SeriousD(km))+1:cusumNcE(SeriousD(km)+1)) =  JEI(cusumNcE(ID(kc))+1:cusumNcE(ID(kc)+1),cusumNcE(SeriousD(km))+1:cusumNcE(SeriousD(km)+1))*((1 + 0.0001*jinc1));
                JIE(cusumNcE(SeriousD(km))+1:cusumNcE(SeriousD(km)+1),cusumNcE(ID(kc))+1:cusumNcE(ID(kc)+1)) =  JIE(cusumNcE(SeriousD(km))+1:cusumNcE(SeriousD(km)+1),cusumNcE(ID(kc))+1:cusumNcE(ID(kc)+1))*((1 + 0.0001*jinc1));
            end
            %             TimeStore(count3) = Tseq(t);
            %             Gstore(count3) = {SeriousG};
            %             Dstore(count3) = {SeriousD};
            %             count3 = count3+1;
        end
        
        %% Weight Decrease for Active DAT and non active Glom
        if ~isempty(ID2)
            for kc = 1: length(ID2)
                JEI(cusumNcE(ID2(kc))+1:cusumNcE(ID2(kc)+1),cusumNcE(SeriousD(km))+1:cusumNcE(SeriousD(km)+1)) =  JEI(cusumNcE(ID2(kc))+1:cusumNcE(ID2(kc)+1),cusumNcE(SeriousD(km))+1:cusumNcE(SeriousD(km)+1))*((1 - 0.0001*jdec1));
                JIE(cusumNcE(SeriousD(km))+1:cusumNcE(SeriousD(km)+1),cusumNcE(ID2(kc))+1:cusumNcE(ID2(kc)+1)) =  JIE(cusumNcE(SeriousD(km))+1:cusumNcE(SeriousD(km)+1),cusumNcE(ID2(kc))+1:cusumNcE(ID2(kc)+1))*((1 - 0.0001*jdec1));
            end
        end
        
        %% Weight Decrease for Active Glom and non active DAT
        if ~isempty(ID3)
            for kc = 1: length(ID3)
                JEI(cusumNcE(ID3(kc))+1:cusumNcE(ID3(kc)+1),cusumNcE(SeriousD(km))+1:cusumNcE(SeriousD(km)+1)) =  JEI(cusumNcE(ID3(kc))+1:cusumNcE(ID3(kc)+1),cusumNcE(SeriousD(km))+1:cusumNcE(SeriousD(km)+1))*((1 - 0.0001*jdec2));
                JIE(cusumNcE(SeriousD(km))+1:cusumNcE(SeriousD(km)+1),cusumNcE(ID3(kc))+1:cusumNcE(ID3(kc)+1)) =  JIE(cusumNcE(SeriousD(km))+1:cusumNcE(SeriousD(km)+1),cusumNcE(ID3(kc))+1:cusumNcE(ID3(kc)+1))*((1 - 0.0001*jdec2));
            end
        end
        
        
        %             TimeStoreDec(count4) = Tseq(t);
        %             GstoreDec(count4) = {SeriousG};
        %             DstoreDec(count4) = {SeriousD};
        %             count4 = count4+1;
        %disp('weight decrease');
    end
    
    J=[J(1:500,1:500) JEI; JIE J(501:1000,501:1000)];
    
    % recurrent synaptic current
    F=syn_evolve(F,fired);
    % integrate
    muRun=mu;
    if strcmp(Stimulus.option,'on')
        for n=1:nstim
            if Tseq(t)>=stim(n).interval(1) && Tseq(t)<=stim(n).interval(2) % once inside cue interval
                if strcmp(stim(n).name,'CSgauss') || strcmp(stim(n).name,'hab')
                    muRun(stim(n).ind)=muRun(stim(n).ind)+stim(n).profile(Tseq(t))*mu(stim(n).ind).*stim(n).gauss; % update with cue
                else
                    muRun(stim(n).ind)=muRun(stim(n).ind)+stim(n).profile(Tseq(t))*mu(stim(n).ind); % update if no cue?
                end
            end
        end
    end
    Isyn=Isyn+J*F.f;
    v=v-v*dt./tau+muRun(:)*dt+Isyn*dt;
    % neurons in refractory state
    refr=max(-0.001,refr-dt);
    v(refr>0)=c(refr>0);
    % store spikes
    if ~isempty(fired)
        % if firings_tmp has no more space, preallocate more memory
        if firings_cnt+numel(fired)>size(firings,1)
            firings=[firings; zeros(10*numel(Tseq),2)];
        end
        firings(firings_cnt+1:firings_cnt+numel(fired),1:2)=[Tseq(t)+0*fired, fired];
        firings_cnt=firings_cnt+numel(fired);
    end
    % store values for plotting, only last Sim.plot_length interval
    if Tseq(t)>Sim.t_End-Sim.plot_length
        vi=vi+1;
        % membrane potential
        PlotData.vplot(1:nplot,vi)=v(ind_plot);
        % input currents
        PlotData.iEplot(1:nplot,vi)=J(ind_plot,1:N_e)*F.f(1:N_e);
        PlotData.iIplot(1:nplot,vi)=J(ind_plot,N_e+1:N_e+N_i)*F.f(N_e+1:N_e+N_i);
        PlotData.iExtplot(1:nplot,vi)=muRun(ind_plot,1);
    end
    
    
    
    
end

% fprintf('--- End of trial...\n');
% toc
%---------------------------------------------------------------------------
if ~any(any(firings))
    fprintf('\n --- NO SPIKES GENERATED... \n');
else
    % find last spike in firings
    IndexEnd=find(firings(:,2)==0,1)-1;
    if isempty(IndexEnd)
        IndexEnd=size(firings,1);
    end
    all_firings=firings(1:IndexEnd,[1 2]);
end


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

    function F=syn_evolve(F,fired)
        
        
        % update synaptic filter
        F.f=F.fexp.*F.f;
        if ~isempty(fired)
            % update of synaptic filter with spikes
            F.f(fired)=F.f(fired)+F.fSpike(fired);
        end
    end


end