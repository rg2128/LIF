
% SIM of one trials given parameters
%
% Luca Mazzucato March 2014

% SET OPTIONS
% ParamsRun = structure containing parameters for simulation


function [all_firings, PlotData, StimulusFirings,J]=fun_LIF_SIM(ParamsRun,Inh,i,m,j,k,rr2, Odor1, Odor2)

Odor = i;
Tone =j;
CondIter = m;
Theta=ParamsRun.Theta; Sim=ParamsRun.Sim; Stimulus=ParamsRun.Stimulus;
Ext=ParamsRun.Ext; J=ParamsRun.J; N_e=ParamsRun.N_e; N_i=ParamsRun.N_i; p=ParamsRun.p; He=ParamsRun.He;
Hi=ParamsRun.Hi; tau_e=ParamsRun.tau_e; tau_i=ParamsRun.tau_i; tausyn_e=ParamsRun.tausyn_e;
tausyn_i=ParamsRun.tausyn_i; tau_arp=ParamsRun.tau_arp;
%
all_firings=[];
dt=Sim.dt_step;            % time step (s)
Tseq=Sim.t_Start:dt:Sim.t_End;

tgain = k;
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
NcUnits = 10;
popsize=repmat(NcUnits,50,1);
cusumNcE=[0 cumsum(popsize)'];


JEI = J(1:500,501:1000);
JIE = J(501:1000,1:500);

% 
jinc1 = 0.05;
taur = 0.01;
Tthresh = 0.002;
Tmax = 0.02;

% switch CondIter
%     case 1
%         Tthresh = 0.01;
%         Tmax = 0.015;
%     case 2
%         Tthresh = 0.01;
%         Tmax = 0.02;
%     case 3
%         Tthresh = 0.01;
%         Tmax = 0.03;
%     case 4
%         Tthresh = 0.01;
%         Tmax = 0.04;
%     case 5
%         Tthresh = 0.01;
%         Tmax = 0.05;
% end
% rr2 = 0;
% while rr2<0.5 && rr2>-0.5
%     rr2 = (-10 + (20)*rand(1,1))/10;
% end

% if CondIter == 1 || CondIter == 2 || CondIter == 3 || CondIter == 4
%     bgain = 1; 
% elseif CondIter == 5 || CondIter == 6 || CondIter == 7 || CondIter == 8  
%     bgain = 1.5; 
% elseif CondIter == 9 || CondIter == 10 || CondIter == 11 || CondIter == 12      
     bgain = 2; 
% elseif CondIter == 13 || CondIter == 14 || CondIter == 15 || CondIter == 16     
%     bgain = 2.5; 
% elseif CondIter == 17 || CondIter == 18 || CondIter == 19 || CondIter == 20  
%     bgain = 3;
% end
    
% else
%     bgain = 4;
% end
%     bgain = 2;
% elseif tgain == 2 || tgain == 6 || tgain == 10 || tgain == 14  
%     bgain = 3;
% elseif tgain == 3 || tgain == 7 || tgain == 11 || tgain == 15  
%     bgain = 4;
% elseif tgain == 4 || tgain == 8 || tgain == 12 || tgain == 16  
%     bgain = 5;
% end
    
% prune = 0.5;
% % bgain = 4;
% tau_i_Thresh = 0.0207;
% tau_i_tc = 6;
% Taufile = 'Taufix.mat';
% formatSpec = 'Simulations/k%d';
% file_name = sprintf(formatSpec,tgain);
% 
% if rem(Odor,2) == 0 %&& Odor > 3
%     for dat = 1: length(stim(2).ind)
%         glom = ceil(stim(2).ind(dat)/50);
%         [xD,~] = find(Inh == glom);
%         for dat2 =1: length(xD)
%             allDAT = 500+(xD(dat2)-1)*10+1: 500 + xD(dat2)*10;
%             a=randperm(numel(allDAT(:)));
%             sparseDAT = allDAT(a(1:ceil(prune*numel(allDAT(:)))));
% %             tau(sparseDAT) = tau_i_Thresh - (tau_i_Thresh-tau_i)*exp((4-Odor)/tau_i_tc);
%             tau(sparseDAT) = tau_i_Thresh - (tau_i_Thresh-tau_i)*exp((2-Odor)/tau_i_tc);
%         end
%     end
%     parsave(fullfile(pwd,file_name,Taufile),tau);
% elseif rem(Odor,2) == 1 && Odor > 2
%     tau = par_load(fullfile(pwd, file_name,Taufile));
%     parsave(fullfile(pwd,file_name,Taufile),tau);
%     tau = tau.x;
% end

switch CondIter
    case 1
        synplus = 10;
        synminus = 26;
    case 2
        synplus = 10;
        synminus = 26;
    case 3
        synplus = 10;
        synminus = 26;
%     case 3
%         synplus = 10;
%         synminus = 26;
end
% 
%     case 3
%           synplus = 0.001;
% %         synminus = 0.001;
% end
%     case 5 
%         synplus = 0.1;
%         synminus = 1;
%     case 6 
%         synplus = 0.1;
%         synminus = 0.1;       
% end
% AllJ = zeros(10,10,12000);
% temptime = 1;
% allfct = cell(1,numel(Tseq) );
refr=zeros(size(mu,1),1);       % neurons in refractory state
for t=2500:numel(Tseq)-2500         % siMulation of 1000 ms
    fired=find(v>VTh); % indices of spikes
    Isyn=zeros(N_e+N_i,1);
    % spikes
    if ~isempty(fired)
        v(fired)=c(fired);
        refr(fired)=tau_arp;
    end
 
    % recurrent synaptic current
    F=syn_evolve(F,fired);
    % integrate
    muRun=mu;
    if strcmp(Stimulus.option,'on')
        for n=1:nstim
            if Tseq(t)>=stim(n).interval(1) && Tseq(t)<=stim(n).interval(2) % once inside cue interval
                if  strcmp(stim(n).name,'hab')
                    muRun(stim(n).ind)=muRun(stim(n).ind)+stim(n).profile(Tseq(t))*mu(stim(n).ind).*stim(n).gauss; % update with cue
                elseif strcmp(stim(n).name,'CSgauss')
                    if rem(Tone,2) == 0
                                            muRun(stim(n).ind)=muRun(stim(n).ind)+stim(n).profile(Tseq(t))*mu(stim(n).ind).*stim(n).gauss; % update with cue
                    end
%                         if t<5000
%                             muRun(stim(n).ind)=muRun(stim(n).ind)+stim(n).profile(rr2*2*3.14*Tseq(t))*mu(stim(n).ind).*stim(n).gauss;
%                         elseif t<10000 && t>4999
%                             muRun(stim(n).ind)=muRun(stim(n).ind)+stim(n).profile(rr2*2*bgain*3.14*Tseq(t))*mu(stim(n).ind).*stim(n).gauss;
%                             elseif t<15000 && t>9999
%                             muRun(stim(n).ind)=muRun(stim(n).ind)+stim(n).profile(rr2*3*bgain*3.14*Tseq(t))*mu(stim(n).ind).*stim(n).gauss;
%                         elseif t>14999
%                             muRun(stim(n).ind)=muRun(stim(n).ind)+stim(n).profile(rr2*2*3.14*Tseq(t))*mu(stim(n).ind).*stim(n).gauss;
%                         end
%                 else
%                         if t<10000
%                             muRun(stim(n).ind)=muRun(stim(n).ind)+stim(n).profile(rr2*2*3.14*Tseq(t))*mu(stim(n).ind).*stim(n).gauss;
%                         elseif t<15000 && t>9999
%                             muRun(stim(n).ind)=muRun(stim(n).ind)+stim(n).profile(rr2*3*bgain*3.14*Tseq(t))*mu(stim(n).ind).*stim(n).gauss;
%                         elseif t>14999
%                             muRun(stim(n).ind)=muRun(stim(n).ind)+stim(n).profile(rr2*2*3.14*Tseq(t))*mu(stim(n).ind).*stim(n).gauss;
%                         end 
%                     end
                elseif strcmp(stim(n).name,'US')
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
%     if Tseq(t)>Sim.t_End-Sim.plot_length
%         vi=vi+1;
%         % membrane potential
%         PlotData.vplot(1:nplot,vi)=v(ind_plot);
%         % input currents
%         PlotData.iEplot(1:nplot,vi)=J(ind_plot,1:N_e)*F.f(1:N_e);
%         PlotData.iIplot(1:nplot,vi)=J(ind_plot,N_e+1:N_e+N_i)*F.f(N_e+1:N_e+N_i);
%         PlotData.iExtplot(1:nplot,vi)=muRun(ind_plot,1);
%     end
    
%     low = [];
%     high = [];
%     temptime = [];
%     if ~isempty(fired)
%         for neuronsfired = 1:length(fired)
%             if fired(neuronsfired)>500
%                 tempf = firings_cnt-1;
% %               % [cl,in] = min(abs(firings(tempf,1) - Tmax, firings(1:tempf,1))) 
% %                 [~,inds]=min(unique((abs(firings(1:tempf,1)-(firings(tempf,1) - Tmax))),'stable'));
% 
%                 low = find((firings(1:tempf,1) + Tthresh)<Tseq(t));
%                 high = find((firings(1:tempf,1) + Tmax)<Tseq(t));
%                 temptime = setdiff(low,high);
%                 
% %                 temptime = inds:tempf;
%                 if ~isempty(temptime)
%                     ExcOnly = firings((firings(temptime,2)<500),:);
%                     dg4 = double(ismember(ceil(ExcOnly(:,2)/10),Inh(ceil(fired(neuronsfired)/10 -50),:)));
%                     for dg3 =1:length(ExcOnly)
%                         if dg4(dg3)>0
%                             JEI(cusumNcE(ceil(ExcOnly(dg3,2)/10))+1:cusumNcE(ceil(ExcOnly(dg3,2)/10)+1),cusumNcE(ceil(fired(neuronsfired)/10)-50)+1:cusumNcE((ceil(fired(neuronsfired)/10)-50)+1)) =  JEI(cusumNcE(ceil(ExcOnly(dg3,2)/10))+1:cusumNcE(ceil(ExcOnly(dg3,2)/10)+1),cusumNcE(ceil(fired(neuronsfired)/10)-50)+1:cusumNcE((ceil(fired(neuronsfired)/10)-50)+1))*((1 + 0.01*jinc1*exp(-abs(Tseq(t)- ExcOnly(dg3,1))/taur)));
%                             JIE(cusumNcE(ceil(fired(neuronsfired)/10)-50)+1:cusumNcE((ceil(fired(neuronsfired)/10)-50)+1),cusumNcE(ceil(ExcOnly(dg3,2)/10))+1:cusumNcE(ceil(ExcOnly(dg3,2)/10)+1)) =  JIE(cusumNcE(ceil(fired(neuronsfired)/10)-50)+1:cusumNcE((ceil(fired(neuronsfired)/10)-50)+1),cusumNcE(ceil(ExcOnly(dg3,2)/10))+1:cusumNcE(ceil(ExcOnly(dg3,2)/10)+1)) *((1 + 0.01*jinc1*exp(-abs(Tseq(t)- ExcOnly(dg3,1))/taur)));
%                         end
%                     end
%                 else
%                         for dg3 = Inh(ceil(fired(neuronsfired)/10 -50),:)
%                                 JEI(cusumNcE(ceil(ExcOnly(dg3,2)/10))+1:cusumNcE(ceil(ExcOnly(dg3,2)/10)+1),cusumNcE(ceil(fired(neuronsfired)/10)-50)+1:cusumNcE((ceil(fired(neuronsfired)/10)-50)+1)) =  JEI(cusumNcE(ceil(ExcOnly(dg3,2)/10))+1:cusumNcE(ceil(ExcOnly(dg3,2)/10)+1),cusumNcE(ceil(fired(neuronsfired)/10)-50)+1:cusumNcE((ceil(fired(neuronsfired)/10)-50)+1))*((1 + 0.01*jinc1*exp(-abs(Tseq(t)- ExcOnly(dg3,1))/taur)));
%                                 JIE(cusumNcE(ceil(fired(neuronsfired)/10)-50)+1:cusumNcE((ceil(fired(neuronsfired)/10)-50)+1),cusumNcE(ceil(ExcOnly(dg3,2)/10))+1:cusumNcE(ceil(ExcOnly(dg3,2)/10)+1)) =  JIE(cusumNcE(ceil(fired(neuronsfired)/10)-50)+1:cusumNcE((ceil(fired(neuronsfired)/10)-50)+1),cusumNcE(ceil(ExcOnly(dg3,2)/10))+1:cusumNcE(ceil(ExcOnly(dg3,2)/10)+1)) *((1 + 0.01*jinc1*exp(-abs(Tseq(t)- ExcOnly(dg3,1))/taur)));
%                            end
%                     end
%                 
%                 end      
%             end           
%         end         
%     end
    
    
    
%     if t>10000 &&t<17001
%         aJ = [];
%         countJ = 1;
%         j1 = 1;
%         for jJ =1:50:500
%             countJ = 1;
%             for iJ =1:50:500
%                 aJ(j1,countJ) = mean(mean(JIE(iJ:iJ+49,jJ:jJ+49)));
%                 countJ = countJ+1;
%             end
%             j1 = j1 +1;
%         end
%         AllJ(:,:,t-10000) = aJ;
%     end
% toc
end
go = find(Odor1==1);
  for q=1:length(go)
        goAll(q,:) = go(q)*10-9:go(q)*10;
  end
    
  nogo = find(Odor2==1);
  for q=1:length(nogo)
        nogoAll(q,:) = nogo(q)*10-9:nogo(q)*10;
  end
  if rem(Odor,2) == 1
    OdorCurrent = go;
  else
    OdorCurrent = nogo;
  end

% if (firings(:,1)>0) && (firings(:,1)<0.6)
% ExcOnly=  firings(intersect(find(firings(:,1)>0),find(firings(:,1)<0.6)),2)<500;
% firingstemp = firings(intersect(find(firings(:,1)>0),find(firings(:,1)<0.6)),:);
% for pqr=1:length(ExcOnly)
%     if ExcOnly(pqr) && firingstemp(pqr,2) ~= 0
%         [datx,~] =  find(Inh == ceil(firingstemp(pqr,2)/10));
%         [~,idz] = (min(abs(firingstemp(pqr,1) + 0.01-firingstemp(:,1))));
%         NoDAT = 0;
%         for qrs = pqr+1:idz
%             if any(ceil(firingstemp(qrs,2)/10)-50 == datx)
%                 JEI(cusumNcE(ceil(firingstemp(pqr,2)/10))+1:cusumNcE(ceil(firingstemp(pqr,2)/10)+1),cusumNcE(ceil(firingstemp(qrs,2)/10)-50)+1:cusumNcE((ceil(firingstemp(qrs,2)/10)-50)+1)) =  JEI(cusumNcE(ceil(firingstemp(pqr,2)/10))+1:cusumNcE(ceil(firingstemp(pqr,2)/10)+1),cusumNcE(ceil(firingstemp(qrs,2)/10)-50)+1:cusumNcE((ceil(firingstemp(qrs,2)/10)-50)+1))*((1 + synplus*jinc1*exp(-abs(firingstemp(qrs,1)- firingstemp(pqr,1))/taur)));
%                 JIE(cusumNcE(ceil(firingstemp(qrs,2)/10)-50)+1:cusumNcE((ceil(firingstemp(qrs,2)/10)-50)+1),cusumNcE(ceil(firingstemp(pqr,2)/10))+1:cusumNcE(ceil(firingstemp(pqr,2)/10)+1)) =  JIE(cusumNcE(ceil(firingstemp(qrs,2)/10)-50)+1:cusumNcE((ceil(firingstemp(qrs,2)/10)-50)+1),cusumNcE(ceil(firingstemp(pqr,2)/10))+1:cusumNcE(ceil(firingstemp(pqr,2)/10)+1)) *((1 + synplus*jinc1*exp(-abs(firingstemp(qrs,1)- firingstemp(pqr,1))/taur)));
%                 NoDAT =1;
%             end
%         end
%         if NoDAT == 0
%             for rst =datx
%                 JEI(cusumNcE(ceil(firingstemp(pqr,2)/10))+1:cusumNcE(ceil(firingstemp(pqr,2)/10)+1),cusumNcE(rst)+1:cusumNcE((rst)+1)) =  JEI(cusumNcE(ceil(firingstemp(pqr,2)/10))+1:cusumNcE(ceil(firingstemp(pqr,2)/10)+1),cusumNcE(rst)+1:cusumNcE((rst)+1))*((1 - synminus*jinc1*exp(-abs(0.0001)/taur)));
%                 JIE(cusumNcE(rst)+1:cusumNcE((rst)+1),cusumNcE(ceil(firingstemp(pqr,2)/10))+1:cusumNcE(ceil(firingstemp(pqr,2)/10)+1)) =  JIE(cusumNcE(rst)+1:cusumNcE((rst)+1),cusumNcE(ceil(firingstemp(pqr,2)/10))+1:cusumNcE(ceil(firingstemp(pqr,2)/10)+1)) *((1 - synminus*jinc1*exp(-abs(0.0001)/taur)));
%             end
%         end
%     end
% end
% end
if rem(Tone,2) == 1
ExcOnly=  firings(intersect(find(firings(:,1)>0.1),find(firings(:,1)<0.5)),2)<500;
firingstemp = firings(intersect(find(firings(:,1)>0.1),find(firings(:,1)<0.5)),:);
addcount=0;
addcounttemp = 0;
whichglom = [];
addcountng = 0;
addcountg = 0;
for pqr=1:length(ExcOnly)
    if ExcOnly(pqr) && firingstemp(pqr,2) ~= 0
        [datx,~] =  find(Inh == ceil(firingstemp(pqr,2)/10));
        alldat = [];
%         for dalength =1:length(datx)
%             alldat = [alldat cusumNcE(datx(dalength))+1:cusumNcE(datx(dalength)+1)];
%         end

        [~,idz] = (min(abs(firingstemp(pqr,1) + 0.05-firingstemp(:,1))));
        for qrs = pqr+1:idz
            if any(ceil(firingstemp(qrs,2)/10)-50 == datx)
                deltat=abs(firingstemp(qrs,1)-firingstemp(pqr,1));
                
                if deltat>0.03
%                 JEI(cusumNcE(ceil(firingstemp(pqr,2)/10))+1:cusumNcE(ceil(firingstemp(pqr,2)/10)+1),alldat) =  JEI(cusumNcE(ceil(firingstemp(pqr,2)/10))+1:cusumNcE(ceil(firingstemp(pqr,2)/10)+1),alldat)*((1 + synplus*jinc1*exp(-abs(firingstemp(qrs,1)- firingstemp(pqr,1))/taur)));
%                 JIE(alldat,cusumNcE(ceil(firingstemp(pqr,2)/10))+1:cusumNcE(ceil(firingstemp(pqr,2)/10)+1)) =  JIE(alldat,cusumNcE(ceil(firingstemp(pqr,2)/10))+1:cusumNcE(ceil(firingstemp(pqr,2)/10)+1)) *((1 + synplus*jinc1*exp(-abs(firingstemp(qrs,1)- firingstemp(pqr,1))/taur)));  
                    
                
%                if ismember(ceil(firingstemp(pqr,2)/10),OdorCurrent)
                JEI(cusumNcE(ceil(firingstemp(pqr,2)/10))+1:cusumNcE(ceil(firingstemp(pqr,2)/10)+1),cusumNcE(ceil(firingstemp(qrs,2)/10)-50)+1:cusumNcE((ceil(firingstemp(qrs,2)/10)-50)+1)) =  JEI(cusumNcE(ceil(firingstemp(pqr,2)/10))+1:cusumNcE(ceil(firingstemp(pqr,2)/10)+1),cusumNcE(ceil(firingstemp(qrs,2)/10)-50)+1:cusumNcE((ceil(firingstemp(qrs,2)/10)-50)+1))*((1 + synplus*jinc1*exp(-abs(firingstemp(qrs,1)- firingstemp(pqr,1))/taur)));
                JIE(cusumNcE(ceil(firingstemp(qrs,2)/10)-50)+1:cusumNcE((ceil(firingstemp(qrs,2)/10)-50)+1),cusumNcE(ceil(firingstemp(pqr,2)/10))+1:cusumNcE(ceil(firingstemp(pqr,2)/10)+1)) =  JIE(cusumNcE(ceil(firingstemp(qrs,2)/10)-50)+1:cusumNcE((ceil(firingstemp(qrs,2)/10)-50)+1),cusumNcE(ceil(firingstemp(pqr,2)/10))+1:cusumNcE(ceil(firingstemp(pqr,2)/10)+1)) *((1 + synplus*jinc1*exp(-abs(firingstemp(qrs,1)- firingstemp(pqr,1))/taur)));  
               if ismember(ceil(firingstemp(pqr,2)/10),nogo)
                    addcountng= addcountng+1;
                elseif ismember(ceil(firingstemp(pqr,2)/10),go)
                    addcountg= addcountg+1;
               end
                
               addcount= addcount+1;
               whichglom(addcount) = firingstemp(pqr,1);
%                break;
               end

            end
        end
    
%     if addcounttemp ~= addcount
%     firingstemp(pqr,2) = 0;
%     end
%         addcounttemp = addcount;
    end
end
% aaaa = ceil(firingstemp(pqr,2)/10)
% for q=1:length(go)
% length(find(go(q) == aaaa))
%     end
% find(17 == aaaa)



InhOnly=  firings(intersect(find(firings(:,1)>0.1),find(firings(:,1)<0.5)),2)>500;
% firingstemp = firings(intersect(find(firings(:,1)>0.1),find(firings(:,1)<0.5)),:);
subtractcount=0;
subcountng = 0;
subcountg = 0;
for pqr=1:length(InhOnly)
    if InhOnly(pqr) && firingstemp(pqr,2) ~= 0
        glomx = Inh((ceil(firingstemp(pqr,2)/10)-50),:);
        glomx(glomx == 0) = [];
        [~,idz] = (min(abs(firingstemp(pqr,1) + 0.05-firingstemp(:,1))));
        for qrs = pqr+1:idz
            if any(ceil(firingstemp(qrs,2)/10) == glomx)
                deltat=abs(firingstemp(qrs,1)-firingstemp(pqr,1));
                if deltat>0.03
%                 if ismember(ceil(firingstemp(qrs,2)/10),OdorCurrent)
                    JEI(cusumNcE(ceil(firingstemp(qrs,2)/10))+1:cusumNcE(ceil(firingstemp(qrs,2)/10)+1),cusumNcE(ceil(firingstemp(pqr,2)/10)-50)+1:cusumNcE((ceil(firingstemp(pqr,2)/10)-50)+1)) =  JEI(cusumNcE(ceil(firingstemp(qrs,2)/10))+1:cusumNcE(ceil(firingstemp(qrs,2)/10)+1),cusumNcE(ceil(firingstemp(pqr,2)/10)-50)+1:cusumNcE((ceil(firingstemp(pqr,2)/10)-50)+1))*((1 - synminus*jinc1*exp(-abs(firingstemp(qrs,1)- firingstemp(pqr,1))/taur)));
                    JIE(cusumNcE(ceil(firingstemp(pqr,2)/10)-50)+1:cusumNcE((ceil(firingstemp(pqr,2)/10)-50)+1),cusumNcE(ceil(firingstemp(qrs,2)/10))+1:cusumNcE(ceil(firingstemp(qrs,2)/10)+1)) =  JIE(cusumNcE(ceil(firingstemp(pqr,2)/10)-50)+1:cusumNcE((ceil(firingstemp(pqr,2)/10)-50)+1),cusumNcE(ceil(firingstemp(qrs,2)/10))+1:cusumNcE(ceil(firingstemp(qrs,2)/10)+1)) *((1 - synminus*jinc1*exp(-abs(firingstemp(qrs,1)- firingstemp(pqr,1))/taur)));
                 
                if ismember(ceil(firingstemp(qrs,2)/10),nogo)
                   subcountng= subcountng+1;
                elseif ismember(ceil(firingstemp(qrs,2)/10),go)
                    subcountg= subcountg+1;
                end
                end
            end
        end
    end
end
end
% go = find(Odor1==1);
%   for q=1:length(go)
%         goAll(q,:) = go(q)*10-9:go(q)*10;
%   end
%     
%   nogo = find(Odor2==1);
%   for q=1:length(nogo)
%         nogoAll(q,:) = nogo(q)*10-9:nogo(q)*10;
%     end
% sum(sum(JIE(:,goAll)))
% sum(sum(JIE(:,nogoAll)))

J=[J(1:500,1:500) JEI; JIE J(501:1000,501:1000)];

%% Choosing stimulus specific firings
% tic
cf = 1;
StimulusFirings = [];
for a1 = 1:length(firings)
    if any(Ext.stim(2).ind == firings(a1,2))
        StimulusFirings(cf,:) = firings(a1,:);
        cf = cf +1;
    end
end
% toc
% J2= J;


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