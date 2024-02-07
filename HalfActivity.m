Params.Sim.plot_length=Params.Sim.plot_length;
% save plots
savedir=Params.savedir;
if ~exist(savedir,'dir'); mkdir(savedir); end
filesave=fullfile(savedir,sprintf('SIM_[Jp%0.03g]_[p%d]',...
    Params.Jplus,Params.p));
%
Ne_plot=data.Ne_plot; % number of exc neuron to plot (typically all)
Ni_plot=data.Ni_plot; % number of inh neurons to plot  (typically all)
N_e=Ne_plot;
N_i=Ni_plot;
ind_plot=data.ind_plot;%[5; N_e+5]; % indices of neurons to plot
vplot=data.vplot; % store membrane potential for plots; rows=neurons, cols=time steps;
iEplot=data.iEplot; % store EPSC for plots; rows=neurons, cols=time steps;
iExtplot=data.iExtplot; % store IPSC for plots; rows=neurons, cols=time steps;
iIplot=data.iIplot; % store IPSC for plots; rows=neurons, cols=time steps;
VTh=data.VTh;
Jee=Params.Jee;
tau=data.tau;
p=data.p;
%
%
Ext=Params.Ext;
Network=Params.Network;
Sim=Params.Sim;
if any(strcmp(fieldnames(Params),'NcE'))
    NcE=Params.NcE;
elseif any(strcmp(fieldnames(Params),'popsize'))
    NcE=Params.popsize;
end
% keep only pops with at least 1 neuron
NcE=NcE(NcE>0);
clustermatrix=Params.clustermatrix;
clustermatrix=clustermatrix(1:numel(NcE),:);
cusumNcE=[0 cumsum(NcE)'];
fprintf('only pops sharing up to %d clusters are populated.\n',sum(clustermatrix(end,:)));

XLIM=[Sim.t_Start Sim.t_End];
if Sim.t_Start<-5
    XLIM(1)=-5;
end
if Sim.t_End>5
    XLIM(2)=5;
end

BinSize=0.05;
bins=Sim.t_End-Sim.plot_length:BinSize:Sim.t_End;

Params.Sim.plot_length=Params.Sim.plot_length;
% save plots
savedir=Params.savedir;
if ~exist(savedir,'dir'); mkdir(savedir); end
filesave=fullfile(savedir,sprintf('SIM_[Jp%0.03g]_[p%d]',...
    Params.Jplus,Params.p));
%

Ne_plot=data.Ne_plot; % number of exc neuron to plot (typically all)
Ni_plot=data.Ni_plot; % number of inh neurons to plot  (typically all)
N_e=Ne_plot;
N_i=Ni_plot;
ind_plot=data.ind_plot;%[5; N_e+5]; % indices of neurons to plot
vplot=data.vplot; % store membrane potential for plots; rows=neurons, cols=time steps;
iEplot=data.iEplot; % store EPSC for plots; rows=neurons, cols=time steps;
iExtplot=data.iExtplot; % store IPSC for plots; rows=neurons, cols=time steps;
iIplot=data.iIplot; % store IPSC for plots; rows=neurons, cols=time steps;
VTh=data.VTh;
Jee=Params.Jee;
tau=data.tau;
p=data.p;
%
%
Ext=Params.Ext;
Network=Params.Network;
Sim=Params.Sim;
if any(strcmp(fieldnames(Params),'NcE'))
    NcE=Params.NcE;
elseif any(strcmp(fieldnames(Params),'popsize'))
    NcE=Params.popsize;
end
% keep only pops with at least 1 neuron
NcE=NcE(NcE>0);
clustermatrix=Params.clustermatrix;
clustermatrix=clustermatrix(1:numel(NcE),:);
cusumNcE=[0 cumsum(NcE)'];
fprintf('only pops sharing up to %d clusters are populated.\n',sum(clustermatrix(end,:)));

XLIM=[Sim.t_Start Sim.t_End];
if Sim.t_Start<-5
    XLIM(1)=-5;
end
if Sim.t_End>5
    XLIM(2)=5;
end

for clu=1:numel(NcE)
    ind_Neu(clu) = round(0.5*rand(1,1)*104 +cusumNcE(clu));
end

for iTrial = 1:20 % pick which trial to plot
    
    data=dataload.PlotData{iTrial}; % membrane potentials
    firings=dataload.firings{iTrial}; % spikes
    
    MinRate=10;
    BinSize=diff(bins(1:2));
    fntsz=15;
    cols=aux.distinguishable_colors(numel(NcE)+1);
    binwidth=8;
    cusumNcE=[0 cumsum(NcE)'];
    rate=zeros(numel(NcE),numel(bins)+1);
    figure;
    for clu=1:numel(NcE)
        ind_sp = firings(find(firings(:,2) == ind_Neu(clu)),1);
        x(:,clu,iTrial)=histc(ind_sp,bins);
        [x_bins, zold(:,clu,iTrial)]=aux.gaussfilt(bins,x(:,clu,iTrial),binwidth,0);
        
        rate(clu,:)=zold(:,clu,iTrial)/(NcE(clu)*BinSize);
        plot(x_bins-diff(x_bins(1:2))/2,rate(clu,:),'color',cols(clu,:),'linewidth',1.5);
        hold on
    end
    hold off
end

for i = 1:37
    x2(i,:,:)  = mean(x(i:i+3,:,:));
end
for iTrial = 1:20
    for clu=1:numel(NcE)
        [tout, z(:,clu,iTrial)]=aux.gaussfilt(bins(1:37),x2(:,clu,iTrial),binwidth,0);
    end
end

ActivationFrame=zeros(4,20);
distances=zeros(1,20);
peak =[];

count = 0;
for i = 1:20
    count = 1;
    for j = 1:9
        if j == 1 || j == 2 || j == 4 || j == 9
            distances=zeros(1,21);
            [peak(count,i),binn(count,i)] = max(z(18:38,j,i));
            slope = gradient(z(18:38,j,i));
            if peak(count,i)>0.2                   %%do a closed range rather than an open range which leads to redundacy
            for k=1:binn(count,i)
                if slope(k)>0
                    distances(k)=abs(peak(count,i)/2-z(k+17,j,i));
                end
            end
            for m =1:length(distances)
                if distances(m)==0
                    distances(m)=10;                                %if distance is 0, something went wrong set distance to high number so that it isn't captured later
                end
            end
            [val,idx]= min(distances(:));               %Selecting the Index for minimum distance only if the difference between actual half peak and expected half peak is less than 1%. Since it isn't a continuous vector...
            %                 if val<0.01
            ActivationFrame(count,i)=idx;
            end
            count=count+1;
            %                 else
            %                     ActivationFrame(j,k)=0;
            %                 end
            %             else
            %                 ActivationFrame(j,k)= 0;
            %             end
        end
    end
end

for i =1:20
    ActivationFrameAvg(i) = mean(nonzeros(ActivationFrame(:,i)));            
end        
            
ActivationFrameAvg2 = ActivationFrameAvg;
% 
[h,p] = ttest(ActivationFrameAvg2', ActivationFrameAvg')

mean(ActivationFrameAvg)