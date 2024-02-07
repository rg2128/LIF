% Creates parameter files from choice of parameter set specified in Opt
% and saves it in .../DATA/Params.mat


% Repeat ketamine injection in this stimulation. Inhibition to initiate Ca2+ level of glomeruli firing of spontaneous
% activity associated
%COMBINED LEARNING CURVES TO SHOW dat inhibition at assessment  still leads to less performance level on average in those trials. but not with the training level of the same mice on saline


function [params2, Inh,rr2, Odor1, Odor2, OverlappedDAT] = create_params_EI(condition,k,j,i,m,imjk)
%----------
% NETWORKS
%----------
% Start and end of trial (in units of seconds)
StimulusSelection = i;
Sim.t_Start=-1;
Sim.t_End=1;
Sim.dt_step=0.0001; % integration step (s)
% StimulusSelection = i;
CountDrac = k;
CountDrac2 = j;

EIC = m;
%------------------
% NETWORK OPTIONS
%------------------
Network.clusters={'EI','IE','II'}; % all weights are clustered
% Network.clust='hom'; % homogeneous EE cluster size
Network.clust='hom'; % heterogeneous EE cluster size
Network.clust_std=0.01; % het cluster size: cluster size is sampled from a gaussian with SD=Ncluster*Network.std
Network.clustEI='hom'; % EI clusters: homogeneous='hom', heterogeneous='het'
Network.clustIE='hom'; % IE clusters: homogeneous='hom', heterogeneous='het'
Network.clustII='hom'; % II clusters: homogeneous='hom', heterogeneous='het'
% Network.clust_syn='';
N=1000; % network size
N_e = N/2; % exc neurons
N_i = N/2; % inh neurons
Scale=(5000/N)^(1/2);
% % global spontaneous firing rates (neeed to fix thresholds)
% ni_e = 2; % 3.6;   %6.6   % AB97: 3.0
% ni_i = 5; % 5.2;   %8.2   % AB97: 4.2
%------------------
% TIME CONSTANTS
%------------------
% tau_arp = .005;  % refractory period
tau_i = .03; %+(EIC-1)*0.002;     % inh membrane time
tau_e = .02;	 % exc membrane time
% tausyn_e=0.0046;  % exc synaptic time
% tausyn_i=0.0035;  % inh synaptic time
% tausyn_e=0.0044;
% tausyn_i=0.0035;
tausyn_e=0.011;
tausyn_i=0.05;
%------------------
% SYNAPTIC WEIGHTS
%------------------
% syn weights are drawn from a gaussian distribution with std delta and
% mean Jab
% Jee = Scale*0.02;
% Jii = 1.*Scale*0.12; %     Jii = Scale*0.06;
% Jie = 2*Scale*0.010; %     Jei 3 Scale*0.045;
% Jei = 3.*Scale*0.02;
  Jee = Scale*0.01;
%   Jee = 0; % mean Jee weight
% Jii = (5/2)^(1/2)*Scale*0.06;
 Jii = Scale*0.02;
Jie = Scale*0.02;
% Jei = (5/2)^(1/2)*Scale*0.045;
Jei = Scale*0.02;
%------------------
% CLUSTER PARAMETERS
%------------------
delta = 0.01;
% delta=0.;% SD of synaptic weights distribution, eg: ~Jee*(1+delta*randn(N_e))
Network.delta = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
Network.deltaEI = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
Network.deltaIE = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
Network.deltaII = delta; % SD of synaptic weights: Jee*(1+delta*randn(N_e)) Larger delta -> clusters look more async in the high state
Jplus = 1; % EE intra-cluster potentiation factor 10
Network.factorEI = 4; %+ CountDrac/2; % EI intra-cluster potentiation factor 4
Network.factorIE = 4; %+ CountDrac/2; % IE intra-cluster potentiation factor 8
% if StimulusSelection==2
%     Network.factorEI = 40; %+ CountDrac/2; % EI intra-cluster potentiation factor 4
% Network.factorIE = 40; %+ C
% end

Network.factorII = 1.6; %-(EIC-1); 
% switch EIC
% case 1
% Network.factorEI = 4; %+ CountDrac/2; % EI intra-cluster potentiation factor 4
% Network.factorIE = 4; %+ CountDrac/2; % IE intra-cluster potentiation factor 8
% Network.factorII = -1.6; % II intra-cluster potentiation factor 8
% Jii = Scale*0.01;
% Jie = Scale*0.01;
% Jei = Scale*0.01;
% case 2 
%     Network.factorEI = 8; %+ CountDrac/2; % EI intra-cluster potentiation factor 4
% Network.factorIE = 8; %+ CountDrac/2; % IE intra-cluster potentiation factor 8
% Network.factorII = -1.6; % II intra-cluster potentiation factor 8
% Jii = Scale*0.01;
% Jie = Scale*0.01;
% Jei = Scale*0.01;
% case 3 
%     Network.factorEI = 8; %+ CountDrac/2; % EI intra-cluster potentiation factor 4
% Network.factorIE = 8; %+ CountDrac/2; % IE intra-cluster potentiation factor 8
% Network.factorII = -1.6; % II intra-cluster potentiation factor 8
% Jii = Scale*0.01;
% Jie = Scale*0.005;
% Jei = Scale*0.005;
% case 4 
%     Network.factorEI = 16; %+ CountDrac/2; % EI intra-cluster potentiation factor 4
% Network.factorIE = 16; %+ CountDrac/2; % IE intra-cluster potentiation factor 8
% Network.factorII = -1.6; % II intra-cluster potentiation factor 8
% Jii = Scale*0.01;
% Jie = Scale*0.005;
% Jei = Scale*0.005;
% end
%
%Original: theta_e=1.8;
%theta_e=2.2;

% theta_i=1.65 + 4*0.02;
theta_e = 1.92;
theta_i=2.5; %changes between breath spontaneous firing
    
%     
tau_arp = .002;

% EIC = 4;
% Network.factorII = -2 + (EIC-1)*0.4;

%     Network.factorEI= 2 + 0.1*round(EIC/2);
% else
%     Network.factorEI= 2 - 0.1*EIC/2;
% end

%% Condition-5


%% Condition-4

% switch CountDrac
%     case 1
%         tausyn_e = tausyn_e + 0.02*CountDrac2*tausyn_e;
%     case 2 
%         tausyn_e = tausyn_e - 0.02*CountDrac2*tausyn_e;
%     case 3
%          tausyn_i = tausyn_i + 0.02*CountDrac2*tausyn_i;
%     case 4 
%         tausyn_i = tausyn_i - 0.02*CountDrac2*tausyn_i;
% end

%% Condition-3

% switch CountDrac
%     case 1
%         Jie = Jie + 0.02*CountDrac2*Jie;
%     case 2 
%         Jie = Jie - 0.02*CountDrac2*Jie;
%     case 3
%          Jei = Jei + 0.02*CountDrac2*Jei;
%     case 4 
%         Jei = Jei - 0.02*CountDrac2*Jei;
% end

%% Condition-2
% No sign change observed
% switch CountDrac
%     case 1
%         Network.factorEI=Network.factorEI + 0.2*CountDrac2;
%     case 2 
%         Network.factorEI=Network.factorEI - 0.2*CountDrac2;
% 
%     case 3
%         Network.factorIE=Network.factorIE + 0.2*CountDrac2;
%     case 4 
%         Network.factorIE=Network.factorIE - 0.2*CountDrac2;
% end

%% Condition-1
%%start with 1.8; 1.7
% switch CountDrac
%     case 1
%         theta_e=theta_e + 0.03*CountDrac2*theta_e;
%     case 2 
%         theta_e=theta_e - 0.03*CountDrac2*theta_e;
% 
%     case 3
%         theta_i=theta_i + 0.03*CountDrac2*theta_i;
%     case 4 
%         theta_i=theta_i - 0.03*CountDrac2*theta_i;
% end


%------------------
% CONNECTIVITY PARAMETERS
%------------------

Cee = N_e*0.2; % # presynaptic neurons
Cie = N_e*0.2; %
Cii = N_i*0.2; %
Cei = N_i*0.2; %


bgr=0.0; % fraction of background excitatory neurons (unclustered)
Network.bgrI=0.0; % fraction of background neurons
Ncluster = 10; % average # neurons per cluster
p = round(N_e*(1-bgr)/Ncluster); % # of clusters
% p = N_i*4;
f = (1-bgr)/p;
% f = Ncluster/N_e;
Network.fI = (1-Network.bgrI)/p;       % 0.09
% gam = 0.5;% % parameter related to inter-cluster depression (see function aux.SynWeights)
rho = 2.75;
%------------------
% THRESHOLDS
%------------------
%
% theta_e=164.938*Jee; % exc threshold potential
% theta_i=169.368*Jee; % inh threshold potentials
% 
% theta_e=2.7;
% theta_i=2.6;

% % theta_e=1; % exc threshold potential
% theta_i=1; % inh threshold potentials
% reset potentials
He = 0;%
Hi = 0;%



%------------------
% EXTERNAL BIAS
%------------------
% external input parameters, eg: external current given by mu_e_ext=Cext*Jee_ext*ni_ext
Cext = (N_e)*0.2; % # presynaptic external neurons
Jie_ext=0.8*Scale*0.0915;% external input synaptic strengths
Jee_ext=0.8*Scale*0.1027; %
% EXTERNAL CURRENT
% default external currents
ni_ext = 5; % 7;
mu_E0=Cext*Jee_ext*ni_ext;
mu_I0=Cext*Jie_ext*ni_ext;
% random ext current, different in each trial
% Mu=[mu_E0*(ones(N_e,1)+(0.1/2)*(2*rand([N_e,1])-1)); ...
%     mu_I0*(ones(N_i,1)+(0.05/2)*(2*rand([N_i,1])-1))];     % bias
% Mu=[mu_E0*(ones(N_e,1)); mu_I0*(ones(N_i,1))];     % bias
Mufile = 'ExternalCurrent.mat';
Mu = load(Mufile);
Mu = Mu.Mu;

%% Connectivity
Inhfile = 'Inhfix.mat';
GoOdor = 'GoOdor.mat';
NoGoOdor = 'NoGoOdor.mat';

OverlapGain = 4;
GlomOverlap = 0;
NcUnits=round(f*N_e); 
Q=p;
popsize=repmat(NcUnits,Q,1);
cusumNcE=[0 cumsum(popsize)'];

if StimulusSelection == 1 && imjk == 1 && CountDrac2 == 1 && EIC == 1
%     disp('flag');
% for k =1:100
    for i = 1:length(cusumNcE)-1
        temp = [i i i i i i i i];
        perm(i,:) = temp;
    end
    
    Inh = zeros(50,8);
    for i =1:size(perm,2)
        Inh(:,i) = randperm(numel(perm(:,i)));
    end
    %% Making connections between DAT and Glom
    for i =1:size(perm,1)
        tb = tabulate(Inh(i,:));
        if any(tb(:,2)>1)
            xm = find(tb(:,2)==2);
            for k21 = 1:length(xm)
                xo = find(Inh(i,:) == xm(k21));
                Inh(i,xo(1)) = 0;
            end
        end
    end
    
    %% Choosing 15% of total glom as Go Odor
    tempOdor = ones(1,50);
    a=randperm(numel(tempOdor));
    tempOdor(a(1:round(0.9*numel(tempOdor))))=0;
    Odor1 = tempOdor;
    
    %% Finding DAT cells that have no connections with any Go Odor glom. 
    alldat = 1:50;
    for i =1:50
        if Odor1(i) ==1
            [x,~] = (find(Inh == i));
            alldat = setdiff(alldat,x,'stable');
        end
    end
%     UsedDat = setdiff(1:50,alldat,'stable');

%% Glom of unused DAT, can stil be connected to other Go Odor DATs depending on connections with unused DAT found. More means less connections with GoDATs
    tempOdor = [];
    for i = alldat
        tempOdor = [tempOdor Inh(i,:)];
    end
    GlomLeft = tabulate(tempOdor);
    
    %% Choosing NoGo-Odor
    UsableGlom = [];
    count = 1;
    for i = 1:length(GlomLeft)
        if GlomLeft(i,2) >= OverlapGain && GlomLeft(i,1) ~= 0
            UsableGlom(count) = GlomLeft(i,1);
            count = count+1;
        end
    end
    a=randperm(numel(UsableGlom));
    UsableGlom = UsableGlom(a);
    tempOdor = zeros(1,50);
    if length(UsableGlom)>nnz(Odor1)
        tempOdor(UsableGlom(1:nnz(Odor1))) = 1;
    else
        tempOdor(UsableGlom(1:length(UsableGlom))) = 1;
    end
    Odor2 = tempOdor;
    
    %% Deleting extra glom
    ngdelete =[];
    ngdelete = find(Odor2==1);
    ngdelete = ngdelete(randperm(numel(ngdelete)));
%     if nnz(Odor2)>10
    if GlomOverlap<=length(ngdelete)
        Odor2(ngdelete(1:GlomOverlap)) = 0;
    else
        Odor2(ngdelete) = 0;
    end
% %         break;
% %     elseif nnz(Odor2)==11
% %         Odor2(ngdelete(1:(EIC-2))) = 0;
% %         break;
% %     elseif nnz(Odor2)==11
% %         Odor2(ngdelete(1:(EIC-2))) = 0;
% %     break;
    %% Overlapping some Glom
    if GlomOverlap>0
    gngGlom = [];
    gngGlom = find(Odor1==1);
    gngGlom = gngGlom(randperm(numel(gngGlom)));
    Odor2(gngGlom(1:(GlomOverlap))) = 1;
    end
     %% making DAT number same
    alldatgo = [];
    go = find(Odor1==1);
    for ii =go
       [x,~] = (find(Inh == ii));
       alldatgo = [alldatgo; x];
    end
    
    alldatnogo = [];
    nogo = find(Odor2==1);
    for ii =nogo
            [x,~] = (find(Inh == ii));
            alldatnogo = [alldatnogo; x];
    end    
          Inh1= Inh;
    Inh = [Inh zeros(50,2)]; 
    a2 = tabulate(alldatnogo);
    a1 = tabulate(alldatgo);
    [xa1,ya1] = max(a1(:,2));
    [xa2,ya2] = max(a2(:,2));

    x12 =setdiff(go,Inh(ya1,:));
    if xa1 == 2
        Inh(ya1,9:10) = (x12(1:2));
    elseif xa1 == 3
        Inh(ya1,9) = (x12(1));
    end
    
    x12 =setdiff(nogo,Inh(ya1,:));
    if xa2 == 2
        Inh(ya2,9:10) = (x12(1:2));
    elseif xa2 == 3
        Inh(ya2,9) = (x12(1));
    end
    
 alldatgo = [];
    go = find(Odor1==1);
    for ii =go
       [x,~] = (find(Inh == ii));
       alldatgo = [alldatgo; x];
    end
    
    alldatnogo = [];
    nogo = find(Odor2==1);
    for ii =nogo
            [x,~] = (find(Inh == ii));
            alldatnogo = [alldatnogo; x];
    end    
    
    
    if length(alldatgo) < length(alldatnogo)
            removeDAT = length(alldatnogo)- length(alldatgo);
            for iii = 1:removeDAT
            aa = tabulate(alldatnogo);
            aa = find(aa(:,2) == 1);
            aaa=randperm(numel(aa(:)));
            aa = aa(aaa);
            for iv = 1:length(nogo)
                [~,y] = find(nogo(iv) == Inh(aa(1),:));
                Inh(aa(1),y) = 0;
            end
            alldatnogo((alldatnogo == aa(1))) = [];
            end
    end
    if length(alldatgo) > length(alldatnogo)
            removeDAT = length(alldatgo)- length(alldatnogo);
            for iii = 1:removeDAT
            aa = tabulate(alldatgo);
            aa = find(aa(:,2) == 1);
            aaa=randperm(numel(aa(:)));
            aa = aa(aaa);
            for iv = 1:length(go)
                [~,y] = find(go(iv) == Inh(aa(1),:));
                Inh(aa(1),y) = 0;
            end
            alldatgo((alldatgo == aa(1))) = [];
            end
    end

 %% increasing dat numbers   
%     while length(alldatgo)<42
%        [x0,y0] = find(Inh == 0);
%        x12 = randperm(numel(go));
%        for x11 = 1:length(x0) 
%         if nnz(ismember(go,Inh(x0(x11),:)))
%             Inh(x0(x11),y0(x11)) = go(x12(1));
%             break;
%         end
%        end
%        alldatgo = [];
%        for ii =go
%        [x,~] = (find(Inh == ii));
%        alldatgo = [alldatgo; x];
%     end
%     end
%      while length(alldatnogo)<42
%        [x0,y0] = find(Inh == 0);
%        x12 = randperm(numel(nogo));
%        for x11 = 1:length(x0) 
%         if nnz(ismember(nogo,Inh(x0(x11),:)))
% %             x12 = randperm(numel(nogo));
%             Inh(x0(x11),y0(x11)) = nogo(x12(1));
%             break;
%         end
%        end
%        alldatnogo = [];
%        for ii =nogo
%        [x,~] = (find(Inh == ii));
%        alldatnogo = [alldatnogo; x];
%     end
%      end   
%         a2 = tabulate(alldatnogo);
%     a1 = tabulate(alldatgo);
%     a1 = tabulate(a1(:,2));
%     a2 = tabulate(a2(:,2));
%     connection(k,:) = [a1(end,1)' a2(end,1)']
% end
    
    formatSpec = 'Simulations/k%d';
    
    file_name = sprintf(formatSpec,CountDrac);
    
    parsave(fullfile(pwd,file_name,Inhfile),Inh);
    parsave(fullfile(pwd,file_name,GoOdor),Odor1);
    parsave(fullfile(pwd,file_name,NoGoOdor),Odor2);

 end
% elseif StimulusSelection ~= 1 && imjk ~= 1 && CountDrac2 ~= 1 && EIC ~= 1
% else
    formatSpec = 'Simulations/k%d';
    file_name = sprintf(formatSpec,CountDrac);
    Inh = par_load(fullfile(pwd,file_name,Inhfile));
    Inh = Inh.x;
    Odor1 = par_load(fullfile(pwd,file_name,GoOdor));
    Odor1 = Odor1.x;
    Odor2 = par_load(fullfile(pwd,file_name,NoGoOdor));
    Odor2 = Odor2.x;


    for i =1:length(Inh)
        if Odor1(i) == 1
             [xO1,~] = (find(Inh == i));
        end
        if Odor2(i) == 1
             [xO2,~] = (find(Inh == i));
        end
    end
    OverlappedDAT = length(intersect(xO1,xO2));
% end
% if StimulusSelection == 1 && imjk ~= 1 && CountDrac2 ==1
%     formatSpec = 'Simulations/k%d';
%     file_name = sprintf(formatSpec,CountDrac);
%     Inh = par_load(fullfile(pwd, file_name,Inhfile));
%     Inh = Inh.x;
%     Odor1 = par_load(fullfile(pwd, file_name,GoOdor));
%     Odor1 = Odor1.x;
%     Odor2 = par_load(fullfile(pwd, file_name,NoGoOdor));
%     Odor2 = Odor2.x;
% end
%     
% if StimulusSelection > 1 && CountDrac2 ==1
%     formatSpec = 'Simulations/k%d';
%     file_name = sprintf(formatSpec,CountDrac);
%     Inh = par_load(fullfile(pwd, file_name,Inhfile));
%     Inh = Inh.x; 
%     Odor1 = par_load(fullfile(pwd, file_name,GoOdor));
%     Odor1 = Odor1.x;
%     Odor2 = par_load(fullfile(pwd, file_name,NoGoOdor));
%     Odor2 = Odor2.x;
% end
% 
% if CountDrac2 ==2
%     formatSpec = 'Simulations/k%d';
%     file_name = sprintf(formatSpec,CountDrac);
%     Inh = par_load(fullfile(pwd, file_name,Inhfile));
%     Inh = Inh.x; 
%     Odor1 = par_load(fullfile(pwd, file_name,GoOdor));
%     Odor1 = Odor1.x;
%     Odor2 = par_load(fullfile(pwd, file_name,NoGoOdor));
%     Odor2 = Odor2.x;
% end
% if EIC == 1 || EIC == 2 || EIC == 3 || EIC == 4
%     bgain = 1.5; 
% elseif EIC == 5 || EIC == 6 || EIC == 7 || EIC == 8  
    bgain = 2; 
% elseif EIC == 9 || EIC == 10 || EIC == 11 || EIC == 12      
%      bgain = 2.5; 
% elseif EIC == 13 || EIC == 14 || EIC == 15 || EIC == 16     
%     bgain = 3; 
% elseif EIC == 17 || EIC == 18 || EIC == 19 || EIC == 20  
%     bgain = 3.5;
% end
%     
    
rr = randi(10,1,1);
Tseq = -1:0.0001:1;
p2 = @(t)(sin(rr-5*t));
rr2 = 0;
while rr2<0.5 && rr2>-0.5
    rr2 = (-10 + (20)*rand(1,1))/10;
end

% for t =1:9999
%      temp2(t) = p(rr2*3.14*Tseq(t));
% end
temp2 = 0;
for t = 10000:14999
     temp2(t) = p2(rr2*3*bgain*3.14*Tseq(t));
end
OdorStart = 0.0001*find(temp2(10000:12000)>0.999,1);
% for t = 15000:20000
%      temp2(t) = p(rr2*3.14*Tseq(t));
% end

%  plot(temp2)


%----------------
% DEFAULT STIMULI
%----------------
% STIMULUS
Stimulus.option='on'; % 'off'=no stimulus
Stimulus.input='Const'; % constant external current
scnt=0;
% TASTE (specific stimulus)
scnt=scnt+1;
feat(scnt).name='US'; % unconditioned stimulus (taste)
% feat(scnt).interval=[0 Sim.t_End]; % stimulus interval
OdorStart = 0;
feat(scnt).interval=[OdorStart 0.5];
%     gain = 0.00000000001;   

    gain = 0.1;   

feat(scnt).gain=gain;
feat(scnt).profile=@(t)t*(0.5-t); % time course of stimulus, eg a linear ramp
feat(scnt).selectivity='exc';
% if StimulusSelection>200 && StimulusSelection<401
%     feat(scnt).selective=rand(1,50)<0.4; % US selective clusters
% else
    if rem(StimulusSelection,2) == 1
        feat(scnt).selective = Odor1;
    else
        feat(scnt).selective = Odor2;     
    end
% end
feat(scnt).connectivity=1;       
    
% ANTICIPATORY CUE
scnt=scnt+1;
feat(scnt).name='CSgauss'; % conditioned stimulus (cue) with "quenched" noise
% dur = 0.05+ 0.1*(CountDrac-1); % 1-20
% start = -0.5 + 0.05*(CountDrac2-1); % 1-11
start = -0.5;
dur = 0.6;
zero = start +dur;
feat(scnt).interval=[start zero]; % stimulus interval
gain=3; % SD of quenched noise across neurons
feat(scnt).gain=gain;
tau_cue=[0.5,1]; % rise and decay time of double exp cue time course
formatSpec = '@(t)(%0.03g-t)';
str = sprintf(formatSpec,start);
feat(scnt).profile=str2func(str);
feat(scnt).selectivity='mixed'; % targets exc neurons only
feat(scnt).selective = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
feat(scnt).connectivity=1; % fraction of neurons targeted within each cluster

scnt=scnt+1;
feat(scnt).name='hab'; % conditioned stimulus (cue) with "quenched" noise
% dur = 0.05+ 0.1*(CountDrac-1); % 1-20
% start = -0.5 + 0.05*(CountDrac2-1); % 1-11
 start2 = 0.1;
if rem(StimulusSelection,2) == 1
    dur2 = 0.4;
else
    dur2 = 0.3;
end

%  dur2 = 0.4;
zero2 = start2+dur2;
feat(scnt).interval=[start2 zero2]; % stimulus interval
gain=3; % SD of quenched noise across neurons
feat(scnt).gain=gain;
tau_cue=[0.5,1]; % rise and decay time of double exp cue time course
formatSpec = '@(t)(%0.03g-t)*(%0.03g-t)';
str = sprintf(formatSpec,start2,zero2);
feat(scnt).profile=str2func(str);
feat(scnt).selectivity='mixed'; % targets exc neurons only
% feat(scnt).selective=ones(1,p); % CS targets all clusters
feat(scnt).selective = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
feat(scnt).connectivity=1; % fraction of neurons targeted within each cluster

scnt=scnt+1;
feat(scnt).name='breath'; % conditioned stimulus (cue) with "quenched" noise
% dur = 0.05+ 0.1*(CountDrac-1); % 1-20
% start = -0.5 + 0.05*(CountDrac2-1); % 1-11
start2 = -1;
zero2 = 1;
feat(scnt).interval=[start2 zero2]; % stimulus interval
gain= 0.1; % SD of quenched noise across neurons
feat(scnt).gain=gain;
formatSpec = '@(t)(sin(%0.03g-5*t))';
% rr = randi(10,1,1);
str = sprintf(formatSpec,rr);
feat(scnt).profile=str2func(str);
feat(scnt).selectivity='exc'; % targets exc neurons only
feat(scnt).selective = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
feat(scnt).connectivity=1; % fra
Stimulus.feat=feat;


%------------------------------------------------------------------------
% PLOT PARAMETERS: ------------------------------------------------
%------------------------------------------------------------------------/
% PLOTS
Sim.Plotf=0;
Sim.plot_length=Sim.t_End-Sim.t_Start; % length of plot intervals
% indices of ensemble units to store
exc=randperm(N_e);
inh=N_e+1:N_e+N_i;
inh=inh(randperm(numel(inh)));
Sim.ind_p=[exc(1) inh(1)]; % choosing neuron index for membrane potential plot (one E and one I)
Sim.weights_save='off'; % save weight matrix: 'Yes'
extra='';
%
%
% paramsfile='params2.mat';

events={'US'};
if strcmp(condition,'UT')
    events={'US'}; % US only
elseif strcmp(condition,'ET')
    events={'CSgauss','US','hab','breath'}; % CS + US
end


params2 = struct('ni_ext',ni_ext,'tau_arp',tau_arp,'tau_i',tau_i,'tau_e',tau_e,'theta_e',theta_e,...
    'theta_i',theta_i,'delta',delta,'f',f,'Jee',Jee,'Jii',Jii,'Jie',Jie,'Jei',Jei,'Jee_ext',Jee_ext,'Jie_ext',Jie_ext,...
    'Jplus',Jplus,'He',He,'Hi',Hi,'N_e',N_e,...
    'N_i',N_i,'Cee',Cee,'Cie',Cie,'Cei',Cei,'Cii',Cii,'Cext',Cext,'p',p,'Sim',Sim,'Network',Network,'Stimulus',Stimulus,...
    'tausyn_e',tausyn_e,'tausyn_i',tausyn_i,'extra',extra,'Mu',Mu,'events',{events});



