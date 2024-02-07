function [stimulus_save, params]=fun_stim(params,Inh,i,j,m,k)
CountDrac = k;
Gaincontrol  = j;
Gain2 =m;
% unpack vars
learn = i;
p=params.p;
Stimulus=params.Stimulus;
events=params.events;
Sim=params.Sim;
BiasDistribution =1;
% for each event, create external current and  stim properties in .Ext
% and add clusters selective to the event .Stimulus.feat(n).StimClust
stimulus_save=struct('Ext',[],'Stimulus',[]);

% select stimuli
stimuli=events;
temp_Stimulus=struct('option',Stimulus.option,'input',Stimulus.input);
indfeat=zeros(1,numel(stimuli));
for ev=1:numel(stimuli)
    %     fprintf('Stimulus %s',stimuli{ev});
    % match current stimuli to features in Stimulus
    indfeat(ev)=find(cell2mat(arrayfun(@(x)strcmp(stimuli{ev},x(:).name),...
        Stimulus.feat(:),'uniformoutput',false)));
end
fprintf('\n');
if ~isempty(indfeat)
    temp_Stimulus.feat(1:numel(indfeat))=Stimulus.feat(indfeat);
    for n=1:numel(indfeat)
        sclust=[];
        if ~isempty(temp_Stimulus.feat(n).selective)
            sclust=find(temp_Stimulus.feat(n).selective(1,:));
        end
        temp_Stimulus.feat(n).StimClust=sclust;
    end
end
Stimulus=temp_Stimulus;
Ext=struct('Mu',[]);

% LOAD PARAMETERS
fieldNames={'Sim','Network','p','popsize','clustermatrix','N_e','N_i','Cext','Jee_ext','Jie_ext','ni_ext','tau_e','tau_i','fieldNames'};
aux.v2struct(params,fieldNames);
Q=p; % number of clusters
cusumNcE=[0 cumsum(popsize)'];
Tseq=Sim.t_Start:Sim.dt_step:Sim.t_End;
% Inhfile = 'Inhfix.mat';
% Inh = load(Inhfile);
% Inh = Inh.Inh;
temp_Inh = Inh;
DATGofile = 'DATGofix.mat';
DATNoGofile = 'DATNoGofix.mat';







if strcmp(Stimulus.option,'on')
    feat=Stimulus.feat;
    nstim=numel(feat); % number of stimuli in current trials
    stim=repmat(struct('profile',[],'ind',[],'interval',[]),1,nstim);
    temp_ind=repmat(struct('ind',[]),1,nstim); % stores indices for mixed cue (see below)
    % selective neurons
        
    for n=1:nstim
        StimClust=Stimulus.feat(n).StimClust; % clusters activated by current stimulus
        % units selective to stimulus
        ind=[];
        switch feat(n).name
            case 'US'
                for c=StimClust
                    pop_ind=find(clustermatrix(:,c));
                    for k=1:numel(pop_ind)
                        ind=[ind cusumNcE(pop_ind(k))+1:cusumNcE(pop_ind(k)+1)]; % stim selective units
                    end % To get all the possible active neurons
                end
                DATGo = [];
                DATNoGo = [];
                if rem(learn,2) == 1
                    for dtg = 1:length(ind)
                        [xdtg,~] = find(Inh == ceil(ind(dtg)/10));
                        DATGo= [DATGo; xdtg];
                    end
                    DATGo = unique(DATGo);
                    DATGoind = [];
                    for dtg = 1:length(DATGo)
                        DATGoind =[DATGoind cusumNcE(DATGo(dtg))+501:cusumNcE(DATGo(dtg)+1) + 500];
                    end
                   formatSpec = 'Simulations/k%d';
                    file_name = sprintf(formatSpec,CountDrac);
                    parsave(fullfile(pwd,file_name,DATGofile),DATGoind);
    
                else
                    formatSpec = 'Simulations/k%d';
                    file_name = sprintf(formatSpec,CountDrac);
                    DATGoind = par_load(fullfile(pwd, file_name,DATGofile));
                    DATGoind = DATGoind.x;
                    
                    for dtg = 1:length(ind)
                        [xdtg,~] = find(Inh == ceil(ind(dtg)/10));
                        DATNoGo= [DATNoGo; xdtg];
                    end
                    DATNoGo = unique(DATNoGo);
                    DATNoGoind = [];
                    for dtg = 1:length(DATNoGo)
                        DATNoGoind =[DATNoGoind cusumNcE(DATNoGo(dtg))+501:cusumNcE(DATNoGo(dtg)+1) + 500];
                    end
                   formatSpec = 'Simulations/k%d';
                    file_name = sprintf(formatSpec,CountDrac);
                    parsave(fullfile(pwd,file_name,DATNoGofile),DATNoGoind);
                end
                allDat = 501:1000; 
                sparseDAT = allDat;
                sparseDAT2 = sparseDAT;
                a1=randperm(numel(sparseDAT));
                sparseDAT(a1(1:round(0.1*(2)*numel(sparseDAT))))=[];
                sparseDAT2(a1(1:round(0.1*(0.0)*numel(sparseDAT))))=[];

                remainingDat = setdiff(allDat,DATGoind);
                sparseDATGoind = DATGoind;
                sparseremainingDat = remainingDat;
                a1=randperm(numel(DATGoind));
                a2 = randperm(numel(remainingDat));
                sparseDATGoind(a1(1:round(0.1*(5-1)*numel(sparseDATGoind))))=[];
                sparseremainingDat(a2(1:round((1-0.1*(5-1))*numel(sparseremainingDat))))=[];
                OverlappedDATGoInd = [sparseDATGoind sparseremainingDat];
    
        end
    end
    
    
    
    
    
    
    for n=nstim:-1:1
        Inh = temp_Inh;
        % stimulus interval
        interv=feat(n).interval;
        Gain=feat(n).gain;
        Profile=feat(n).profile;
        start = interv(1);
        stop = interv(2);
        Profile=@(t)Profile(t-interv(1));
        MaxThInput = 0;
        if ~isempty(strfind(feat(n).name,'gauss' ))
            if rem(Gaincontrol,2) == 1
                Gain = 0.000000001; % (learn-1)*0.5 + with gaussian stim set profile to peak at 1, then multiply each profile by gaussian with SD feat(n).gain for each neuron in feat(n).gauss
            else
                Gain = 0.25;
            end
            Profile=feat(n).profile;
            
        elseif ~isempty(strfind(feat(n).name,'hab'))
%             if learn<201
                if rem(learn,2) == 1
                    if rem(Gaincontrol,2) == 1
                        Gain = 1.0; 
                    else
                        Gain= 1.0;
                     end
%                 else
%                     Gain= 0.5;
                
            else
                Gain = 0.5;
            end
            Profile=feat(n).profile;
            
            
        elseif ~isempty(strfind(feat(n).name,'breath'))
                 Gain = 0.1;         
            Profile=feat(n).profile;
        end
        
        temp = [];
        
        for k = 1: length(Tseq)
            if Tseq(k)>interv(1) && Tseq(k)<interv(2)
                temp = Profile(Tseq(k));
                if temp<0
                    temp = abs(temp);
                end
                if temp> MaxThInput
                    MaxThInput= temp;
                    %save = k
                end
            end
        end
    
        Profile=@(t)Gain*Profile(t)/MaxThInput;
        stim(n).profile=@(t)Profile(t); % fraction increase above baseline
        % selective neurons
        StimClust=Stimulus.feat(n).StimClust; % clusters activated by current stimulus
        % units selective to stimulus
        ind=[]; % indices of stim sel units
        switch feat(n).name
            case 'US'
                for c=StimClust
                    pop_ind=find(clustermatrix(:,c));
                    for k=1:numel(pop_ind)
                        ind=[ind cusumNcE(pop_ind(k))+1:cusumNcE(pop_ind(k)+1)]; % stim selective units
                    end % To get all the possible active neurons
                end
               
                
                
            case 'CSgauss'
                %                 Inh = temp_Inh;
                %                 for c=StimClust
                %                     for nrow =1: length(Inh(:,1))
                %                         if (find(Inh(nrow,:) == c))
                %                             if sum(Inh(nrow,:)) ~= 0
                %                                 ind =[ind cusumNcE(nrow)+501:cusumNcE(nrow+1) + 500];
                %                                 Inh(nrow,:) = zeros(1,size(Inh,2));
                %                             end
                %                         end
                %                     end
                %                 end
                
%                 ind = OverlappedDATGoInd;
%                 ind = allDat;
                ind = sparseDAT;
                  
            case 'hab'
%                 Inh = temp_Inh;
%                 for c=StimClust
%                     for nrow =1: length(Inh(:,1))
%                         if (find(Inh(nrow,:) == c))
%                             if sum(Inh(nrow,:)) ~= 0
%                                 ind =[ind cusumNcE(nrow)+501:cusumNcE(nrow+1) + 500];
%                                 Inh(nrow,:) = zeros(1,size(Inh,2));
%                             end
%                         end
%                     end
%                 end
                ind = sparseDAT2;

%                 ind = OverlappedDATGoInd;
%                 if rem(learn,2) == 1
%                     ind = DATGoind;
%                 else
%                     ind = DATNoGoind;
%                 end
            otherwise
                ind=1:N_e;
        end
        % sparsify
%         a=randperm(numel(ind));
        temp_ind(n).ind=ind;
%         ind=ind(a(1:round(feat(n).connectivity*numel(ind))));
        % gaussian stimulus, draw from randn
        if ~isempty(strfind(feat(n).name,'gauss'))
%             switch Gain2
%                 case 1
%                     aa =  randn(numel(allDat),1);
%                     aa(aa<0) = 0;
%                     %             aa = (randn(numel(allDat),1));
%                     aaa = randperm(numel(aa(:)));
%                     aa(DATGoind-500) = aa(DATGoind-500) + abs(randn(numel(DATGoind),1));
%                 case 2
%                     aa =  1+randn(numel(allDat),1);
%                     %              aa(aa<0) = 0;
%                     %             aa = (randn(numel(allDat),1));
%                     aaa = randperm(numel(aa(:)));
%                     aa(DATGoind-500) = aa(DATGoind-500) + abs(randn(numel(DATGoind),1));
%                 case 3
%                     aa =  abs(randn(numel(allDat),1));
%                     %              aa(aa<0) = 0;
%                     %             aa = (randn(numel(allDat),1));
%                     aaa = randperm(numel(aa(:)));
%                     aa(DATGoind-500) = aa(DATGoind-500) + abs(randn(numel(DATGoind),1));
%                 case 4
%                     aa =  abs(randn(numel(allDat),1));
%                     %              aa(aa<0) = 0;
%                     %             aa = (randn(numel(allDat),1));
%                     aaa = randperm(numel(aa(:)));
%                     aa(DATGoind-500) = aa(DATGoind-500) + 1+abs(randn(numel(DATGoind),1));
%                 case 5
%                     aa =  randn(numel(allDat),1);
%                     %              aa(aa<0) = 0;
%                     %             aa = (randn(numel(allDat),1));
%                     aaa = randperm(numel(aa(:)));
%                     aa(DATGoind-500) = aa(DATGoind-500) + 1 + abs(randn(numel(DATGoind),1));
%             end

%             if Gain2>1
%                 aa(DATGoind-500) = aa(DATGoind-500) + 2 + randn(numel(DATGoind),1);
%                 
%             else
%                 aa(DATGoind-500) = aa(DATGoind-500);
%             end
%              if Gain2>1
                
%              else
%                 aa(DATGoind-500) = aa(DATGoind-500);
%             end
              aa = abs(randn(numel(sparseDAT),1));
              aaa=randperm(numel(aa(:)));
              aa(aaa(1:ceil(0.0*numel(aa(:)))))=0;

              stim(n).gauss = feat(n).gain*aa; % Only Positve
            
            
            
            
%                         stim(n).gauss=feat(n).gain*(randn(numel(ind),1));
%             switch Gain2
%                 case 1
%                     Gainfile = 'ExternalGainLow1.mat';
%                     aa = load(Gainfile);
%                     aa = aa.aa;
%                     stim(n).gauss=feat(n).gain*aa;
%                 case 2
%                     Gainfile = 'ExternalGainHigh1.mat';
%                     aa = load(Gainfile);
%                     aa = aa.aa;
%                     stim(n).gauss=feat(n).gain*aa;
%                 case 3
%                     Gainfile = 'ExternalGainHigh2.mat';
%                     aa = load(Gainfile);
%                     aa = aa.aa;
%                     stim(n).gauss=feat(n).gain*aa;
%             end
        elseif ~isempty(strfind(feat(n).name,'hab')) 
%                 if rem(learn,2) == 1
%                    aa = abs(randn(numel(DATGoind),1));
% 
%                 else
%                     aa = abs(randn(numel(DATNoGoind),1));
% 
%                 end
%             aa = abs(randn(numel(DATGoind),1));
            aa = abs(randn(numel(sparseDAT2),1));
            
            aaa=randperm(numel(aa(:)));
            
            aa(aaa(1:ceil(0.0*numel(aa(:)))))=0;

            stim(n).gauss=feat(n).gain*aa;

        elseif ~isempty(strfind(feat(n).name,'breath')) 

            aa = abs(randn(numel(ind),1));
%              aa(aa<0) = 0.5+aa;
 aaa=randperm(numel(aa(:)));
            
            aa(aaa(1:ceil(0.0*numel(aa(:)))))=0;
            stim(n).gauss=feat(n).gain*aa; % Only 
        end
        %
        %         ind = load('InputUnitsfix.mat');
        %         ind = ind.a;
        stim(n).ind=ind;
        stim(n).interval=interv;
        stim(n).name=feat(n).name;
        stim(n).StimClust=StimClust;
        stim(n).selectivity=feat(n).selectivity;
        
        
    end
    Ext.stim=stim;
end
Ext.Mu=params.Mu;



stimulus_save.Ext=Ext;
stimulus_save.Stimulus=temp_Stimulus;



