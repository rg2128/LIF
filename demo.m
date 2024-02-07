% Demo script for running spiking network simulations and analyses
%
% by Luca Mazzucato 2019
%
% ----------------------------------------
% Please cite:
% L. Mazzucato, G. La Camera, A. Fontanini
% Expectation-induced modulation of metastable activity underlies faster coding of sensory stimuli,
% Nat. Neuro. 22, 787-796 (2019).
% ----------------------------------------
% This script run simulations of LIF networks with excitatory (E) and inhibitory (I) spiking neurons.
% You may run 2 different network architectures:
% 1) A network with E clusters only (ClustersOption='E') [This part reproduces results from L. Mazzucato et al., 2019]
% 2) A network with E and I clusters (ClustersOption='EI')  [This part generates unpublished results (manuscript in preparation)]
%
%  maxRT = ones(25,20);
% % Iteration = 0 ;
% tic
% for h=1
% % maxRTtemp = ones(5,1);
tic
 delete(gcp('nocreate'))
 parpool(16)
AllTrace = cell(1,100,2,1,16);
AllTrace2 = cell(1,100,2,1,16);
AllTrace3 = cell(1,100,2,1,16);
AllTrace4 = cell(1,100,2,1,16);
% AllTrace = cell(1,2,2,16);
% AllTrace2 = cell(1,2,2,16);
% AllTrace3 = cell(1,2,2,16);
% AllTrace4 = cell(1,2,2,16);
%  AllTrace = cell(20,11);
% close all

parfor k =1:16
% for k =1
%     SomeTrace(imjk,i,j,m)
    SomeTrace = cell(1,100,2,1);
    SomeTrace2 = cell(1,100,2,1);
    SomeTrace3 = cell(1,100,2,1);
    SomeTrace4 = cell(1,100,2,1);
%     SomeTrace = cell(1,2,2);
%     SomeTrace2 = cell(1,2,2);
%     SomeTrace3 = cell(1,2,2);
%     SomeTrace4 = cell(1,2,2);
    %% CHOOSE CLUSTERED ARCHITECTURE
    for m =1
%         for j = 1:2      %             tempc2 = 10000;
            %             tempc = 1000;
            %             tempgauss = zeros(500,1);
            for imjk = 1
                
                for i =1:100
                    for j = 1:2
%                     tic
%                          rng('shuffle');
                    %     maxRTtemp = ones(1,20);
                    
                    %         Iteration = Iteration+1;
                    %         fprintf('Iteration: %d\n',Iteration);
                    
                    %             switch j1
                    %                 case 1
                    %                     k = 4;
                    %                     j=4;
                    %                 case 2
                    %                     k=11;
                    %                     j= 2;
                    %                 case 3
                    %                     k=11;
                    %                     j = 6;
                    %                 case 4
                    %                     k = 18;
                    %                     j = 4;
                    %                 case 5
                    %                     k = 20;
                    %                     j = 3;
                    %             end
                    
                    % close all;
                    CountDrac2 = j;
                    CountDrac = k;
                    Gaincontrol = k;
                    % CHOOSE CONDITION: EITHER ET OR UT
                    %             if rem(i,2) == 1
                    %                 condition='ET';
                    %             else
                    %                 condition='UT';
                    %             end
                    
                    condition='ET';
                    
                    % aux.create_params(condition);
                    %------------------------
                    ClustersOption='EI';%
                    % ClustersOption='E';%
                    %------------------------
                    % LOAD PARAMETERS
                    %------------------------
                    %         paramsfile='params2.mat'; % file where all network parameters are saved
                    if strcmp(ClustersOption,'E'); aux.create_params(condition);
                    elseif strcmp(ClustersOption,'EI'); [paramsfile, Inh,rr2, Odor1, Odor2, OverlappedDAT]=aux.create_params_EI(condition,CountDrac,CountDrac2,i,m,imjk);
                    end
                    
                    %% CHOOSE STIMULI (OR CREATE YOUR OWN)
                    % You have 3 built-in stimulus templates:
                    % 0) no stimuli: run the network during ongoing activity
                    % 1) 'US' stimulus onset at t=0s (stimulus is on until the end of the trial);
                    %   a linearly ramping external current delivered to 50% of E clusters (chosen randomly,
                    %   half of the neurons in each stimulus-selective clusters receive stimulus) with a slope of
                    %   gain=0.2*mu_ext, where mu_ext is the baseline external current (bias)
                    % 2) 'CSgauss' stimulus onset at t=-0.5s (stimulus is on until the end of
                    %   the trial); double exponential profile with rise and decay times
                    %   [0.5,1]s; the stimulus targets all E neurons. For each neuron, the
                    %   stimulus peak is drawn from a gaussian distribution with mean 0 and
                    %   standard deviation 0.2*mu_ext.
                    % CUSTOM OPTIONS:
                    % 1) edit your custom-made stimuli inside aux.create_params. Available options: 'US', 'CSgauss' from L. Mazzucato et al., Nat. Neuro. 2019
                    % 2) list which stimuli to run in current trials in the cell array 'stimuli'
                    % 3) if stimuli={}, run trial with spontaneous activity only, without any stimulus
                    %------------------------
                    % NOTE:
                    % This demo runs the unexpected (UT) and expected (ET) conditions from L. Mazzucato et al., 2019, depending on the following settings:
                    % % 1) with condition='UT' you simulate an unexpected trial, with a taste
                    % stimulus delivered at t=0s and not anticipatory cue (Fig. 3a right panel in the paper)
                    % 2) with condition='ET' you simulate an expected trial, with a taste stimulus
                    % delivered at t=0, preceded by an anticipatory cue at t=-0.5s (Fig. 3a left panel in the paper).
                    
                    %-----------------
                    
                    
                    %-----------------
                    
                    %-----------------
                    % SELECT STIMULUS
                    %-----------------
                    params2 = paramsfile;
                    % stimuli={}; % ongoing activity
                    % stimuli={'US'}; % stimulus evoked-activity targeting selective clusters
                    % stimuli={'CSgauss'}; % anticipatory cue speeds up network dynamics
                    stimuli={'US','CSgauss','hab','breath'}; % anticipatory cue preceeds stimulu delivery
                    %         save(paramsfile,'stimuli','-append');
                    [paramsfile(:).stimuli] = stimuli;
                    %          savedir=fullfile('Simulations'); if ~exist(savedir,'dir'); mkdir(savedir); end % setup directory for saving HMM data
                    
                    %% RUN SIMULATION
                    ntrials=1; % number of trials
                    %         file_sim=fullfile(savedir,'results2.mat');  % file where simulation results are saved
                    %---------------------------
                    % GENERATE SYNAPTIC WEIGHTS
                    %---------------------------
                    if i == 1
                        CountGo = -1;
                        CountNeut = -1;
                    end
                    
                    if rem(i,2) ==1
                        CountGo = CountGo+1;
                    else
                        CountNeut = CountNeut+1;
                    end
                    
                    
                    % J = N x N matrix of synaptic weights
                    % params = structure containing all network parameters
                    if strcmp(ClustersOption,'E'); [J, params]=aux.fun_SynWeights(paramsfile);
                    elseif strcmp(ClustersOption,'EI'); [J, params]=aux.fun_SynWeights_EI(paramsfile,i,CountGo,CountNeut,j,k,imjk, Inh, m);
                    end
                    [stimulus_save, params]=aux.fun_stim(params,Inh,i,j,m,k); % STIMULUS
                    %------------------------
                    % SIMULATION
                    %------------------------
                    %             if i>1
                    %                 Jfile = 'Jfix.mat';
                    %                 J = load(fullfile(pwd, 'Simulations', Jfile));
                    %                 J = J.J2;
                    %             end
%                     imagesc(J)
                    
                    firings=cell(1,ntrials); % cell array with all spike times in each trial
                    PlotData=cell(1,ntrials); % cell array with data for plotting
                    % parfor iTrial=1:ntrials % uncomment this line if you have a multi-core  machine with 4 or more cores
                    for iTrial=1:ntrials
                        ParamsRun=params;
                        ParamsRun.Ext=stimulus_save.Ext;
                        ParamsRun.Stimulus=stimulus_save.Stimulus;
                        ParamsRun.J=J;
                        %             fprintf('--- Start SIM ...\n');, Gstore{iTrial}, GstoreDec{iTrial}
                        %                     [firings{iTrial}, PlotData{iTrial}, J2, AllJ, allfct]=aux.fun_LIF_SIM(ParamsRun,Inh,i,m,k);
                        [firings{iTrial}, PlotData{iTrial}, StimulusFirings,J2]=aux.fun_LIF_SIM(ParamsRun,Inh,i,m,j,k,rr2, Odor1, Odor2);
                    end
                    
                    %% Synaptic update
% JIE1 = J(501:1000,1:500);
%   go = find(Odor1==1);
%   for q=1:length(go)
%         goAll(q,:) = go(q)*10-9:go(q)*10;
%   end
%     
%   nogo = find(Odor2==1);
%   for q=1:length(nogo)
%         nogoAll(q,:) = nogo(q)*10-9:nogo(q)*10;
%     end
 
                                
                                if i==1
                                JEIinitial =J;
%                                 JIE1 = J(501:1000,1:500);
%                                  Weightchange(1,:) = [sum(sum(JIE1(:,goAll))) sum(sum(JIE1(:,nogoAll)))   ];   
                                end
                                
%                                   JIE1 = J2(501:1000,1:500);
%                               
                                JEI = J2(1:500,501:1000);
%                                 Weightchange(i+1,:) = [sum(sum(JIE1(:,goAll))) sum(sum(JIE1(:,nogoAll)))   ];   

                JIE = J2(501:1000,1:500);
                JIEfile = 'JIEfix2.mat';
                JEIfile = 'JEIfix2.mat';
%              
                if Gaincontrol ==1
                    parsave(fullfile(pwd, 'Simulations/k1',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k1',JEIfile),JEI);
                elseif Gaincontrol ==2
                    parsave(fullfile(pwd, 'Simulations/k2',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k2',JEIfile),JEI);
                elseif Gaincontrol ==3
                    parsave(fullfile(pwd, 'Simulations/k3',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k3',JEIfile),JEI);
                elseif Gaincontrol ==4
                    parsave(fullfile(pwd, 'Simulations/k4',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k4',JEIfile),JEI);
                elseif Gaincontrol ==5
                    parsave(fullfile(pwd, 'Simulations/k5',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k5',JEIfile),JEI);
                elseif Gaincontrol ==6
                    parsave(fullfile(pwd, 'Simulations/k6',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k6',JEIfile),JEI);
                elseif Gaincontrol ==7
                    parsave(fullfile(pwd, 'Simulations/k7',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k7',JEIfile),JEI);
                elseif Gaincontrol ==8
                    parsave(fullfile(pwd, 'Simulations/k8',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k8',JEIfile),JEI);
                elseif Gaincontrol ==9
                    parsave(fullfile(pwd, 'Simulations/k9',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k9',JEIfile),JEI);
                elseif Gaincontrol ==10
                    parsave(fullfile(pwd, 'Simulations/k10',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k10',JEIfile),JEI);
                elseif Gaincontrol ==11
                    parsave(fullfile(pwd, 'Simulations/k11',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k11',JEIfile),JEI);
                elseif Gaincontrol ==12
                    parsave(fullfile(pwd, 'Simulations/k12',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k12',JEIfile),JEI);
                elseif Gaincontrol ==13
                    parsave(fullfile(pwd, 'Simulations/k13',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k13',JEIfile),JEI);
                elseif Gaincontrol ==14
                    parsave(fullfile(pwd, 'Simulations/k14',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k14',JEIfile),JEI);
                elseif Gaincontrol ==15
                    parsave(fullfile(pwd, 'Simulations/k15',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k15',JEIfile),JEI);
                elseif Gaincontrol ==16
                    parsave(fullfile(pwd, 'Simulations/k16',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k16',JEIfile),JEI);
                elseif Gaincontrol ==17
                    parsave(fullfile(pwd, 'Simulations/k17',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k17',JEIfile),JEI);
                elseif Gaincontrol ==18
                    parsave(fullfile(pwd, 'Simulations/k18',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k18',JEIfile),JEI);
                elseif Gaincontrol ==19
                    parsave(fullfile(pwd, 'Simulations/k19',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k19',JEIfile),JEI);
                elseif Gaincontrol ==20
                    parsave(fullfile(pwd, 'Simulations/k20',JIEfile),JIE);
                    parsave(fullfile(pwd, 'Simulations/k20',JEIfile),JEI);
                    
                end

                    %
                    Params=params;
                    Params.Ext=stimulus_save.Ext;
                    %         Params.savedir=[];
                    
                    iTrial = 1;
                    data=PlotData{iTrial}; % membrane potentials
                    firings=firings{iTrial}; % spikes
                    
                    
                    %
                    [Trace1] = aux.fun_PlotTrial(data,firings,Params,J2,i,j,k,m,imjk);
                    
                    
                      SomeTrace(imjk,i,j,m) = {Inh};
                     SomeTrace2(imjk,i,j,m) = {[Odor1 Odor2]};
                    SomeTrace3(imjk,i,j,m) = {Trace1};
                     SomeTrace4(imjk,i,j,m) = {J2};
                    X = [imjk i j k m];
                    disp(X);
%                     toc
%                     close all;
                    
                end
            end
            %             SomeTrace4(m,j) = {tempgauss};
        end
    end

      AllTrace(:,:,:,:,k) = SomeTrace;
    AllTrace2(:,:,:,:,k) = SomeTrace2;
    AllTrace3(:,:,:,:,k) = SomeTrace3;
    AllTrace4(:,:,:,:,k) = SomeTrace4;
end
%
% for ii=1:100
%     aaa(ii,:) = Weightchange(ii+1,:)-Weightchange(ii,:);
% end
% aaa(99,:)
% % figure
% imagesc(J2- JEIinitial)
% 
file_name = 'AllTraceDynamicSyn1.mat';
save(file_name,'AllTrace','-v7.3');
file_name = 'AllTraceDynamicSyn2.mat';
save(file_name,'AllTrace2','-v7.3');
file_name = 'AllTraceDynamicSyn3.mat';
save(file_name,'AllTrace3','-v7.3');
file_name = 'AllTraceDynamicSyn4.mat';
save(file_name,'AllTrace4','-v7.3');

toc
% if rt<1
%     run 'demo.m'
% end
% Maxslopefile = 'MaxSlopeRatio-22-FixedV-Mu-J.mat';
% save(Maxslopefile, 'maxRT');
%  toc
% maxRT = permute(maxRT, [1 3 2])
% figure;
% plot(maxRT(1,:))
%   AvgMaxRatio = median(maxRT,3);
%  plot(AvgMaxRatio(:))
% for i=1:25
%     for j=1:10
%         count=0;
%         for k=1:1000
%             if maxRT(i,j,k)>1
%                 count = count+1;
%             end
%         end
%         CountmaxRT(i,j) = count;
%     end
% end
% figure;
% plot(CountmaxRT(:))
% [M,I]= max(AvgMaxRatio(:))
% [x,y] = find(AvgMaxRatio == 1.9615)
%
% maxRT = permute(maxRT,[3 1 2]);
% a = var(maxRT);
% a = permute(a,[2 3 1]);
% for i = 1:length(x)
%     a(x(i),y(i))
% end
% for i = 1:length(x)
%     AvgMaxRatio(x(i),y(i))
% end
%
%
% plot(AvgMaxRatio(x,y))
%
% [x,y] = find(AvgMaxRatio > 1.6)
% th = find(AvgMaxRatio > 1.6)
% %
% % mth = min(a(th))
%
% [x1,y1] = find(a(th)<0.9)
% AvgMaxRatio(x(x1(1)),y(x1(1)))
% AvgMaxRatio(x(x1(2)),y(x1(2)))
% a(x(x1(1)),y(x1(1)))
% a(x(x1(2)),y(x1(2)))
% figure 1 - rasterplot of all neurons in trial
% figure 2 - time course of membrane potential and PSC traces for E and I representative neurons
% figure 3 - time course of firing rate in clusters
% figure 4 - time course of stimuli (with CSgauss stimulus, the cue profile should be multiplied by a factor drawn from figure 5, one for each neuron
% figure 5 (with CSgauss stimulus only) - across-neurons distribution of cue peak values
