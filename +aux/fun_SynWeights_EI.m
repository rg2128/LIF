function [J, params]=fun_SynWeights_EI(paramsfile,i,CountGo,CountNeut,j,k,imjk,Inh,EIC)

% LOAD PARAMETERS
Gaincontrol = k;
CondIter = j;
ImageCount = i;
params=paramsfile;
aux.v2struct(params);
Network=params.Network;
aux.v2struct(Network);
Sim=params.Sim;
%-----------------------
% PARAMETERS VALUES
%-----------------------
numfig=1;
Next=N_e; % external units
% CLUSTERS
Q=p; % number of clusters
%-----------------------
% SYNAPTIC WEIGHTS
%-----------------------
% WEIGHTS
% depression
gam=1/(2-f*(Q+1));%
% Jminus = 1.-gam*f*(Jplus-1.);
Jminus = 1;
params.Jminus = Jminus;
if any(strcmp(clusters,'EI'))
    %     JminusEI = 1.-gam*fI*(JplusEI-1.);
    JplusEI = 1/(1/p+(1-1/p)/factorEI);
    JminusEI = JplusEI/factorEI;
    params.JminusEI = JminusEI;
    jei_out=-JminusEI*Jei; % intra-cluster
    jei_in=-JplusEI*Jei; % inter-cluster
    %     fprintf('JplusEI=%0.03g, JminusEI=%0.03g\n',JplusEI,JminusEI);
end
if any(strcmp(clusters,'IE'))
    %     JminusIE = 1.-gam*fI*(JplusIE-1.);
    JplusIE = 1/(1/p+(1-1/p)/factorIE);
    JminusIE = JplusIE/factorIE;
    params.JminusIE = JminusIE;
    jie_out=JminusIE*Jie; % intra-cluster
    jie_in=JplusIE*Jie; % inter-cluster
    %     fprintf('JplusIE=%0.03g, JminusIE=%0.03g\n',JplusIE,JminusIE);
end
if any(strcmp(clusters,'II'))
    JplusII=factorII;
    JminusII = 1.-gam*fI*(JplusII-1.);
    params.JminusII = JminusII;
    jii_out=-JminusII*Jii; % intra-cluster
    jii_in=-JplusII*Jii; % inter-cluster
    %     fprintf('JplusII=%0.03g, JminusII=%0.03g\n',JplusII,JminusII);
end

%
jee=Jee;
jee_out=Jminus*Jee; % intra-cluster potentiation
jee_in=Jplus*Jee; % inter-cluster depression
jei=-Jei;
jie=Jie;
jii=-Jii;
% connection probability
pee=Cee/N_e;
pie=Cie/N_e;
pei=Cei/N_i;
pii=Cii/N_i;
pext=Cext/Next;
peeout=pee;
peein=pee;
peiout=pei;
peiin=pei;
pieout=pie;
piein=pie;
piiout=pii;
piiin=pii;

% SYNAPTIC MATRIX
%----------------------------
% if delta>0, generate a distribution of synaptic weights with mean J and
% variance delta^2 J^2
peeout=pee;
peein=pee;
% check #clusters and coding level are consistent
NcUnits=round(f*N_e);    %  number of Exc units per cluster
Numbg=round(N_e*(1-f*p));% number of background (i.e. non-selective) Exc units

switch Network.clust
    case 'hom'
        popsize=repmat(NcUnits,Q,1);
    case 'het'
        Nc=[];
        clust_std=Network.clust_std;
        while (sum(Nc)-(N_e-Numbg))~=0 || any(Nc<0)
            Nc=round(NcUnits+(NcUnits*clust_std)*randn(Q,1));
        end
        popsize=Nc; % array of cluster sizes
        if any(sum(popsize)-(N_e-Numbg))
            fprintf('\n---ERROR: Heterogeneous clusters: Problem with cluster numbers\n');
        end
end

cusumNcE=[0 cumsum(popsize)'];
% background units (if present), if not, override in next line
JEE=(jee*(ones(N_e)+delta*randn(N_e,N_e))).*(rand([N_e,N_e])<peeout);
JEI=(jei*(ones(N_e,N_i)+deltaEI*randn(N_e,N_i))).*(rand([N_e,N_i])<pei);
JIE=(jie*(ones(N_i,N_e)+deltaIE*randn(N_i,N_e))).*(rand([N_i,N_e])<pie);
JII=(jii*(ones(N_i)+delta*randn(N_i,N_i))).*(rand([N_i,N_i])<pii);
clustermatrix=eye(Q);



popsizeI=repmat(NcUnits,Q,1);

cusumNcI=[0 cumsum(popsizeI)'];


if any([strcmp(clusters,'EI'), strcmp(clusters,'EI'), strcmp(clusters,'II')])
    NcUnits=round(fI*N_i);    %  number of Exc units per cluster
%     fprintf('  --- Synaptic weights: %d units/cluster \n',NcUnits);
    Numbg=round(N_i*(1-fI*p)); % number of background (i.e. non-selective) Exc units
%     fprintf('  --- fraction of bg Inh units: %0.03g',Numbg/N_i);
    if any(strcmp(Network.clusters,'II'))
        % background units (if present), if not, override in next line
        if strcmp(Network.clustII,'het') || strcmp(Network.clustII,'hom')
            % clustered units: inter-cluster weights
            JII(1:cusumNcI(Q+1),1:cusumNcI(Q+1))=...
                (jii_out*(ones(cusumNcI(Q+1),cusumNcI(Q+1))+deltaII*randn(cusumNcI(Q+1),...
                cusumNcI(Q+1)))).*(rand([cusumNcI(Q+1),cusumNcI(Q+1)])<piiout); % inter-cluster weights
            for clu=2:Q+1 % intra-cluster higher weights
                JII(1+cusumNcI(clu-1):cusumNcI(clu),1+cusumNcI(clu-1):cusumNcI(clu))=...
                    (jii_in*(ones(popsizeI(clu-1),popsizeI(clu-1))+deltaII*randn(popsizeI(clu-1),popsizeI(clu-1)))).*...
                    (rand([popsizeI(clu-1),popsizeI(clu-1)])<piiin);
            end
        end
    end
    params.popsizeI=popsizeI;
end


for j=1:length(cusumNcE)-1
    for k = 1:size(Inh,2)
        if Inh(j,k)>0
            JEI(cusumNcE(Inh(j,k))+1:cusumNcE(Inh(j,k)+1),cusumNcE(j)+1:cusumNcE(j+1)) =  (jei_in*(ones(popsize(1),popsizeI(1))+deltaEI*randn(popsize(1),popsizeI(1)))).*...
                (rand([popsize(1),popsizeI(1)])<peiin);
            JIE(cusumNcE(j)+1:cusumNcE(j+1),cusumNcE(Inh(j,k))+1:cusumNcE(Inh(j,k)+1)) = (jie_in*(ones(popsizeI(1),popsize(1))+deltaIE*randn(popsizeI(1),popsize(1)))).*...
                (rand([popsizeI(1),popsize(1)])<piein);
        end
    end
end
% if ImageCount ==2
% 
%     for j=1:length(cusumNcE)-1
%         if params.Stimulus.feat(1).selective(j)==1
%         for k = 1:size(Inh,2)
%             if Inh(j,k)>0
%                 JEI(cusumNcE(Inh(j,k))+1:cusumNcE(Inh(j,k)+1),cusumNcE(j)+1:cusumNcE(j+1)) =  1*(jei_in*(ones(popsize(1),popsizeI(1))+deltaEI*randn(popsize(1),popsizeI(1)))).*...
%                     (rand([popsize(1),popsizeI(1)])<peiin);
%                 JIE(cusumNcE(j)+1:cusumNcE(j+1),cusumNcE(Inh(j,k))+1:cusumNcE(Inh(j,k)+1)) = 1*(jie_in*(ones(popsizeI(1),popsize(1))+deltaIE*randn(popsizeI(1),popsize(1)))).*...
%                     (rand([popsizeI(1),popsize(1)])<piein);
%             end
%         end
% 
%     end
% end
% end

Inhfile = 'Inhfix.mat';
JIEfile = 'JIEfix2.mat';
JEIfile = 'JEIfix2.mat';

JEEfile = 'JEEfix2.mat';
JIIfile = 'JIIfix2.mat';


if ImageCount ==1 && imjk == 1 && CondIter == 1 && EIC == 1
    formatSpec = 'Simulations/k%d';
    file_name = sprintf(formatSpec,Gaincontrol);
    parsave(fullfile(pwd,file_name,JEEfile),JEE);
    parsave(fullfile(pwd,file_name,JIIfile),JII);
end    



if Gaincontrol ==1 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k1',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k1',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==2 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k2',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k2',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==3 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k3',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k3',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==4 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k4',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k4',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==5 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k5',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k5',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==6 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k6',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k6',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==7 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k7',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k7',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==8 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k8',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k8',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==9 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k9',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k9',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==10 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k10',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k10',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==11 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k11',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k11',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==12 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k12',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k12',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==13 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k13',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k13',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==14 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k14',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k14',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==15 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k15',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k15',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==16 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k16',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k16',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==17 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k17',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k17',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==18 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k18',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k18',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
elseif Gaincontrol ==19 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k19',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k19',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
elseif Gaincontrol ==20 && ImageCount>1
    JIE = par_load(fullfile(pwd, 'Simulations/k20',JIEfile));
    JEI = par_load(fullfile(pwd, 'Simulations/k20',JEIfile));
    JIE= JIE.x;
    JEI= JEI.x;
    
end





params.popsizeI=popsizeI;

J=[JEE JEI; JIE JII];
% 
% figure;
% % subplot(2,1,1);
% colormap(parula); %xy=J; fun_colormapLim;
% % imagesc(J(500:1000,1:500));
% aux.figset(gca,'neurons','neurons','weights',10);
% colorbar;
% % subplot(2,1,2);
% % lambda=eig(J);
% % plot(real(lambda),imag(lambda),'.');
% % aux.figset(gca,'Re(\lambda)','Im(\lambda)','eig(weights)',10);
% formatSpec = 'Weights-%d--%d--%d--%d.fig';
% file_name = sprintf(formatSpec,ImageCount,k1,j1,m1);
% saveas(gcf,fullfile('Simulations',file_name),'fig');

params.popsize=popsize;
params.clustermatrix=clustermatrix;

theta=[params.theta_e, params.theta_i];
Theta=zeros(1,numel(params.popsize)+2);
Theta(1:numel(params.popsize)+1)=theta(1);
Theta(end) = theta(2);
params.Theta=Theta;