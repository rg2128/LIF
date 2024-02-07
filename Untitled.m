function [J, params]=fun_SynWeights_EI(paramsfile)

% LOAD PARAMETERS

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
% gam=1/(2-f*(Q+1));%
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
NcUnits=1;    %  number of Exc units per cluster
Numbg=0; % number of background (i.e. non-selective) Exc units
popsize=repmat(NcUnits,Q,1);

cusumNcE=[0 cumsum(popsize)'];

JEE = zeros(N_e,N_e);
JII = zeros(N_i,N_i);
JEI = zeros(N_i,N_e);
JIE = zeros(N_e,N_i);

clustermatrix=eye(Q);


for i = 1:N_e
    temp = [i i i i]
    perm(i,:) = temp
end
Inh = zeros(N_i,Q/N_i);
for j = 1:N_i
    for k =1:Q/N_i
        pos = randi(length(perm(:)));
            while perm(pos) == 0 || any(Inh(j,:) == perm(pos))
                pos = randi(length(perm(:)));
                if nnz(Inh(j,:))>2
                    break;
                end
            end
            Inh(j,k) = perm(pos);
            perm(pos) = 0;
    end
end

G= graph;
for j = 1:N_i
    s = [];
    t = [];
    s(1,1:nnz(Inh(j,:)))  = j+20;
    t(1,:) = Inh(j,1:length(s));
    G = addedge(G,s,t)
%     p = plot(G)
%     if j<=20
        
end
p = plot(G);
% a(1,1:20) = 1;
% a(1,21:40) = 5;
% G.Nodes.NodeColors = a';
% p.NodeCData = G.Nodes.NodeColors;
% colorbar
highlight(p,[21:40], 'NodeColor' , 'r');


for j=1:N_i
    for k = 1:Q/N_i
        if Inh(j,k)>0
            JEI(j,Inh(j,k)) =  1;
        end
    end
end
JEI = Jei*JEI;  
JIE = -transpose(JEI);

params.popsizeI=popsizeI;

J=[JEE JEI; JIE JII];

figure(1); clf; 
% subplot(2,1,1);
colormap(parula); %xy=J; fun_colormapLim;
imagesc(J);
aux.figset(gca,'neurons','neurons','weights',10);
colorbar;
% subplot(2,1,2);
% lambda=eig(J);
% plot(real(lambda),imag(lambda),'.');
% aux.figset(gca,'Re(\lambda)','Im(\lambda)','eig(weights)',10);
% saveas(gcf,fullfile('data','weights.pdf'),'pdf');
        
params.popsize=popsize;
params.clustermatrix=clustermatrix;

theta=[params.theta_e, params.theta_i];
Theta=zeros(1,numel(params.popsize)+2);
Theta(1:numel(params.popsize)+1)=theta(1);
Theta(end) = theta(2);
params.Theta=Theta;