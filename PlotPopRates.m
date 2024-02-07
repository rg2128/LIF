function [rate, x_bins, Trace]=PlotPopRates(bins,firings,N_e,N_i,NcE,J,i,ExcCells,Tone,breathgain,inhgain,copies)
ActiveDAT = J(1:500,501:1000);
ImageCount = i;
MinRate=10;
BinSize=diff(bins(1:2));
fntsz=15;
cols=aux.distinguishable_colors(numel(NcE)+1);
binwidth=8;
cusumNcE=[0 cumsum(NcE)'];
rate=zeros(numel(NcE),numel(bins)+1);
count = 1;
PeakGlom = [];
PeakDAT = [];
TraceGlom =[];
TraceDAT = [];
 figure;  
for clu=1:numel(NcE)
%     if ExcCells(clu)
    ind_clu=(firings(:,2)>cusumNcE(clu) & firings(:,2)<=cusumNcE(clu+1));
    spikes=firings(ind_clu,1);
    x=histc(spikes,bins);
%     x = x/max(x(:));
    [tout, z]=aux.gaussfilt(bins,x,binwidth,0);
    x_bins=tout;
    rate(clu,:)=z/(NcE(clu)*BinSize);
%     rate(clu,:) = rate(clu,:)/max(rate(clu,:));
%     PeakGlom(count) =  max(rate(clu,:));
    TraceGlom(count,:) =  (rate(clu,:));
    count = count +1;
    a1 = rate(clu,:);
    plot(x_bins-diff(x_bins(1:2))/2,rate(clu,:),'color',cols(clu,:),'linewidth',1.5);
    hold on;
%     end
    
end
xlim([-1 1])
%          if copies == 1
%                 formatSpec = 'TraceGlom-%d--%d--%d--%d----%d.fig';
%                 file_name = sprintf(formatSpec,ImageCount,Tone,breathgain,inhgain,copies);
%                 saveas(gcf,fullfile('Simulations',file_name),'fig');
%          end
% hold off;
% formatSpec = 'TraceGlom-%d.fig';
% file_name = sprintf(formatSpec,ImageCount);
% saveas(gcf,fullfile('Simulations',file_name),'fig');

% DAT = zeros(length(ActiveDAT),10);
% for j = 1:10
%     for i= 1:length(ActiveDAT)
%         
%         if ActiveDAT(j,i)<0
%             DAT(i,j) = 1;
%         end
%     end
% end
count =1;
 figure;
for clu=1:numel(NcE)
    ind_clu=(firings(:,2)>cusumNcE(clu)+500 & firings(:,2)<=cusumNcE(clu+1)+500);
    spikes=firings(ind_clu,1);
    x=histc(spikes,bins);
%     x = x/max(x(:));
    [tout, z]=aux.gaussfilt(bins,x,binwidth,0);
    x_bins=tout;
    rate(clu,:)=z/(NcE(clu)*BinSize);
%     rate(clu,:) = rate(clu,:)/max(rate(clu,:));
%     PeakDAT(count) =  max(rate(clu,:));
    TraceDAT(count,:) =  (rate(clu,:));
    count = count +1;
    a1 = rate(clu,:);
    plot(x_bins-diff(x_bins(1:2))/2,rate(clu,:),'color',cols(clu,:),'linewidth',1.5);
    hold on;
    xlim([-0.5 0.5])
end
xlim([-1 1])

% hold off;
%          if copies == 1
%                 formatSpec = 'TraceDAT-%d--%d--%d--%d----%d.fig';
%                 file_name = sprintf(formatSpec,ImageCount,Tone,breathgain,inhgain,copies);
%                 saveas(gcf,fullfile('Simulations',file_name),'fig');
%          end

% Peaks = {PeakGlom PeakDAT};
Trace = {TraceGlom TraceDAT};
%         mid = round(length(bins)/2);
%          SlopeGlom1 = gradient(a1(mid+1:mid+4));
%          SlopeGlom2 = gradient(a1(mid+5:mid+8));
%          SlopeDAT1 = gradient(a2(mid+1:mid+4));
%          SlopeDAT2 = gradient(a2(mid+5:mid+8));
%          x1 = [SlopeGlom1;SlopeDAT1];
%          x2 = [SlopeGlom2;SlopeDAT2];
%          rt = pdist(x2)/pdist(x1);
%          fprintf(' %0.03g\n',rt);
%          slopefile = 'SlopeRatio2.mat';
%         save(slopefile, 'rt');
% unclustered subpop
%         if N_e-cusumNcE(end)>0 % if there's an unclustered population
%             ind_unclu=(firings(:,2)>cusumNcE(end) & firings(:,2)<N_e);
%             spikes=firings(ind_unclu,1);
%             x=histc(spikes,bins);
%             [tout, z]=aux.gaussfilt(bins,x,binwidth,0);
%             x_bins=tout;
%             rate=[rate; z'/((N_e-cusumNcE(end))*BinSize)];
%             plot(x_bins-diff(x_bins(1:2))/2,rate(clu+1,:),'color',cols(clu+1,:),'linewidth',1.5,'linestyle','--');
%         end
%         % inhibitory neurons
%         ind_inh=(firings(:,2)>N_e);
%         spikes=firings(ind_inh,1);
%         x=histc(spikes,bins);
%         [tout, z]=aux.gaussfilt(bins,x,binwidth,0);
%         x_bins=tout;
%         rate=[rate; z'/((N_i)*BinSize)];
%         plot(x_bins-diff(x_bins(1:2))/2,rate(clu+2,:),'color','r','linewidth',2,'linestyle',':');
%
%         line([x_bins(1)-diff(x_bins(1:2))/2 x_bins(end)],[MinRate MinRate],'linestyle','--','color','k');
%         hold on;
%         xlim([bins(1) bins(end)]);
%         xlab='Time [s]';
%         ylab='Firing rate [spks/s]';
%         tt='';
%         aux.figset(gca,xlab,ylab,tt,fntsz);
%         ylim([0 max(max(rate))+10]);

end
