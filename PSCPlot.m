
    function PSCPlot(Sim,Ext,iIplot,iEplot,iExtplot,hicutoff,srate,VTh,ind_plot,tau)
        stimcolors={'m',[0.5,0.5,0.5]};
        nplot=min(2,size(iIplot,1))+1;
        dt=Sim.dt_step;
        hE=[]; hI=[]; hext=[]; htot=[]; hth=[];
        t_plot=max([Sim.t_End-Sim.plot_length,Sim.t_Start])+dt:dt:Sim.t_End;
        ExcInh={'Exc','Inh'};
        for j=1:nplot-1
            subplot(nplot,1,j+1)
            hold on;
            [iISmooth]=aux.eegfilt(iIplot(j,:),srate,0,hicutoff);
            [iESmooth]=aux.eegfilt(iEplot(j,:),srate,0,hicutoff);
            [iExtSmooth]=aux.eegfilt(iExtplot(j,:),srate,0,hicutoff);
            hE(j)=plot(t_plot, iESmooth,'color','b','linewidth',1); % EPSC
            hI(j)=plot(t_plot, iISmooth,'color','r','linewidth',1); % IPSC
            hext(j)=plot(t_plot, iExtplot(j,:),'color','g','linewidth',1); % EPSC
            itot=iEplot(j,:)+iIplot(j,:)+iExtplot(j,:);
            hth(j)=plot([t_plot(1) t_plot(end)],ones(1,2)*VTh(ind_plot(j))/tau(ind_plot(j)),'color','b','linewidth',1,'linestyle','-.');
            [itotSmooth]=aux.eegfilt(itot,srate,0,hicutoff);
            htot(j)=plot(t_plot,itotSmooth,'color','k','linewidth',1); % membrane potentials
            grid on
            if any(strcmp(fieldnames(Ext),'stim'))
                for n=1:numel(Ext.stim)
                    time_bin=[Ext.stim(n).interval(1) Ext.stim(n).interval(2)];
                    lower=min(min(iISmooth))*ones(1,numel(time_bin));
%                     upper=max(max(iESmooth))*ones(1,numel(time_bin));
%                     [~,~]=aux.jbfill(time_bin,lower,upper,'m',0,0,0.2);
                    hstim(n)=line(time_bin([1,end]),(1+0.1*n)*lower([1,end]),'color',stimcolors{n},'linewidth',5);
                end
            end
            xlab='Time (s)';
            ylab='mV/s';
            tt=sprintf('PSC to %s unit',ExcInh{j});
            fntsz=15;
            aux.figset(gca,xlab,ylab,tt,fntsz);
            ind_lim=(t_plot>Sim.t_Start+0.1 & t_plot<Sim.t_End-0.1);
            legend([hE(1) hI(1) hext(1) htot(1) hth(1), hstim(1),hstim(2)],'i_E','i_I','i_{ext}','i_{tot}','V_{th}/\tau','CS','US')
            xlim([t_plot(1) t_plot(end)]);
            ylim([1.5*min(min(iISmooth(:,ind_lim)))...
                1.5*max([max(max(iESmooth(:,ind_lim))) max(max(iExtSmooth(:,ind_lim)))])]);
        end
        
    end