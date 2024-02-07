
    function PlotMembrane(Sim,nplot,vplot,firings,ind_plot,VTh)
        
        dt=Sim.dt_step;
        colors=['b'; 'r'];
        hh=[];
        bins=min([Sim.t_End-Sim.plot_length,Sim.t_Start])+dt:dt:Sim.t_End;
        for j=1:nplot
            hh(j)=plot(bins, vplot(j,1:numel(bins)),'color',colors(j,:)); % membrane potentials
            hold on
            if any(any(firings))
                ind=find(firings(:,2)==ind_plot(j));
                if ~isempty(ind) % plot spike overshoots
                    h=line([firings(ind,1), firings(ind,1)]',...
                        [VTh(ind_plot(j))*ones(size(ind)) VTh(ind_plot(j))*ones(size(ind))*2]');
                    set(h,'color',colors(j,:));
                    hold on
                end
            end
            plot([bins(1) bins(end)],[VTh(ind_plot(j)) VTh(ind_plot(j))],'color',colors(j,:),'linestyle','--');
            hold on
        end
        
        legend(hh,'E','I')
        ylim([1.2*min(min(vplot(:,0.2/dt:end))) 2*max(VTh(ind_plot))]);
        xlim([bins(1) bins(end)]);
        xlab='Time (s)';
        ylab='mV';
        tt='';
        fntsz=15;
        aux.figset(gca,xlab,ylab,tt,fntsz);
        
    end