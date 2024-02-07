    function PlotAllNeurons(Ext,firings,IN,imagecount,Tone,breathgain,inhgain,copies)
%         i,j,k,m,imjk
        hold on;
        % unpack
        aux.v2struct(IN);
        %
        cusumNcE=[0 cumsum(NcE)'];
        plot(firings(firings(:,2)<=Ne_plot,1),firings(firings(:,2)<=Ne_plot,2),'k.','markersize',1);
        plot(firings(firings(:,2)>N_e & firings(:,2)<N_e+Ni_plot,1),firings(firings(:,2)>N_e & firings(:,2)<N_e+Ni_plot,2)-N_e+Ne_plot,'r.','markersize',1);
%         plot(firings(firings(:,2)<=100,1),firings(firings(:,2)<=100,2),'k.','markersize',1);
%         hold on
%          plot(firings(firings(:,2)>1500 & firings(:,2)<1600,1),firings(firings(:,2)>1500 & firings(:,2)<1600,2)-N_e+Ne_plot,'r.','markersize',1);
%      xlim([0, 0.5])
if any(strcmp(fieldnames(Ext),'stim'))
            for n=1:numel(Ext.stim)
                time_bin=[Ext.stim(n).interval(1) Ext.stim(n).interval(2)];
                if ~isempty(Ext.stim(n).ind) && strcmp(Ext.stim(n).name,'US')
                    % find pops with target neurons
                    pop_ind=zeros(1,numel(cusumNcE)-1);
                    for k=1:numel(cusumNcE)-1
                        pop_ind(k)=any(Ext.stim(n).ind>=cusumNcE(k)+1 & Ext.stim(n).ind<=cusumNcE(k+1));
                    end
                    pop_ind=find(pop_ind);
                    for k=1:numel(pop_ind)
                        lower=(cusumNcE(pop_ind(k)));
                        upper=cusumNcE(pop_ind(k)+1);
                        indspikes=firings(:,2)>=lower & firings(:,2)<=upper;
                        indstim=firings(:,1)>=Ext.stim(n).interval(1) & firings(:,1)<Ext.stim(n).interval(2);
                        indspikes=indspikes & indstim;
                        plot(firings(indspikes,1),firings(indspikes,2),'m.','markersize',1);
                        hold on;
                    end
                end
            end
        end
        line([0 0],[0 Ne_plot+Ni_plot],'color','b');
        hold on;
        ylim([0 Ne_plot+Ni_plot]);
        xlab='Time (s)';
        ylab='Neurons';
        tt='All neurons ordered by population';
        aux.figset(gca,xlab,ylab,tt,15);
        
        xlim(XLIM);
        hold off
         if copies == 1
                formatSpec = 'Raster-%d--%d--%d--%d----%d.fig';
                file_name = sprintf(formatSpec,imagecount,Tone,breathgain,inhgain,copies);
                saveas(gcf,fullfile('Simulations',file_name),'fig');
         end
    end