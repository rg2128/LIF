
    function PlotPopRaster(rate,x_bins,bins)
        
        
        Cmin=min(min(rate));
        Cmax=max(max(rate));
        imagesc(x_bins,1:size(rate,1),rate); axis xy;
        caxis([Cmin, Cmax]);
        % colormap gray;
        % colormap(1-colormap);
        xlim([bins(1) bins(end)]);
        t=colorbar; get(t,'ylabel');
        set(get(t,'ylabel'),'String', 'Firing rate [spks/s]');
        hold off
        xlab='Time [s]';
        ylab='Population index';
        tt='';
        fntsz=15;
        aux.figset(gca,xlab,ylab,tt,fntsz);
        
    end