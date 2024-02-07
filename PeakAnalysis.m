
% 
% for i = 1 : length(AllPeaks)
%     
%     for j = 1 : length(AllPeaks{1,1}{1,1}{1})
%         Glom(i,j) = AllPeaks{1,i}{1,1}{1}(j);
%     end
%     for j = 1 : length(AllPeaks{1,i}{1,1}{2})
%         DAT(i,j) = AllPeaks{1,i}{1,1}{2}(j);
%     end
% end  

for i = 1 : length(AllPeaks(1,:,1))
    for j = 1 : length(AllPeaks(1,1,:))
        temp = AllPeaks(:,i,j);
        for k = 1 : length(AllPeaks)
            for m = 1 : length(temp{k,1}{1,1})
                Glom(k,m) = temp{k,1}{1,1}(m);
            end
            for m = 1 : length(temp{k,1}{1,2})
                DAT(k,m) = temp{k,1}{1,2}(m);
            end
        end
        
        dur = 0.05+ 0.05*(j-1);
        start = -0.5 + 0.05*(i-1);
        stop = start+dur;
        AvgGlom(:,i,j) = mean(Glom,2);
        AvgGlomGo(:,i,j) = AvgGlom(1:2:80,i,j);
        AvgGlomNeut(:,i,j) = AvgGlom(2:2:80,i,j);
        AvgGlomGo(:,i,j) = AvgGlomGo(:,i,j)/max(AvgGlomGo(:,i,j));
        AvgGlomNeut(:,i,j) = AvgGlomNeut(:,i,j)/max(AvgGlomNeut(:,i,j));
        figure
        plot(AvgGlomGo(:,i,j))
        hold on 
        plot(AvgGlomNeut(:,i,j))
        hold off
        title (['start= ' num2str(start)]) 
        xlabel (['stop= ' num2str(stop)])
%         formatSpec = 'GoNeutCompareStart%0.03gStop%0.03g.fig';
%         file_name = sprintf(formatSpec,start,stop);
%         saveas(gcf,fullfile('Simulations',file_name),'fig');
%         close all;
    end
end






  

plot(AvgGlomGo)
hold on 
plot(AvgGlomNeut)