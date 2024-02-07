i = 1;

%for i = 1 : length(AllTrace(1,:,1))
    for j = 14 : length(AllTrace(1,1,:))
        temp = AllTrace(:,i,j);
        for k = 1 : length(AllTrace)
            for m = 1 : length(temp{k,1}{1,1}(:,1))
                Glom(m,:,k) = temp{k,1}{1,1}(m,:);
            end
        end
    end
%end

for i =1:2
    for j = 1:80
Glom(i,:,j) = Glom(i,:,j)/max(Glom(i,:,j));
    end
end


ActivationFrame=[];
distances=[];
PeakAmplitude=[];
activeglomeruli=[];
for j=1:length(Glom(1,1,:))               %For the number of trials
    count=1;
    Grad = [];
    for k=1:length(Glom(:,1,1))           %For the number of glomeruli in that particular trial
        distances=[];
        [val,idx]= max(Glom(k,1:68,j));% Calculate peak in this stringent range. Stringent because the range is more selective than GIA
        PeakAmplitude(j,k)= val;
        for i=20:idx
            distances(i)=abs(val/2-Glom(k,i,j));
        end
        for m =1:length(distances)
            if distances(m)==0
                distances(m)=10;                                %if distance is 0, something went wrong set distance to high number so that it isn't captured later
            end
        end
        [val,idx2]= min(distances(:));               %Selecting the Index for minimum distance only if the difference between actual half peak and expected half peak is less than 1%. Since it isn't a continuous vector...
        ActivationFrame(j,k)=idx;
%         hm = ActivationFrame(j,:);
%         hm = mean(hm);
           
        Grad(k) = mean(gradient(Glom(k,31:33,j)));
    end
    FinGrad(j) = mean(Grad);

end

ActivationGo = ActivationFrame(1:2:80,:);
ActivationNoGo = ActivationFrame(2:2:80,:);

% FinGradGo = FinGrad(1:2:80);
% FinGradNoGo = FinGrad(2:2:80);
for m =1:2
FinGradGo(m,:,:) = (Glom(m,30:40,1:2:80));
FinGradNoGo(m,:,:) = (Glom(m,30:40,2:2:80));
end
FinGradGoAvg  = mean(FinGradGo,1);
FinGradNoGoAvg= mean(FinGradNoGo,1);
for k =1:40
    X(1,:) = FinGradGoAvg(1,:,k);
    X(2,:) = FinGradNoGoAvg(1,:,k);
D(k)=pdist(X)
end
plot(D)

DiffAct = ActivationNoGo - ActivationGo;
DiffGrad = FinGradGoAvg - FinGradNoGoAvg
plot(DiffGrad)
plot(mean(DiffAct,2))











%         dur = 0.05+ 0.05*(j-1)
%         start = -0.5 + 0.05*(i-1)
%         stop = start+dur;
        for n = 1:2:40
        glomGo = mean(Glom(:,:,n));
        glomNoGo = mean(Glom(:,:,n+1));
        glomGo = glomGo/max(glomGo);
        glomNoGo = glomNoGo/max(glomNoGo);
        figure;
        plot(glomGo)
        hold on
        plot(glomNoGo)
        hold off
        end
%
% 


         figure;
        plot(Glom(2,:,33))
        hold on
        plot(Glom(2,:,34))
        
%
%         figure;
%         plot(Glom(1,:,30))
%         hold on
%         plot(Glom(2,:,30))
%         ylim([0 25])
%         AvgGlom(:,i,j) = mean(Glom,2);
%         AvgGlomGo(:,i,j) = AvgGlom(1:2:80,i,j);
%         AvgGlomNeut(:,i,j) = AvgGlom(2:2:80,i,j);
%         AvgGlomGo(:,i,j) = AvgGlomGo(:,i,j)/max(AvgGlomGo(:,i,j));
%         AvgGlomNeut(:,i,j) = AvgGlomNeut(:,i,j)/max(AvgGlomNeut(:,i,j));
%         figure
%         plot(AvgGlomGo(:,i,j))
%         hold on
%         plot(AvgGlomNeut(:,i,j))
%         hold off
%         title (['start= ' num2str(start)])
%         xlabel (['stop= ' num2str(stop)])
%         formatSpec = 'GoNeutCompareStart%0.03gStop%0.03g.fig';
%         file_name = sprintf(formatSpec,start,stop);
%         saveas(gcf,fullfile('Simulations',file_name),'fig');
%         close all;
end
end








plot(AvgGlomGo)
hold on
plot(AvgGlomNeut)