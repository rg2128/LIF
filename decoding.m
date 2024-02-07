%% imjk,i,m,j,k
%% JEI = J2(1:500,501:1000);
%%   SomeTrace2(imjk,i,j,m) = {[Odor1 Odor2]};
for j=1:2
count = 1;
AvgAUC = [];
i=1;

% AvgAUCgoAll = [];
% AvgAUCnogoAll = [];
% diffAUCavgAll = [];
 aaaa=1:16;
    
for k =aaaa
AvgAUCgoAll2 = [];
AvgAUCnogoAll2 = [];
diffAUCavgAll2 = [];PeakAmplitude =[];
    for m=1:2:100
    tracedat = [];
    traceglom = [];
    tracedat2 = [];
    traceglom2 = [];
   
        
        % tracedat = [];
        % traceglom = [];
        % tracedat2 = [];
        % traceglom2 = [];
        aa= AllTrace2(1,m,j,i,k);
        goAll = [];
        nogoAll = [];
        go=find(aa{1,1}(1:50)>0);
        nogo=find(aa{1,1}(51:100)>0);
        
        for q=1:length(go)
            goAll(q,:) = go(q)*10-9:go(q)*10;
        end
        for q=1:length(nogo)
            nogoAll(q,:) = nogo(q)*10-9:nogo(q)*10;
        end
        
        traceGlomdat = AllTrace3{1,m,j,i,k};
        traceGlomdat2 = AllTrace3{1,m+1,j,i,k};
        traceglom = (traceGlomdat{1,1});
        tracedat = (traceGlomdat{1,2});
        traceglom2= (traceGlomdat2{1,1});
        tracedat2= (traceGlomdat2{1,2});
        
        [~,LatencyGo] = max(mean(traceglom(:,100:150)));
        [~,LatencyNoGo] = max(mean(traceglom2(:,100:150)));
        for indiglom = 1:size(traceglom,1)
            PeakAmplitude(m,indiglom) = max((traceglom(indiglom,100:150)));
            PeakAmplitude(m+1,indiglom) = max((traceglom2(indiglom,100:150)));
        end
        
        
         Avg_traceglom(k,:) = mean(traceglom);
        Avg_tracedat(k,:) = mean(tracedat(:,:));
        Avg_traceglom2(k,:)= mean(traceglom2);
        Avg_tracedat2(k,:)= mean(tracedat2(:,:));
        %         figure;
        %     plot(traceglom(k,:));
        %     hold on
        %     plot(traceglom2(k,:));
        %         figure;
        %     stdshade(traceglom,[],'b');
        %     hold on
        %     stdshade(traceglom2);
        %     xlim([90 150])
        %
        %     figure;
        %     stdshade(tracedat,[],'b');
        %     hold on
        %     stdshade(tracedat2);
        %         xlim([90 150])
        
    end
    
    count = count+1;
    



for slide = 1:81
PeakAmplitudeGo = PeakAmplitude(1 + 1*(slide-1) : 2 : 20 + 1*(slide-1),:);
PeakAmplitudeNoGo = PeakAmplitude(2 + 1*(slide-1) : 2 : 20 + 1*(slide-1),:);
PeakAmplitudetemp = [PeakAmplitudeGo; PeakAmplitudeNoGo];
% AvgAUCgo(count,:) =AUC(:,1);
Index.Go = 1:10;
Index.NoGo =11:20;
% Index = [Index.Go Index.NoGo];
Name = 'gain1';
c1 = [ 0 0 1];
c2 = [ 0 1 0];
c = [repmat(c1,10,1) ;repmat(c2,10,1)];

Output2 = PCA_Peak(PeakAmplitudetemp,c,Index,Name);

if j == 1
    OutputNT(k,slide) = Output2{1,2};
else
    OutputT(k,slide) = Output2{1,2};
end

end
end
end
figure;
stdshade(1-OutputNT,[],'g')
hold on
stdshade(1-OutputT)
% plotweights(1,:) = mean(AvgAUCgoAll2');
% plotweights(2,:) = std(AvgAUCgoAll2')/sqrt(48);
% plotweights(3,:) = mean(AvgAUCnogoAll2');
% plotweights(4,:) = std(AvgAUCnogoAll2')/sqrt(48);
% plotweights(5,:) = mean(diffAUCavgAll2');
% plotweights(6,:) = std(diffAUCavgAll2')/sqrt(48);