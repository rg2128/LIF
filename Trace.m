%% imjk,i,m,j,k
%% JEI = J2(1:500,501:1000);
%%   SomeTrace2(imjk,i,j,m) = {[Odor1 Odor2]};


AvgAUC = [];
i=1;
AvgAUCgoAll2 = [];
AvgAUCnogoAll2 = [];
diffAUCavgAll2 = [];
AvgAUCgoAll = [];
AvgAUCnogoAll = [];
diffAUCavgAll = [];
AvgAUCnogo = [];
AvgAUCgo = [];
for j= 1:2
    count = 1;
 for m=1:2:100
     tracedat = [];
traceglom = [];
tracedat2 = [];
traceglom2 = [];
aaaa=1:16;
for k =aaaa
%     tracedat = [];
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
%     traceglom(k,:) = mean(traceGlomdat{1,1});
%     tracedat(k,:) = mean(traceGlomdat{1,2});
%     traceglom2(k,:) = mean(traceGlomdat2{1,1});
%     tracedat2(k,:) = mean(traceGlomdat2{1,2});
        traceglom = (traceGlomdat{1,1});
        
    tracedat = (traceGlomdat{1,2});
    traceglom2= (traceGlomdat2{1,1});
    tracedat2= (traceGlomdat2{1,2});
%         traceglom = traceglom(go,:);
%         traceglom2 = traceglom2(nogo,:);
    [~,LatencyGo] = max(mean(traceglom(:,100:150)));
[~,LatencyNoGo] = max(mean(traceglom2(:,100:150)));

%    AUC(k,:)= [max(mean(traceglom(:,100:150))) max(mean(tracedat(:,110:150))) max(mean(traceglom2(:,100:150))) max(mean(tracedat2(:,110:150)))];
%        AUC(k,:)= [LatencyGo max(mean(tracedat(:,100:150))) LatencyNoGo max(mean(tracedat2(:,100:150)))];

   AUC(k,:)= [max(mean(traceglom(:,100:150))) max(mean(tracedat(:,100:150))) max(mean(traceglom2(:,100:150))) max(mean(tracedat2(:,100:150)))];
    
%    AUC(k,:)= [trapz(mean(traceglom(:,100:175))) trapz(mean(tracedat(:,100:175))) trapz(mean(traceglom2(:,100:175))) trapz(mean(tracedat2(:,100:175)))];
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
    AvgAUC = [AvgAUC AUC];

% AvgAUC(count,:) = AUC;
AvgAUCgo(count,:) =AUC(:,1);
AvgAUCnogo(count,:) =AUC(:,3);

%     AvgAUC(count,:) = [mean(AUC) std(AUC)/2];
%     figure;
%     stdshade(traceglom,[],'g');
%     hold on
%     stdshade(traceglom2);
%     
%     xlim([0 200])
%     figure;
%     stdshade(tracedat,[],'g');
%     hold on
%     stdshade(tracedat2);
%             xlim([0 200])
% figure;
%     stdshade(Avg_traceglom,[],'g');
%     hold on
%     stdshade(Avg_traceglom2);
    
    
%     plotweights(1,:) = mean(Avg_traceglom);
% plotweights(2,:) = std(Avg_traceglom)/sqrt(16);
% plotweights(3,:) = mean(Avg_traceglom2);
% plotweights(4,:) = std(Avg_traceglom2)/sqrt(16);
% % %     
% % % %     xlim([0 200])
%     figure;
%     stdshade(Avg_tracedat,[],'g');
%     hold on
%     stdshade(Avg_tracedat2);
    
%      figure;
%     plot(Avg_tracedat);
%     hold on
%     plot(Avg_tracedat2);
%     
%     figure;
%     plot(Avg_traceglom);
%     hold on
%     plot(Avg_traceglom2);
%             xlim([0 200])

count = count+1;

 end
% AvgAUCgo(count,:) =AUC(:,1);

% figure;
% plot(AvgAUC(:,1),'g','LineWidth', 2)
% hold on
% plot(AvgAUC(:,3),'r','LineWidth', 2)

% plotmayb= [mean(AvgAUCgo(1:20,:)); mean(AvgAUCnogo(1:20,:))]

% for z = 1:30
% %         mean((AvgAUC(:,1+z*4-4) - AvgAUC(:,3+z*4-4)))
% %         mean(AvgAUC(:,201+z*4-4) - AvgAUC(:,203+z*4-4))
% 
%     [h,p(z)] = ttest2((AvgAUC(:,1+z*4-4) - AvgAUC(:,3+z*4-4)), (AvgAUC(:,201+z*4-4) - AvgAUC(:,203+z*4-4)))
% end
%     [h,p(z)] = ttest2((AvgAUC(:,5) - AvgAUC(:,7)), (AvgAUC(:,9) - AvgAUC(:,11)))
% 
%    plotmayb= [ (mean(AvgAUC(:,[1 5 9 13 17 21 25 29 33 37]),2)-mean(AvgAUC(:,[3 7 11 15 19 23 27 31 35 39]),2))]
   
   
%    (mean(AvgAUC(:,[201 205 209 213 217]),2)-mean(AvgAUC(:,[203 207 211 215 219]),2)) ]
% 
%    z=24;
%    plotmayb=[(AvgAUC(:,5) - AvgAUC(:,7)) (AvgAUC(:,9) - AvgAUC(:,11))]
%         xlim([90 150])

% end
    
%         xlim([90 150])

% end
for ii=1:length(AvgAUCgo)-5
AvgAUCgo(ii,:) = mean(AvgAUCgo(ii:ii+5,:));
AvgAUCnogo(ii,:) = mean(AvgAUCnogo(ii:ii+5,:));
end
diffAUCavg = AvgAUCgo-AvgAUCnogo;
figure;
stdshade(AvgAUCgo',[],'g')
hold on
stdshade(AvgAUCnogo')
hold on
stdshade(diffAUCavg',[],'b')
% diffAUCTNT = diffAUCavg1- diffAUCavg;
% hold on
% stdshade(diffAUCTNT',[],'k')

%  hold on
% stdshade(diffAUCavg',[],'g')
% diffAUCavg1 = diffAUCavg;

% plotmax(1,:) = mean(AvgAUCgo(20:30,:));
% plotmax(2,:) = mean(AvgAUCnogo(20:30,:));
if j == 2
AvgAUCgoAll2 = [AvgAUCgoAll2 AvgAUCgo];
AvgAUCnogoAll2 = [AvgAUCnogoAll2 AvgAUCnogo];
diffAUCavgAll2 = [diffAUCavgAll2 diffAUCavg];
else
    
AvgAUCgoAll = [AvgAUCgoAll AvgAUCgo];
AvgAUCnogoAll = [AvgAUCnogoAll AvgAUCnogo];
diffAUCavgAll = [diffAUCavgAll diffAUCavg];

end
end
figure;
stdshade(diffAUCavgAll',[],'g')
hold on
stdshade(diffAUCavgAll2')
% plotweights(1,:) = mean(AvgAUCgoAll2');
% plotweights(2,:) = std(AvgAUCgoAll2')/sqrt(48);
% plotweights(3,:) = mean(AvgAUCnogoAll2');
% plotweights(4,:) = std(AvgAUCnogoAll2')/sqrt(48);
% plotweights(5,:) = mean(diffAUCavgAll2');
% plotweights(6,:) = std(diffAUCavgAll2')/sqrt(48);