%% imjk,i,m,j,k
%% JEI = J2(1:500,501:1000);
%%   SomeTrace2(imjk,i,j,m) = {[Odor1 Odor2]};
GoOnly = [];
NoGoOnly = [];
deltaW = [];
diffwG = [];
diffwNG = [];
diffw = [];
% aaaa=[16,15,13,12,10,9,8];
aaaa=1:16;
j=1;
for k = aaaa
     for m = 1
    %     for j =1
    DATaa= AllTrace(1,1,j,m,k);
    DATaa = DATaa{1,1};
    
    aa= AllTrace2(1,1,j,m,k);
    Weightei = AllTrace4(1,1:100,j,m,k);
    goAll = [];
    nogoAll = [];
    go=find(aa{1,1}(1:50)>0);
    nogo=find(aa{1,1}(51:100)>0);
    
    alldatgo = [];
    for ii =go
            [x,~] = (find(DATaa == ii));
            alldatgo = [alldatgo; x];
    end
    
    alldatnogo = [];
    for ii =nogo
            [x,~] = (find(DATaa == ii));
            alldatnogo = [alldatnogo; x];
    end    
    DATnumber(k,:) = [length(alldatgo) length(alldatnogo)];
    a2 = tabulate(alldatnogo);
    a1 = tabulate(alldatgo);
    a1 = tabulate(a1(:,2));
    a2 = tabulate(a2(:,2));
    connection(k,:) = [a1(end,1)' a2(end,1)'];
    for q=1:length(go)
        goAll(q,:) = go(q)*10-9:go(q)*10;
    end
    for q=1:length(nogo)
        nogoAll(q,:) = nogo(q)*10-9:nogo(q)*10;
    end
    count1 = 1;
    

    for w = 1:99
%         diffw = Weightei{1,w+1}(1:500,501:1000)-Weightei{1,1}(1:500,501:1000);
        
        diffw = Weightei{1,w+1}(501:1000,1:500)-Weightei{1,1}(501:1000,1:500);
        deltaW(k,count1) = sum(sum(diffw(:,nogoAll(:))))-sum(sum(diffw(:,goAll(:))));
        NoGoOnly(k,count1) = mean(mean(diffw(:,nogoAll(:))));
        GoOnly(k,count1) = mean(mean(diffw(:,goAll(:))));
        count1 = count1+1;
    end
    
    
    
    count1 = 1;
    for w = 2:2:98
%         diffw = Weightei{1,w+1}(501:1000,1:500)-Weightei{1,1}(501:1000,1:500);        
%         diffwG(k,count1) = sum(sum(Weightei{1,w+1}(501:1000,1:500)-Weightei{1,w}(501:1000,1:500)));
%         diffwNG(k,count1) = sum(sum(Weightei{1,w+2}(501:1000,1:500)-Weightei{1,w+1}(501:1000,1:500)));
        diffwG(k,count1) = mean(mean(Weightei{1,w+1}(501:1000,goAll(:))-Weightei{1,w}(501:1000,goAll(:))));
        diffwNG(k,count1) = mean(mean(Weightei{1,w+2}(501:1000,nogoAll(:))-Weightei{1,w+1}(501:1000,nogoAll(:))));
        count1 = count1+1;
    end
    
    
     end
end
% J3 = Weightei{1,4}(501:1000,1:500)-Weightei{1,1}(501:1000,1:500);
% J3 = Weightei{1,1}(501:1000,1:500);
% for ii =1:50
%     J3(ii,:) = mean(Weightei{1,49}(501+ii*10-10:500+ii*10,:));
% end
% for ii =1:50
%     JJ3(:,ii) = mean(J3(:,1+ii*10-10:ii*10),2);
% end
% figure;
% imagesc(JJ3(:,[go nogo]))
% caxis([0 0.05])
% figure;
% imagesc(JJ3)
% caxis([0 0.1])

% figure;
% % % subplot(2,1,1);
% colormap(parula); %xy=J; fun_colormapLim;
% imagesc(J3);
% aux.figset(gca,'neurons','neurons','weights',10);
% colorbar;
% 
% figure;
% stdshade(diffwG,[],'g')
% hold on
% stdshade(diffwNG)
%  xlim([0 50])
% 
% figure;
% stdshade(deltaW,[],'g')
% alld = median(deltaW);
% for k=aaaa
% figure;
% plot(GoOnly(k,:),'b','LineWidth', 2)
% hold on
% plot(NoGoOnly(k,:),'r','LineWidth' , 2)
% 
% figure;
% plot(diffwG(k,:),'b','LineWidth', 2)
% hold on
% plot(diffwNG(k,:),'r','LineWidth' , 2)
% % formatSpec = 'WeightChange-%d.jpg';
% %                 file_name = sprintf(formatSpec,k);
% %                 saveas(gcf,fullfile('Simulations',file_name),'jpg');
% % xlim([1 5])
% end
% figure;
% plot(alld)
%         figure;
%         imagesc(Go2)
%         caxis([-0.5 0.5])
%     end
% end
figure;
stdshade(GoOnly(:,1:2:99),[],'g')
hold on
stdshade(NoGoOnly(:,2:2:99))
%   xlim([0 40])
GoOnlyAll = [GoOnlyAll; GoOnly];
NoGoOnlyAll = [NoGoOnlyAll; NoGoOnly];


plotweights(1,:) = mean(GoOnlyAll(:,3:2:99));
plotweights(2,:) = std(GoOnlyAll(:,3:2:99))/sqrt(48);
plotweights(3,:) = mean(NoGoOnlyAll(:,2:2:99));
plotweights(4,:) = std(NoGoOnlyAll(:,2:2:99))/sqrt(48);

% Diffweight  = GoWeight{1} - GoWeight{5};
% deltaW(m,j) = sum(sum(Diffweight([21:30 41:50 71:80 171:180 281:290], :)))-sum(sum(Diffweight([11:20 241:250 291:300 371:380 471:480], :)));
%     end
% end


% sum(sum(Diffweight([91:100 131:140 231:240 391:410], :)))
% sum(sum(Diffweight([61:70 161:170 271:280 331:340 351:360], :)))
