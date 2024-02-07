for learn= 1:3
    Thresh = 10;
    formatspec = 'Learn-%d.mat';
    filename = sprintf(formatspec,learn);
    AllTrace3 = load(filename);
    AllTrace3 = AllTrace3.a;
    AllClusterTime = [];
    AllClusterTime2 = [];
    AllTrace3 = permute(AllTrace3, [1 3 4 5 2]);

    for h =1:16
        for i =1:2
            for m2 =1:10
                for j=1:size(AllTrace3,1)               %For the number of trials
                    temp2 = AllTrace3{j,m2,i,h};
                    temp = temp2;
                    start = find(temp(:,1)> 0, 1 );
                    stop = find(temp(:,1) < 0.5 , 1, 'last' );
                    for mj =start:stop
                        temp(mj,2) = ceil(temp(mj,2)/10);
                    end        
                    FindActiveClusters = tabulate(temp(start:stop,2));
                    miCount = 1;
                    ActiveClusters = [];
                    for mi = 1:size(FindActiveClusters,1)
                        if FindActiveClusters(mi,2) >Thresh-1
                            ActiveClusters(miCount) = FindActiveClusters(mi,1);
                            miCount =miCount +1;
                        end
                    end
                    IndClusterTime = [];
                    for mk =1: length(ActiveClusters)
                        tp = start-1 +find(temp(start:stop,2) == ActiveClusters(mk));
                        IndClusterTime(mk) = mean(temp(tp(1:Thresh),1));
                    end
                    AllClusterTime(j) = mean(IndClusterTime);
                end
                AllClusterTime2(h,i,m2) = mean(AllClusterTime);
            end
        end
    end
    LearnTime(learn,:,:,:) = AllClusterTime2;
end    
                           
                                
Learn2 = mean(LearnTime,2);    
Learn2Go = Learn2(1,:,:,:);
Learn2NoGo = Learn2(2,:,:,:);

% Learn2Go = LearnTime(3,:,:,:,:);
% Learn2NoGo = LearnTime(2,:,:,:,:);

% Learn2Go = Learn2(1,:,:,:);
% Learn2NoGo = Learn2(2,:,:,:);
% figure;
% for i =5:2:20
% DiffLearn = Learn2NoGo(:,:,2,i)-Learn2Go(:,:,2,i);
% scatter(i,sum(DiffLearn))
% hold on 
% end
figure;
% DiffLearn=[];
for i =1:10
    DiffLearn=[];

DiffLearn = Learn2NoGo(:,:,1,i)-Learn2Go(:,:,1,i);
scatter(i,mean(DiffLearn))
% hold on 

% NTGNG(i,:) =  Learn2NoGo(:,:,1,i)-Learn2Go(:,:,1,i);
% plot(NTGNG(i,:))
hold on
end

DiffLearn=[];

figure;
for i =1:10
%     DiffLearn=[];
% 
DiffLearn = Learn2NoGo(:,:,2,i)-Learn2Go(:,:,2,i);
scatter(i,mean(DiffLearn))
hold on 
% TGNG(i,:) =  Learn2NoGo(:,:,2,i)-Learn2Go(:,:,2,i);
% plot(TGNG(i,:))
% hold on
end

% plot(TGNG(1,:))

DiffLearn=[];

figure;
for i =1:10
    DiffLearn=[];

DiffLearn = Learn2Go(:,:,1,i)-Learn2Go(:,:,2,i);
scatter(i,mean(DiffLearn))
hold on 

% TNTG(i,:) =  Learn2Go(:,:,1,i)-Learn2Go(:,:,2,i);;
% plot(TNTG(i,:))
% hold on

end

% figure;
% for i =1:20
%     DiffLearn=[];
% 
% DiffLearn = Learn2NoGo(:,:,1,i)-Learn2NoGo(:,:,2,i);
% scatter(i,mean(DiffLearn))
% hold on 
% end
% 
% 
% figure;
% for i =1:20
%     DiffLearn=[];
% 
% DiffLearn =  (Learn2NoGo(:,:,2,i)-Learn2Go(:,:,2,i))- (Learn2NoGo(:,:,1,i)-Learn2Go(:,:,1,i));
% scatter(i,mean(DiffLearn))
% hold on 
% end



% figure;
% plot(Learn2Go(:,:,2,2))
% hold on 
% plot(Learn2Go(:,:,1,2))
% %                       
%                         
                        