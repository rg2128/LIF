for i =1:size(AllTrace,2)
    figure;
    plot(AllTrace2{41,1,2,1}{1,1}(1,:))
    hold on
    plot(AllTrace2{41,1,2,1}{1,1}(2,:))
    hold on
    plot(AllTrace2{42,1,2,1}{1,1}(1,:))
    hold on
    plot(AllTrace2{42,1,2,1}{1,1}(2,:))
    hold off
end
count = 1;
for i =1:2:size(AllTrace,1)
    a(count,:) = mean(AllTrace{i,1,1}{1,1}(1:2,:),1);
    b(count,:) = mean(AllTrace{i+1,1,1}{1,1}(1:2,:),1);
    %     figure;
    %     plot(a)
    %     hold on
    %     plot(b)
    %     hold off
    count = count+1;
    %     plot(AllTrace{i,1,12}{1,1}(1,:))
    %     hold on
    %     plot(AllTrace{i,1,12}{1,1}(2,:))
    %     hold on
    %     plot(AllTrace{i+1,1,12}{1,1}(1,:))
    %     hold on
    %     plot(AllTrace{i+1,1,12}{1,1}(2,:))
    %     hold off
end
figure;
stdshade(a, 0.1,'r')
hold on
stdshade(b, 0.1,'y')


%
% AllTrace2 = AllTrace21;
% AllTrace21 = AllTrace2;
% AllTrace2 = AllTrace2(1,:,:,:,:);
% AllTrace2 = permute(AllTrace2, [2 3 4 5 1]);

for learn= 1:30
    
    formatspec = 'Learn-%d.mat';
    filename = sprintf(formatspec,learn);
    AllTrace2 = load(filename);
    AllTrace2 = AllTrace2.a;
    distances=[];
    ActivationFrame=[];
    Onset =[];
    PeakAmplitude=[];
    activeglomeruli=[];
    ConsistentGlom = [];
    meanOnset = [];
    difOnset = [];
    AllTrace2 = permute(AllTrace2, [1 3 4 5 2]);

    for h =1:16
        for i =1:2
            for m2 =1:4
                for j=1:size(AllTrace2,1)               %For the number of trials
                    for k=1:size(AllTrace2{j,m2,i,h}{1,1},1)
                        distances=[];
                        %                 if i ==1
                        %                     sdval= std(AllTrace2{j,m2,i+1,h}{1,1}(k,1:60));
                        %                 elseif i == 2
                        sdval= std(AllTrace2{j,m2,i,h}{1,1}(k,1:60));
                        %                 end
                        for m =1:30
                            if AllTrace2{j,m2,i,h}{1,1}(k,101+m)> 4*sdval && AllTrace2{j,m2,i,h}{1,1}(k,101+m+1)> 4*sdval && AllTrace2{j,m2,i,h}{1,1}(k,101+m+2)> 4*sdval %&& Trace1(35+k+4,j,i)> 3*sdval
                                Onset(h,i,m2,j,k) = m;
                                break;
                            end
                        end
                        
                    end
                    meanOnset(h,i,m2,j) = nanmean(nonzeros(Onset(h,i,m2,j,:)));
                end
            end
        end
    end
%     for h =1:16
%         for i=1:2
%             for m2 = 1:2
%                 count = 1;
%                 for j =1:2:size(AllTrace2,1)-1
%                     difOnset(count,i,h,m2) = -meanOnset(h,i,m2,j) + meanOnset(h,i,m2,j+1);
%                     count = count+1;
%                 end
%             end
%         end
%     end
    
    
%     for i =1:16
%         for j =1:2
%             for m2 =1:4
%                 countp(i,j,m2) = length(nonzeros(difOnset(:,j,i,m2)>0));
%             end
%         end
%     end
    
%     
%     a = mean(countp(:,1,:),1)
%     b = mean(countp(:,2,:),1);
%     
%     for i =1:20
%         [h,p(i)] = ttest(countp(:,1,i),countp(:,2,i));
%     end

       TrialMean = mean(meanOnset,4);
        
        LearnTrials(learn,:,:,:) = TrialMean;
end


    
    
    
    for i = 1:4
        a1 = meanOnset(:,1,i,1:2:40);
        b1 = meanOnset(:,1,i,2:2:40);
        
        a1 = permute(a1,[4 1 2 3]);
        a1 = nanmean(a1,1);
        
        b1 = permute(b1,[4 1 2 3]);
        b1 = nanmean(b1,1);
        c1(i,:) = b1-a1;
        c1 = b1-a1;
        
        
        a2 = meanOnset(:,2,i,1:2:40);
        b2 = meanOnset(:,2,i,2:2:40);
        
        a2 = permute(a2,[4 1 2 3]);
        a2 = nanmean(a2,1);
        b2 = permute(b2,[4 1 2 3]);
        b2 = nanmean(b2,1);
        c2(i,:) = b2-a2;
        c2 = b2-a2;
        
        for j =1:length(c1(i,:))
            if c2(i,j)>0
                c3(i,j) = c2(i,j)/c1(i,j);
            end
        end
    end
    
    
    for i =1:7
        count(i) = 0;
        for j =1:16
            if c3(i,j)<0
                count(i) = count(i) +1;
            elseif c3(i,j)>0 && abs(c3(i,j))>1
                count(i) = count(i) +1;
                
            end
        end
    end
    
    
    for i=1:16
        aaaa = meanOnset(i,1,1,:);
        
        aaaa = squeeze(aaaa);
        aaaa1 = aaaa(1:2:24);
        aaaa2 = aaaa(2:2:24);
        [h,p2(i)] = ttest2(aaaa1,aaaa2);
    end
    for j =1:16
        for i = 1:20
            cc(j,i) = length(nonzeros(c1(j,:,1,i)>0));
        end
    end
    
    meanO1 = meanOnset(:,:,4,:);
    meanTO1 = meanO1(:,2,:,:);
    meanNTO1 = meanO1(:,1,:,:);
    for i =1: size(meanTO1,1)
        %     for j =1:2:24
        figure;
        adiff = meanTO1(i,:,1,2:2:50)-meanTO1(i,:,1,1:2:50);
        scatter(1:25,adiff);
        %     hold on
        %     scatter(1:25,meanNTO1(i,:,1,:));
        %         adiff(:,j) = meanNTO1(:,:,1,j)-meanNTO1(:,:,1,j+1);
        
    end
    
    for i =1: size(meanTO1,1)
        figure;
        
        scatter(1:25,meanTO1(i,:,1,1:2:50));
    end
    
    
    
    
    %       distances=[];
    % ActivationFrame=[];
    % Onset =[];
    % PeakAmplitude=[];
    % activeglomeruli=[];
    % ConsistentGlom = [];
    % for h =1:16
    %     for i =1:2
    %
    %         for j=1:size(AllTrace,1)               %For the number of trials
    %             for k=1:size(AllTrace{j,i,h}{1,1},1)
    %                 distances=[];
    %                 sdval= std(AllTrace{j,i,h}{1,1}(k,1:25));
    %                 for m =1:10
    %                     if AllTrace{j,i,h}{1,1}(k,30+m)> 3*sdval && AllTrace{j,i,h}{1,1}(k,30+m+1)> 3*sdval && AllTrace{j,i,h}{1,1}(k,30+m+2)> 3*sdval %&& Trace1(35+k+4,j,i)> 3*sdval
    %                         Onset(h,i,j,k) = m;
    %                     break;
    %                 end
    %                 end
    %
    %             end
    %             meanOnset(h,i,j) = mean(nonzeros(Onset(h,i,j,:)));
    %         end
    %
    %     end
    % end
    % for h =1:16
    %     for i=1:2
    %
    %         count = 1;
    %         for j =1:2:size(AllTrace,1)
    %             difOnset(count,i,h) = meanOnset(h,i,j) - meanOnset(h,i,j+1);
    %             count = count+1;
    %         end
    %         end
    %
    % end
    %
    %
    % for i =1:16
    %     for j =1:2
    %
    %         countp(i,j) = length(nonzeros(difOnset(:,j,i)>0));
    %
    %     end
    % end
    %
    %
    % a = mean(countp(:,1),1)
    % b = mean(countp(:,2),1)
    %
    
    % [h,p] = ttest(countp(:,1),countp(:,2));
