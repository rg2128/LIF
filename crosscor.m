AllTrace = AllTrace3;
AllTrace = permute(AllTrace, [1 2 4 3 5]);
tempmaxc=0;
for j = 1:2
    for i = 1:16 %k
        for mm =1:100 
            allsp = AllTrace(:,:,:,j,i);
            allsp = allsp(:,:,mm);
            for k =1:length(allsp)
                for od = 1:2
                time = allsp{k,od}(:,1);
                spikes = allsp{k,od}(:,2);
                tsp = allsp{k,od};
                count = zeros(100,length(spikes));
                temp = -0.5;
                p = 1;
                for m =1: length(spikes)
                    if time(m,1) > -0.5 && time(m,1) <0.1
                        for n=1:100
                            if spikes(m,1)< (n*10) + 1 && spikes(m,1)> (n-1)*10
                                for q =1:20
                                    if abs(time(m,1) - temp) < 0.005*q && abs(time(m,1) - temp) > 0.005*(q-1)
                                        p = q+p-1;
                                        break;
                                    end
                                end
                                count(n,p) = count(n,p) + 1;
                            end
                        end
                    end
                    temp = time(m,1);
                end
                glom = count(1:50,:);
                dat = count(51:100,:);
                glom = glom ~= 0;
                 dat = dat ~= 0;
                glom = sum(glom,1);
                dat = sum(dat,1);
%                 [c,lags] = xcov(glom');
 [c,lags] = xcov(dat','unbiased');
%  [c,lags] = xcov(dat');
%  stem(c,lags)
%  maxlag = 10;
% [c,lags] = xcov(dat,maxlag,'normalized');
%                 [c,lags] = xcov(glom,dat,'scaleopt','normalized');
%                 stem(c,lags)
% maxc = max(c);
% if maxc>tempmaxc
%     tempmaxc = maxc;
% end

                cgng(od) = max(c);
                glome(od) = sum(glom);
                date(od) = sum(dat);
%                             stem(lags,c)
                %             lags(1,find(c==max(c)));
                end
                cc(:,k) = cgng;
                tglom(:,k) = glome;
                tdat(:,k) = date;
                
            end
            cccG (mm,:) = cc(1,:);
            cccNG(mm,:) = cc(2,:);
            Gglom(mm,:) = tglom(1,:);
            NGglom(mm,:) = tglom(2,:);
            Gdat(mm,:) = tdat(1,:);
            NGdat(mm,:) = tdat(2,:);
        end
        if j ==1
            NTG(:,:,i) = cccG;
            NTNG(:,:,i) = cccNG;
            NTGglom(:,:,i) =Gglom;
            NTNGglom(:,:,i) = NGglom;
            NTGdat(:,:,i) =Gdat;
            NTNGdat(:,:,i) =NGdat;
        else
             TG(:,:,i) = cccG;
            TNG(:,:,i) = cccNG;
             TGglom(:,:,i) =Gglom;
            TNGglom(:,:,i) = NGglom;
            TGdat(:,:,i) =Gdat;
            TNGdat(:,:,i) =NGdat;
        end
        
    end
end
%   GoT = T(1,1:2:40,1:16); 
%         NoGoT = T(1,2:2:40,1:16);

AvgTG = mean(TG,2);
AvgTG = permute(AvgTG,[1 3 2]);
AvgTNG = mean(TNG,2);
AvgTNG = permute(AvgTNG,[1 3 2]);
AvgNTG = mean(NTG,2);
AvgNTG = permute(AvgNTG,[1 3 2]);
AvgNTNG = mean(NTNG,2);
AvgNTNG = permute(AvgNTNG,[1 3 2]);

AvgTGglom = mean(TGglom,2);
AvgTGglom = permute(AvgTGglom,[1 3 2]);
AvgTNGglom = mean(TNGglom,2);
AvgTNGglom = permute(AvgTNGglom,[1 3 2]);
AvgNTGglom = mean(NTGglom,2);
AvgNTGglom = permute(AvgNTGglom,[1 3 2]);
AvgNTNGglom = mean(NTNGglom,2);
AvgNTNGglom = permute(AvgNTNGglom,[1 3 2]);

AvgTGdat = mean(TGdat,2);
AvgTGdat = permute(AvgTGdat,[1 3 2]);
AvgTNGdat = mean(TNGdat,2);
AvgTNGdat = permute(AvgTNGdat,[1 3 2]);
AvgNTGdat = mean(NTGdat,2);
AvgNTGdat = permute(AvgNTGdat,[1 3 2]);
AvgNTNGdat = mean(NTNGdat,2);
AvgNTNGdat = permute(AvgNTNGdat,[1 3 2]);
% 
% 
% for bgf = 1:25
% figure;
% scatter(ones(16,1),AvgTG(bgf,:))
% hold on 
% scatter(2*ones(16,1),AvgTNG(bgf,:))
% hold on 
% scatter(3*ones(16,1),AvgNTG(bgf,:))
% hold on 
% scatter(4*ones(16,1),AvgNTNG(bgf,:))
% hold off
% end
% 
% for bgf = 1:25
% figure;
% scatter(ones(16,1),AvgTGdat(bgf,:))
% hold on 
% scatter(2*ones(16,1),AvgTNGdat(bgf,:))
% hold on 
% scatter(3*ones(16,1),AvgNTGdat(bgf,:))
% hold on 
% scatter(4*ones(16,1),AvgNTNGdat(bgf,:))
% hold off
% end

AllStimAvgTG = mean(AvgTG,2);
AllStimAvgTG = reshape(AllStimAvgTG,10,10);

AllStimAvgNTG = mean(AvgNTG,2);
AllStimAvgNTG = reshape(AllStimAvgNTG,10,10);

clims = [1 4]
figure;
imagesc(AllStimAvgTG,clims)
figure;
imagesc(AllStimAvgNTG,clims)



figure;
imagesc(AvgTG([ 1 6 11 16 21],:))
figure;
imagesc(AvgTNG)
for bgf = 1:25
figure;
scatter(ones(16,1),AvgTGglom(bgf,:))
hold on 
scatter(2*ones(16,1),AvgTNGglom(bgf,:))
hold on 
scatter(3*ones(16,1),AvgNTGglom(bgf,:))
hold on 
scatter(4*ones(16,1),AvgNTNGglom(bgf,:))
hold off
end