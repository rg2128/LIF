%% Distribution of top down signal onto inhibitory population

aa = abs(randn(numel(OverlappedDATGoInd),1));
aaa=randperm(numel(aa(:)));
aa(aaa(1:ceil(0.2*numel(aa(:)))))=0;

figure;
histogram(aa)


aa = (randn(numel(1:500),1));
figure;
histogram(aa)

aa =1+ (randn(numel(1:500),1));
%              aa<0;
figure;
histogram(aa)


aa = 10+(randn(numel(1:500),1));
%              aa(aa<0) = 0.5+aa;
            aaa=randperm(numel(aa(:)));
            aa(aaa(1:ceil(0.2*numel(aa(:)))))=0;
            
            figure;
h=histogram(aa);
hold on
plot(conv(h.BinEdges, [0.5 0.5], 'valid'), h.BinCounts)


for i =1:20
    aa = randn(numel(1:500),1);
    aa1=aa;
%     Allaas (:,i) = aa;
    aaa = randperm(numel(aa(:)));
    aa(1:300) = aa(1:300) + 5 + randn(numel(1:300),1);
    aa1 = aa(1:300);
    aa(aaa(1:ceil(0.2*numel(aa(:))))) = 0;
     figure;
    histogram(aa)
    figure;
    histogram(aa(1:300))
end