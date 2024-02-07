 
learningTone = [];
learningNoTone = [];
learningToneLatency = [];
learningNoToneLatency = [];

learningTone =  LearnTrials(:,:,2,1);
 learningNoTone =  LearnTrials(:,:,1,1);
 
 count = 1;
for i =2:2:size(learningTone,1)
    learningToneLatency(count,:) = learningTone(i,:)-learningTone(i-1,:);
    learningNoToneLatency(count,:) = learningNoTone(i,:)-learningNoTone(i-1,:);
    count = count+1;
end
% 
% % P = polyfit(learningToneLatency,1);
%     yfit = P(1)*length(learningToneLatency)+P(2);
learningToneGo = mean(learningTone(1:2:30,:),2);
learningToneNoGo = mean(learningTone(2:2:30,:),2);
learningNoToneGo = mean(learningNoTone(1:2:30,:),2);
learningNoToneNoGo = mean(learningNoTone(2:2:30,:),2);

figure;
plot(learningToneGo, 'o')
hold on
plot(learningToneNoGo, 'o')

figure;
plot(learningNoToneGo, 'o')
hold on
plot(learningNoToneNoGo, 'o')

figure;
learningToneLatency = mean(learningToneLatency,2);
plot(learningToneLatency, 'o')

learningNoToneLatency = mean(learningNoToneLatency,2);
hold on
plot(learningNoToneLatency, 'o')
