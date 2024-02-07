%% function-1
% tau_cue=[0.5,5];
% temp2 = [];
% % p = @(t)(-1/(tau_cue(2)-tau_cue(1)))*(exp(-t/tau_cue(2))-exp(-t/tau_cue(1)));
% p = @(t)(-1) % or -t
% Tseq = -1:0.0001:1;
% p=@(t)p(t+0.5);
% temp2 = zeros(20000,1);
% for t = 5000:10000
%      temp2(t) = p(Tseq(t));
% end
% % hold on
% plot(temp2)


%%Function-2
temp2 = [];
p = @(t)(t-0.5)*(t)
Tseq = -1:0.0001:1;
p=@(t)p(t+0.5);
for t = 5000:6000
     temp2(t) = p(Tseq(t));
end
% hold on
plot(temp2)



%%Function-3

temp2 = zeros(20000,1);
start = -0.5;
stop = 0.2;
% p = @(t)(t)
p= @(t)(-0.5-t)
Tseq = -1:0.0001:1;
% p=@(t)p(t+0.5);
for t = 5000:12000
     temp2(t) = p(Tseq(t));
end
% hold on
plot(temp2)


%%Function-4

temp2 = zeros(629,1);
temp2 = [];
start = -0.5;
stop = -0.4;
Tseq = -1:0.0001:1;
% for t = 1:4999+(200*rr)
%      temp2(t) = p(3.14*Tseq(t));
% end
% for t = 5000+(200*rr):14999
%      temp2(t) = p(6.28*Tseq(t));
% end
% for t = 15000:20000
%      temp2(t) = p(3.14*Tseq(t));
% end
% 
% for t = 1:20000
%      temp2(t) = p(6.28*Tseq(t));
% end
figure;
for i =1:30
    figure;
rr2 = 0;
rr = randi(10,1,1);

p = @(t)(sin(rr-5*t));

while rr2<0.5 && rr2>-0.5
    rr2 = (-10 + (20)*rand(1,1))/10;
end

for t =1:9999
     temp2(t) = p(rr2*3.14*Tseq(t));
end

for t = 10000:14999
     temp2(t) = p(rr2*9.42*Tseq(t));
end

for t = 15000:20000
     temp2(t) = p(rr2*3.14*Tseq(t));
end

plot(temp2)
% hold on
end

temp2 = [];
start = -0.5;
stop = -0.4;
% p = @(t)(cos(5*t));
Tseq = -1:0.0001:1;
figure;
% for i =1:30
rr2 = 0;
rr = randi(10,1,1);

p = @(t)(sin(rr-5*t));

while rr2<0.5 && rr2>-0.5
    rr2 = (-10 + (20)*rand(1,1))/10;
end

for t =1:4999
     temp2(t) = p(rr2*3.14*Tseq(t));
end

for t = 5000:9999
     temp2(t) = p(rr2*2*3.14*Tseq(t));
end

for t = 10000:14999
     temp2(t) = p(rr2*9.41*Tseq(t));
end

for t = 15000:20000
     temp2(t) = p(rr2*3.14*Tseq(t));
end

plot(temp2)
hold on
% end

prune = 0.5;
% bgain = 4;
tau_i_Thresh = 0.0207;
tau_i_tc = 3;
tau_i = 0.02;
for Odor = 2:2:10
aaaa(Odor) = tau_i_Thresh - (tau_i_Thresh-tau_i)*exp((2-Odor)/tau_i_tc)
end

