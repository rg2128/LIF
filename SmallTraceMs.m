load('AllTraceDynamicSyn3');
for i = 1:2
    
    a =AllTrace3(:,i,:,:,:);
    formatspec = 'Learn-%d.mat';
    filename = sprintf(formatspec,i);
    save(filename,'a');
    
end



% SomeTrace2(imjk,i,m,j)