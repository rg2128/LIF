load('AllTraceDynamicSyn4');
for i = 1:16
    
    a =AllTrace4(:,:,:,:,:);
    formatspec = 'Learn2-%d.mat';
    filename = sprintf(formatspec,i);
    save(filename,'a','-v7.3');
    
end