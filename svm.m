function Output3 = svm(score2,c,Index,Name2,d2)
UnsortedTrials=[];
UnsortedTrials=score2(:,1:3);
data1 = UnsortedTrials([Index.Go],:);
data2 = UnsortedTrials([Index.NoGo],:);
data3= [data1; data2];
theclass = ones(20,1);
theclass(1:10) = -1;

% data1 = UnsortedTrials([Index.O1],:);
% data2 = UnsortedTrials([Index.O2],:);
% data3= [data1; data2];
% theclass = ones(78,1);
% theclass(1:39) = -1;

cl=[];
cl = fitcsvm(data3, theclass, 'BoxConstraint', 100, 'KernelFunction', 'linear', 'Standardize',true);
% cl = fitcsvm(data3,theclass,'Standardize',true,'KernelFunction','RBF','KernelScale','auto');
sv = cl.SupportVectors;
% figure
% scatter3(score1(:,1), score1(:,2),  score1(:,3) , 50, c, 'filled')
% hold on
% plot3(sv(:,1),sv(:,2),sv(:,3),'ko','MarkerSize',10)
% legend('versicolor','virginica','Support Vector')
% hold off
CVSVMModel = crossval(cl);
classLoss = kfoldLoss(CVSVMModel);
d = d2;
x1Grid =[];
x2Grid =[];
x3Grid =[];

[x1Grid,x2Grid,x3Grid] = meshgrid(min(data3(:,1)):d:max(data3(:,1)),min(data3(:,2)):d:max(data3(:,2)), min(data3(:,3)):d:max(data3(:,3)));
xGrid = [x1Grid(:),x2Grid(:),x3Grid(:)];
[~,scores] = predict(cl,xGrid);
zGrid = reshape(scores(:,2), size(x1Grid));
% figure,
% scatter3(score2(:,1), score2(:,2),  score2(:,3) , 50, c, 'filled');
% xlabel 'PC1', ylabel 'PC2', zlabel 'PC3'
% hold on
% plot3(data3(cl.IsSupportVector,1),data3(cl.IsSupportVector,2), data3(cl.IsSupportVector,3), 'ko','MarkerSize',10);
% [faces,verts,~] = isosurface(x1Grid, x2Grid, x3Grid, zGrid, 0, x1Grid);
% patch('Vertices', verts, 'Faces', faces, 'FaceColor','k','edgecolor','none', 'FaceAlpha', 0.2);
% grid on
% box on
% view(3)
% hold off
% Name3 = strcat(Namex,'svm');
% savefig(Name3);
Output3= classLoss;
end




% figure;
% h(1:2) =  scatter3(score1(:,1), score1(:,2),  score1(:,3) , 50, c, 'filled');
% hold on
% h(3) = plot3(data3(cl.IsSupportVector,1),data3(cl.IsSupportVector,2),data3(cl.IsSupportVector,3),'ko');
% contour3(x1Grid, x2Grid, x2Grid(1:2), reshape(scores(:,2),size(x1Grid)),[0 0],'k');
% legend(h,{'-1','+1','Support Vectors'});
% axis equal
% hold off