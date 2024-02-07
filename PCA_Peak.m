function Output2 = PCA_Peak(PeakAmplitude,c,Index,Name)
[coeff2,score2,latent2,tsquared2,explained2,mu2] = pca(PeakAmplitude);
% figure;
% scatter3(score2(:,1), score2(:,2),  score2(:,3) , 50, c, 'filled')
d2=0.5;
%scatter(score1(:,1), score1(:,2), 50, c, 'filled')
% axis equal
% xlabel('1st Principal Component')
% ylabel('2nd Principal Component')
% zlabel('3rd Principal Component')
Name2 = strcat(Name,'PCA_Peak');
% savefig(Name2);
Output3 = svm(score2,c,Index,Name2,d2);
Output2 = {explained2, Output3};
end
