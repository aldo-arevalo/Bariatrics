function [AUROC] = treePerformanceKFold(...
    CMP,Xtrain,Ytrain,Wtrain,Xtest,Ytest,Wtest)

[Ypredict, testScores] = predict(CMP,Xtest);

[~,~,Ttree,~] = perfcurve(Ytest,testScores(:,end),1);
[~ , AUROC ] = getPerformanceMetrics(Ytest,Ypredict,Ttree);

[ACC,pr,speci,sens,Fs,~,] = confusion_matrix_([Ytest Ypredict]);

fprintf('Accuracy %0.2f \n',ACC);
fprintf('AUROC %0.2f \n', AUROC);
fprintf('Precision %0.2f \n', pr);
fprintf('Specificity %0.2f \n', speci);
fprintf('Sensitivity %0.2f \n', sens);
fprintf('F1 score %0.2f \n', Fs);
fprintf('Next fold \n');

end