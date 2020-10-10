function [threshold, threshold_AUC] = getPerformanceMetrics(realtarget,predictedtarget,list_threshold)

res = [];

% Calculate the performance for each threshold.
for f = 1 : size(list_threshold,1)
    
    alfa = list_threshold(f,1);
    
    p = predictedtarget;
    
    % Classification
    p( p < alfa ) = 0; % inferior
    p( p ~= 0 ) = 1; % sup ou igual

    [~,~,~,AUC] = perfcurve(realtarget,p,1);
    
    compare_class_crisp = [realtarget p];
    
    % Evaluation
    [accuracy,precision,specificity,sensitivity,Fscore] = confusion_matrix_(compare_class_crisp);
    
    % matrix res contains the results for each threshold in this order:
    % threshold, accuracy, precision, specificity, sensitivity
    res = [res; list_threshold(f) accuracy precision specificity sensitivity AUC Fscore];
end

% The next measures are not independent from the threshold.
% we select the threshold that minimizes the difference between sensitivity
% and specificity.

[~, I] = nanmin(abs((res(:,4)-res(:,5))));
otherID=find(res(:,5)==res(I,5));
bestI=find(res(:,4)==max(res(otherID,4)) & res(:,5)==res(I,5));
threshold=res(bestI(1),1);
threshold_AUC=res(bestI(1),6);

if isempty(threshold)
    if isnan(AUC) == 0
        [~, I] = nanmin(abs((res(:,4)-res(:,5))));
        otherID=find(res(:,5)==res(I,5));%
        bestI=find(res(:,4)==max(res(otherID,4)) & res(:,5)==res(I,5));
        threshold=res(bestI(1),1);
    else
        threshold = 0;
        AUC = 0;
    end
end
