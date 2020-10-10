function [accuracy,precision,specificity,sensitivity,Fscore,npv,confusion_matrix] = confusion_matrix_(compare_class)

% Build the confusion matrix for the expected class vs. predicted class and calculate accuracy, precision, sensitivity and specificity.

tn=[];
tp=[];
fn=[];
fp=[];

for i=1:length(compare_class)
    if((compare_class(i,1)==compare_class(i,2)) && compare_class(i,1)==0)
        tn=[tn;i];
    elseif((compare_class(i,1)==compare_class(i,2)) && compare_class(i,1)==1)
        tp=[tp;i];
    elseif((compare_class(i,1)~=compare_class(i,2)) && compare_class(i,2)==1)
        fp=[fp;i];
    elseif((compare_class(i,1)~=compare_class(i,2)) && compare_class(i,2)==0)
        fn=[fn;i];
    end
end
% [True Negative    False Positive]
% [False Negative   True Positive ]
tn=size(tn,1);
fp=size(fp,1);
fn=size(fn,1);
tp=size(tp,1);

confusion_matrix(1,1)=tn;
confusion_matrix(1,2)=fp;
confusion_matrix(2,1)=fn;
confusion_matrix(2,2)=tp;

accuracy = (tp+tn)/(tp+tn+fp+fn);

precision = tp/(tp+fp);

specificity = tn/(tn+fp);

sensitivity = tp/(tp+fn);

Fscore=(tp+tp)/(tp+tp+fp+fn);

npv=tn/(tn+fn);

end