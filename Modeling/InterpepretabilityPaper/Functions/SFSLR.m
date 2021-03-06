function [selected_set, AUC_selected_set, results] = ...
    SFSLR(data_train,data_val,choice,responseName,predictorNames,UF)

% Function to run Sequential Forward Selection (SFS) for Logistic
% Regression models
%
% INPUTS:
%
% OUTPUTS:
%
%
% Dependencies: M1_Modeling > runSFS
%
% Author: Aldo Arévalo
% Date: 10/11/2020

nr_features=length(predictorNames);
stop=0;
count=1;
set=ones(1,nr_features);
results=[];
id_best=[];
previous_best_perf=0;
count_min=1;
stop_criterion=choice(1);
plot_results=choice(2);
previous_best_set=[];

if stop_criterion==0  % don't stop when there is no improvement
    while stop~=1
        for i=1:length(set)
            if set(i)~=0
                
                if isempty(id_best)
                    inmodel=i;
                    feature_subset = predictorNames{inmodel};
                    var_names = horzcat({feature_subset}, {responseName});
                else
                    inmodel=[find(best==1) i];
                    feature_subset = predictorNames(inmodel);
                    var_names = feature_subset;
                    var_names{end+1} = responseName;
                end
                
                if choice(4) == 1
                    % BALANCE DATASET
                    Class_0 = find(data_train.TWL == 0);
                    Class_1 = find(data_train.TWL == 1);

                    if (length(Class_1) < length(Class_0)) == 1
                        trainInd = vertcat(Class_1, ...
                            randsample(Class_0,round(size(Class_1,1)*UF),...
                            false));
                    else
                        trainInd = vertcat(Class_0, ...
                            randsample(Class_1,round(size(Class_0,1)*UF),...
                            false));
                    end
                    
                    % CREATE MODEL
                    mdl = fitglm(data_train(trainInd,var_names), ...
                        'Distribution', 'binomial', ...
                        'link', 'logit',...
                        'ResponseVar',responseName);
                else
                   % CREATE MODEL
                    mdl = fitglm(data_train(:,var_names), ...
                        'Distribution', 'binomial', ...
                        'link', 'logit',...
                        'ResponseVar',responseName); 
                end
                
                % TEST MODEL
                yfit = predict(mdl,data_val(:,var_names));
                
                [~,~,Thresholds,~] = perfcurve(data_val.(responseName), yfit,1);  
                [~ , AUC ] = getPerformanceMetrics(data_val.(responseName),...
                    yfit,Thresholds);
                
                res(count,inmodel)=1;
                res(count,nr_features+1)=AUC;
                count=count+1;
            end
        end
        id_best=find(res(:,nr_features+1)==max(res(:,nr_features+1)));
        if length(id_best)>1
            idd=randperm(length(id_best));
            id_best=id_best(idd(1));
            disp('Attention, more than one set giving the same result. In order to proceed a random set was selected')
        end
        best=res(id_best,1:end-1);
        %     id_best=id_best(1);
        set(1,find(best==1))=0;
        results=[results;res];
        if res(id_best,end)<previous_best_perf
            localMinPerf(count_min)=previous_best_perf;
            count_min=count_min+1;
            stop=0;
        end
        if size(res,1)==1
            stop=1;
        end
        previous_best_perf=res(id_best,end);
        clear res
        count=1;
    end
    id=find(results(:,end)==max(results(:,end)));
    if length(id)>1
        idd=randperm(length(id_best));
        id=id(idd(1));
        disp('Attention, more than one set giving the same result. In order to proceed a random set was selected')
    end
    AUC_selected_set=results(id(1),end);
    selected_set=find(results(id(1),1:end-1)==1);
    
elseif stop_criterion==1 % stop when there is no improvement
    while stop~=1
        for i=1:length(set)
            if set(i)~=0
                if isempty(id_best)
                    inmodel=i;
                    feature_subset = predictorNames{inmodel};
                    var_names = horzcat({feature_subset}, {responseName});
                else
                    try
                        inmodel=[find(best==1) i];
                        feature_subset = predictorNames(inmodel);
                        var_names = feature_subset;
                        var_names{end+1} = responseName;
                    catch
                        inmodel=[find(best==1) i];
                        feature_subset = predictorNames(inmodel);
                        var_names = feature_subset;
                        var_names{end+1} = responseName;
                    end
                end
                
                if choice(4) == 1
                    % BALANCE DATASET
                    Class_0 = find(data_train.TWL == 0);
                    Class_1 = find(data_train.TWL == 1);

                    if (length(Class_1) < length(Class_0)) == 1
                        trainInd = vertcat(Class_1, ...
                            randsample(Class_0,round(size(Class_1,1)*UF),false));
                    else
                        trainInd = vertcat(Class_0, ...
                            randsample(Class_1,round(size(Class_0,1)*UF),false));
                    end
                    
                    % CREATE MODEL
                    mdl = fitglm(data_train(trainInd,var_names), ...
                        'Distribution', 'binomial', ...
                        'link', 'logit',...
                        'ResponseVar',responseName);
                else
                   % CREATE MODEL
                    mdl = fitglm(data_train(:,var_names), ...
                        'Distribution', 'binomial', ...
                        'link', 'logit',...
                        'ResponseVar',responseName); 
                end
                
                % TEST MODEL
                yfit = predict(mdl,data_val(:,var_names));
                
                [~,~,Thresholds,~] = perfcurve(data_val.(responseName), yfit,1);  
                [~ , AUC ] = getPerformanceMetrics(data_val.(responseName),...
                    yfit,Thresholds);
                
                res(count,inmodel)=1;
                res(count,nr_features+1)=AUC;
                count=count+1;
            end
        end
        id_best=find(res(:,nr_features+1)==max(res(:,nr_features+1)));
        if length(id_best)>1
            idd=randperm(length(id_best));
            id_best=id_best(idd(1));
            disp('Attention, more than one set giving the same result. In order to proceed a random set was selected')
        end
        best=res(id_best,1:end-1);
        set(1,find(best==1))=0;
        results=[results;res];
        e = res(id_best,end)-previous_best_perf;
        %if res(id_best,end)<previous_best_perf || res(id_best,end)>0.95
        if res(id_best,end)<previous_best_perf || e<0.01
            localMinPerf(count_min)=previous_best_perf;
            count_min=count_min+1;       
            if isempty(previous_best_set) % first iteration already giving the best result
                AUC_selected_set=res(id_best,end);
                selected_set=find(best==1);
            else
                AUC_selected_set=previous_best_perf;
                selected_set=previous_best_set;
            end
            stop=1;
        elseif size(res,1)==1
            stop=1;
        end
        previous_best_perf=res(id_best,end);
        previous_best_set=find(best==1);
        clear res
        count=1;
    end
end

for i=1:size(results,1)
    nr_feat(i,1)=sum(results(i,1:end-1));
end

if plot_results==1
    id=find(results(:,end)==localMinPerf(1));
    
    figure
    plot(results(:,end),'Color','k')
    hold on
    scatter(id(1),localMinPerf(1),'filled');
    legend('All','First local maximum','Location','Best')
    xlabel('Iteration')
    ylabel('AUC')
    yyaxis right
    plot(nr_feat)
        hold on
        scatter(id,nr_feat(id),'filled');
    xlabel('Iteration')
    ylabel('Number of features')
    grid on
end