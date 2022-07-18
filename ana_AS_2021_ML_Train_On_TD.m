clearvars
load EEG_ML_vars_LMM_MWU % Makes no difference whether we load results from LMM or PCA
% remove principal component columns
T(:,contains(T.Properties.VariableNames,'PC')) = []
rng(8792431) % There is another rng seed before the machine learning, so changing this one probably doesn't matter

%%% Do we want to find different features with TD, or just apply the same
%%% features that we found using AS?

% variable job tells us whether we are looking for features or we are
% already doing the machine learning
%job = 'FindFeatures';
job = 'MachineLearning';

% Choose features based on LMMs or PCA?
%method = 'LMM'
method = 'PCA'

varthresh = 90; % use enough PCs to explain at least this much variance
betathresh = 0.5; % minimum effect size for feature selection with LMMs

%% Do PCA on features

switch job
    case 'FindFeatures'
        
        % Add PCs to table based on wake - sleep TD data
        
        feats = unique(vartype);
        EEGvars = find(~contains(vartype,'IV'));
        Tpc = table(); % new table
        PC_coefs = cell(2,length(feats)-1);
        
        % Wake minus sleep
        wms = T{~logical(T.Sleep2),EEGvars} - T{logical(T.Sleep2),EEGvars};
        Td = table();
        for icol = 1:size(wms,2)
            eval(sprintf('Td.%s = wms(:,icol);',varnames{EEGvars(icol)}));
        end
        
        
        for ift = 2:length(feats)
            % do PCA just on TD data
            [COEFF, ~, LATENT] = pca(Td{TDidx(1:length(TDidx)/2),contains(vartype(EEGvars),feats{ift})});
            VE = LATENT./sum(LATENT).*100;
            PC_coefs{1,ift-1} = feats{ift};
            PC_coefs{2,ift-1} = COEFF;
            PC_coefs{3,ift-1} = VE;
            score = (T{:,contains(vartype,feats{ift})}-...
                nanmean(T{:,contains(vartype,feats{ift})}))*COEFF;
            eval(sprintf('Tpc.PC1_%s_VE%i = score(:,1);',replace(feats{ift},' ','_'),round(VE(1))));
            eval(sprintf('Tpc.PC2_%s_VE%i = score(:,2);',replace(feats{ift},' ','_'),round(VE(2))));
            eval(sprintf('Tpc.PC3_%s_VE%i = score(:,3);',replace(feats{ift},' ','_'),round(VE(3))));
        end
        save('./TrainOnTD/TD_PCAtable.mat', 'Tpc', 'PC_coefs')
    case 'MachineLearning'
        load('./TrainOnTD/TD_PCAtable.mat')
        %load PCtable PC_coefs
        for ift = 1:size(PC_coefs,2)
            cVar = 0; % cumulative variance explained
            pccnt = 0;
            while cVar < varthresh % explain at least x% of the variance
                pccnt = pccnt + 1;
                score = (T{:,contains(vartype,PC_coefs{1,ift})} - nanmean(T{:,contains(vartype,PC_coefs{1,ift})}))*PC_coefs{2,ift};
                eval(sprintf('T.PC%i_%s = score(:,%i);',pccnt,replace(PC_coefs{1,ift},' ','_'),pccnt));
                VaEx = PC_coefs{3,ift};
                cVar = cVar + VaEx(pccnt);
                % print out which variables contribute most to this PC
                v = varnames(contains(vartype,PC_coefs{1,ift}));
                contp = v(PC_coefs{2,ift}(:,1)>=0.5);
                contn = v(PC_coefs{2,ift}(:,1)<=0.5);
                for ict = 1:length(contp)
                    if ict == 1
                        fprintf('     Largest positive contributions to %s PC%i from %s',PC_coefs{1,ift},pccnt,contp{ict})
                    else
                        fprintf(' and %s',contp{ict})
                    end
                end
                fprintf('\n')
                for ict = 1:length(contn)
                    if ict == 1
                        fprintf('     Largest negative contributions to %s PC%i from %s',PC_coefs{1,ift},pccnt,contn{ict})
                    else
                        fprintf(' and %s',contn{ict})
                    end
                end
                fprintf('\n')
            end  
            fprintf('%i total PCs used for %s, Var. explained = %2.2f\n',pccnt,PC_coefs{1,ift},cVar)
        end
end

%%

switch job
    case 'FindFeatures'
        
        %% Table and stats with only TD data
        
        Tas = T(T.Group == categorical({'TD'}),:);
        
        TDsleepF = nan(1,length(varnames));
        TDsleepES = nan(3,length(varnames));
        TDsleepP = nan(1,length(varnames));
        
        predictors = '1 + Sleep + (1|Subject)';
        exclude = 'Sleep Sleep2 Group Group2 Conscious Conscious2 Subject Age CopyNumber';
        
        for ivar = 1:length(varnames)
            if ~contains(exclude,varnames(ivar)) % if the response isn't already a predictor
                mdlstr = sprintf('%s ~  %s',varnames{ivar},predictors)
                lme = fitlme(Tas,mdlstr)
                lme_anova = anova(lme); % anova test on fixed effect terms
                lme_anova.FStat
                
                sIDX = find(strcmpi(lme.Coefficients.Name,'Sleep_Yes'));
                
                
                TDsleepF(ivar) = lme_anova.FStat(2);
                TDsleepP(ivar) = lme_anova.pValue(2);
                
                % Beta coef as effect size
                TDsleepES(1,ivar) = lme.Coefficients.Estimate(2);
                TDsleepES(2,ivar) = lme.Coefficients.Lower(2);
                TDsleepES(3,ivar) = lme.Coefficients.Upper(2);
            end
        end
        %%
        
        switch plotting
            case true
                
                % Main effect of sleep p-values
                
                myfigure2
                stem(-log10(TDsleepP),'filled','linewidth',2)
                xticks(1:length(varnames))
                xticklabels(varnames)
                %yticks([-12:4:12])
                xlim([8.5 length(varnames)+0.5])
                %ylim([0 65])
                xtickangle(45)
                xlabel('EEG features')
                ylabel('TD sleep -log_{10}(p-value)')
                title('TD Biomarkers of conscious state','fontsize',18)
                
                box off
                set(gca,'linewidth',3)
                xAX = get(gca,'XAxis');
                set(xAX,'FontSize', 20)
                set(xAX,'color','k')
                yAX = get(gca,'YAxis');
                set(yAX,'FontSize', 20)
                set(yAX,'color','k')
                set(gca, 'TickDir', 'out')
                set(gcf,'color','w')
                set(gca,'Layer','top')
                axis normal
                print('-dpng','./Figures/SleepLMMMainEffectsPvalue.png')
                print('-dsvg','./Figures/SleepLMMMainEffectsPvalue.svg')
                
                %% Main effect of sleep F-value
                
                myfigure2
                stem(log10(TDsleepF),'filled','linewidth',2)
                xticks(1:length(varnames))
                xticklabels(varnames)
                %yticks([-12:4:12])
                xlim([8.5 length(varnames)+0.5])
                %ylim([-25 25])
                xtickangle(45)
                xlabel('EEG features')
                ylabel('TD sleep log_{10}(F-value)')
                title('TD Biomarkers of conscious state','fontsize',18)
                
                box off
                set(gca,'linewidth',3)
                xAX = get(gca,'XAxis');
                set(xAX,'FontSize', 20)
                set(xAX,'color','k')
                yAX = get(gca,'YAxis');
                set(yAX,'FontSize', 20)
                set(yAX,'color','k')
                set(gca, 'TickDir', 'out')
                set(gcf,'color','w')
                set(gca,'Layer','top')
                axis normal
                print('-dpng','./Figures/SleepLMMMainEffectsFvalue.png')
                print('-dsvg','./Figures/SleepLMMMainEffectsFvalue.svg')
                
                %% Main effect of sleep beta coefs
                
                myfigure2
                stem(TDsleepES(1,:),'filled','linewidth',2)
                pltidx = find(abs(TDsleepES(1,:))>=betathresh);
                scatter(pltidx,TDsleepES(1,pltidx),130,	'k','filled','d')
                scatter(1:length(TDsleepES),TDsleepES(2,:),100,'r','filled','s')
                scatter(1:length(TDsleepES),TDsleepES(3,:),100,'r','filled','s')
                makefighandsome
                plot(1:100,ones(1,100).*betathresh,'k:')
                plot(1:100,ones(1,100).*-betathresh,'k:')
                xticks(1:length(varnames))
                xticklabels(varnames)
                %yticks([-12:4:12])
                xlim([8.5 length(varnames)+0.5])
                ylim([-2 2])
                xtickangle(45)
                xlabel('EEG features')
                ylabel('TD sleep-wake beta effect sizes')
                title('TD Biomarkers of conscious state','fontsize',18)
                
                print('-dpng','./Figures/SleepLMMMainEffectsBetas.png')
                print('-dsvg','./Figures/SleepLMMMainEffectsBetas.svg')
            case false
                %Do nothing
        end
        save('./TrainOnTD/TD_sleep_wake_stats.mat','TDsleepP', 'TDsleepF', 'TDsleepES')
        
    case 'MachineLearning'
        
        
        %% Regularized logistic regression
        load TD_sleep_wake_stats
        
        TDidx2 = find(TDidx);
        ASidx2 = find(ASidx);
        DSidx2 = find(DSidx);
                
        % Notes from JEFF
        % Do CV for regularization ONLY using training data, and when you do your ...
        % split for CV, make sure that both samples from the same subjects are in ...
        % either the test or the validation set. Don't spread between both.
               
        rng(58980423) % NOTE: The random number generator seed affects the training performance. Do not change!
        Nboot = 1e4; % bootstrapped resamples
        % results is group x [AUC, accuracy, PPV, RECALL, SPEC] x feature type
        % x [real value, CI lower bound, CI upper bound]
        results = nan(3,5,5,3);
        ccnt = 0; % counter for classifications
                
        % ground truth (is this person conscious? frame it this way so we are talking about a consciousness detector)
        ASgt = T.Conscious(ASidx) == categorical({'Yes'});
        DSgt = T.Conscious(DSidx) == categorical({'Yes'});
        TDgt = T.Conscious(TDidx) == categorical({'Yes'});
               
        feats = unique(vartype);
        LAMBDA = nan(1,length(feats)-1);
        
        % group x resample x feature x measure (AUC, accuracy, sensitivity,
        % specificity) 
        % bootAUC = nan(3,Nboot,length(feats)-1);
                 
        % Use n for each bootstrapped resample,
        TDboot = sum(TDidx);
        ASboot = sum(ASidx);
        DSboot = sum(DSidx);

        %%% TRAIN ON TD %%%
        fprintf('TD TRAINING\n')
        
        % ROC curve preallocation
        xTD = cell(5,1);
        xAS = cell(5,1);
        xDS = cell(5,1);
        yTD = cell(5,1);
        yAS = cell(5,1);
        yDS = cell(5,1);
        
        % Choose regularization parameter using 10-fold cross validation on
        % training data
        kfold = 10; % number of CV folds
        S = T.Subject(TDidx); % unique subjects

        us = unique(S);
        cvid = 1 + mod((1:length(us))',kfold);
        indices = randperm(length(us)); % randperm
        cvid = cvid(indices);

        mytest = nan(size(predTD,1),kfold);
        mytraining = nan(size(predTD,1),kfold);
        for icv = 1:kfold
            sidx = cvid == icv; % subjects belonging to test data for this CV
            mytest(:,icv) = ismember(S,us(sidx));
            mytraining(:,icv) = ~mytest(:,icv); % the ones not in the test data are training data
        end

        assert(all(sum(mytest,2) == ones(size(mytest,1),1)),'Each dataset must be used only once')
        assert(all(mytest + mytraining,[1 2]),'Some data are used both as test and training data')
                
        for ift = 2:length(feats)
     
            switch method
                case 'LMM'
                    feat  = feats{ift};
                    fprintf('\n%s FEATURES ONLY using %s method\n',upper(feat),method)
                    fts = find(ismember(vartype,feat));
                    ranks = find(abs(TDsleepES(1,ismember(vartype,feat))) >= betathresh); % use features with abs(beta) >= 0.5

                    tmp = T(:,fts(ranks));
                    predALL = T{:,fts(ranks)};
                    predTD = T{TDidx,fts(ranks)};
                    predAS = T{ASidx,fts(ranks)};
                    predDS = T{DSidx,fts(ranks)};
                    fprintf('     Using these features:')
                    for jft = 1:length(ranks)
                        fprintf(' %s ',tmp.Properties.VariableNames{jft})
                    end
                    fprintf('\n')
                case 'PCA'
                    feat  = feats{ift};
                    fprintf('\n%s FEATURES ONLY using %s method\n',upper(feat),method)
                    fts = find(ismember(vartype,feat));
                    feat  = replace(feats{ift},' ','_');
                    pcidx = contains(T.Properties.VariableNames,'PC') & contains(T.Properties.VariableNames,feat);
                    predTD = T{TDidx,pcidx};
                    predAS = T{ASidx,pcidx};
                    predDS = T{DSidx,pcidx};
            end
            
            %lambtest = logspace(-8,4,1e3); % test values for lambda
            [B0,FitInfo,CV] = mylassoglm(predTD,T.Conscious(TDidx),'binomial',...
               'Link','logit','CV',kfold,'Balanced',cat(3,mytraining,mytest));
 
           
            % Check that data from the same subjects are on same side of the CV partion
            IS = nan(1,kfold);
            for ifold = 1:kfold
                St = S(logical(CV.training(:,ifold))); % training subjects
                Sv = S(~logical(CV.training(:,ifold))); % test subjects
                IS(ifold) = length(intersect(St,Sv));
            end
            
            assert(max(IS)==0,'At least one fold had data from same subject on different sides of CV partition')
            
            % Examine the cross-validation plot to see the effect of the Lambda regularization parameter.
            lassoPlot(B0,FitInfo,'plottype','CV');
            legend('show') % Show legend
            
            [~,useme] = min(FitInfo.Deviance); % find optimal regularization
            lambda = FitInfo.Lambda(useme); % regularization parameter
            LAMBDA(ift-1) = lambda;
            
            % Do logistic regression just on TD data using the chosen lambda value
            [B,LassoFit] = lassoglm(predTD,T.Conscious(TDidx),'binomial','Link','logit','Lambda',lambda);
            
            coef = [LassoFit.Intercept; B];
            
            % apply model fit with optimal regularization to validation sets
            TDfit = glmval(coef,predTD,'logit','Constant','on');
            ASfit = glmval(coef,predAS,'logit','Constant','on');
            DSfit = glmval(coef,predDS,'logit','Constant','on');
            
            % Compare TD training performance with 10-fold cross validation
            k2fold = kfold; % use same number of CV folds as for the hyperparameter fitting
            
            cvid2 = 1 + mod((1:length(us))',k2fold);
            indices = randperm(length(us)); % randperm
            cvid2 = cvid2(indices);
            
            mytest2 = nan(sum(TDidx),k2fold);
            mytraining2 = nan(sum(TDidx),k2fold);
            for icv = 1:k2fold
                sidx = cvid2 == icv; % subjects belonging to test data for this CV
                mytest2(:,icv) = ismember(S,us(sidx));
                mytraining2(:,icv) = ~mytest2(:,icv); % the ones not in the test data are training data
            end
            
            assert(all(sum(mytest2,2) == ones(size(mytest2,1),1)),'Each dataset must be used only once')
            assert(all(mytest2 + mytraining2,[1 2]),'Some data are used both as test and training data')
            
            [B,LassoFit] = lassoglm(predTD,T.Conscious(TDidx),'binomial','Link','logit','Lambda',lambda);
            coef = [LassoFit.Intercept; B];
            TDfitcvx = [];
            TDgtcvx = [];
            theGT = logical(T.Conscious2(TDidx));
            for icv = 1:k2fold
                [Bx,LassoFitx] = lassoglm(predTD(logical(mytraining2(:,icv)),:),...
                    theGT(logical(mytraining2(:,icv))),'binomial','Link','logit','Lambda',lambda);
                tmpcoef = [LassoFitx.Intercept; Bx];
                TDfitcvx = [TDfitcvx; glmval(tmpcoef,predTD(logical(mytest2(:,icv)),:),'logit','Constant','on')];
                TDgtcvx = [TDgtcvx; TDgt(logical(mytest2(:,icv)))];
            end
            
            [~,~,~,TD_AUCcvx] = perfcurve(TDgtcvx,TDfitcvx,1,'XCrit','FPR',...
                'YCrit','TPR','Nboot',Nboot,'Alpha',0.05);
            
            %%% TD  %%%
            % Test to see if the curve has the correct number of points, 
            % N + 1, where N is the number of data sests
            
            [X_TD,Y_TD,T_TD,AUC_TD,opt_TD] = perfcurve(TDgt,TDfit,1,...
                'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Alpha',0.05);
            
            % Compare training AUC to k-fold CV AUC
            assert(size(TDgt,1)==size(TDgtcvx,1),'Not all subjects were tested with crossvalidation!?')
            assert(size(TDfit,1)==size(TDfitcvx,1),'Not all subjects were tested with crossvalidation!?')
            Nsamp = size(TDgt,1)/2;
            [pvl,u1,u2] = AUCMannWhitney(TD_AUCcvx(1),Nsamp,Nsamp,'right',AUC_TD(1),Nsamp,Nsamp)
            
            if pvl < 0.05
                fprintf('%i-fold cross-validation AUC (%1.2f) is significantly different (p = %1.3f) from the training AUC (%1.2f) for %s\n',...
                    k2fold,TD_AUCcvx(1),pvl,AUC_TD(1),feats{ift})
            else
                fprintf('No sig. difference between CV (AUC = %1.2f) and training (AUC = %1.2f) AUCs (p = %1.3f) for %s\n',...
                    TD_AUCcvx(1),AUC_TD(1),pvl,feats{ift})
            end
            
            pcv(ift-1) = pvl;
            u1cv(ift-1) = u1;
            u2cv(ift-1) = u2;
            AUC1cv(ift-1,:) = TD_AUCcvx; % cross-validated
            AUC2cv(ift-1,:) = AUC_TD; % training
            
              % Leave this commented out for now--it seems to do more harm
              % than good if we interpolate curves
%             if size(X_TD,1) ~= sum(TDidx) + 1
%                 warning('Bad number of points! See message below')
%                 fprintf('    Wrong number of test points, interpolating ROC curve\n')
%                 X_TD = interp1(linspace(0,1,length(X_TD)),X_TD,linspace(0,1,sum(TDidx)+1),'linear');
%                 tstvals = X_TD(:,1);
%                 [~,Y_TD,T_TD,AUC_TD,opt_TD] = perfcurve(TDgt,TDfit,1,...
%                     'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Xvals',tstvals,'Alpha',0.05);
%             end
            
            thresh = unique(T_TD(intersect(find(X_TD(:,1) == opt_TD(1)),find(Y_TD(:,1) == opt_TD(2)))));
            assert(length(thresh)==1,'More than one threshold found')
            TP = sum(TDgt & (TDfit >= thresh));
            FP = sum(~TDgt & (TDfit >= thresh));
            TN = sum(~TDgt & (TDfit < thresh));
            FN = sum(TDgt & (TDfit < thresh));
            sens_TD = TP/(TP+FN);
            spec_TD = TN/(TN+FP);
            ppv_TD = TP/(TP+FP);
            acc_TD = (TP+TN)/(TP+FP+TN+FN);

            % now do perfcurve again to get 95% CIs for accuracy,
            % sensitivity, specificity, etc.
            
            [FALL_TD,ACCU_TD] = perfcurve(TDgt,TDfit,1,'XCrit','FPR','YCrit','accu','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_TD(1));
            %if we had a wrong accuracy because of the weird operating point
            if opt_TD(1) == 0 
                [FALL_TD,ACCU_TD] = perfcurve(TDgt,TDfit,1,'XCrit','FPR','YCrit','accu','Nboot',Nboot,'Alpha',0.05);
                ACCU_TD = ACCU_TD(find(ACCU_TD(:,1) == acc_TD,1),:);
                [FALL_TD,PPV_TD] = perfcurve(TDgt,TDfit,1,'XCrit','FPR','YCrit','PPV','Nboot',Nboot,'Alpha',0.05);
                PPV_TD = PPV_TD(find(PPV_TD(:,1) == ppv_TD,1),:);
                [FALL_TD,RECALL_TD] = perfcurve(TDgt,TDfit,1,'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Alpha',0.05);
                RECALL_TD = RECALL_TD(find(RECALL_TD(:,1) == sens_TD,1),:);
                SPEC_TD = 1-Y_TD;
                SPEC_TD = SPEC_TD(find(SPEC_TD(:,1) == spec_TD,1),:);
                SPEC_TD = SPEC_TD([1 3 2]); % the CI bounds are swapped, flip them
            else
                try % this one will sometimes fail if AUC is very small
                    [FALL_TD,PPV_TD] = perfcurve(TDgt,TDfit,1,'XCrit','FPR','YCrit','PPV','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_TD(1));
                    assert(~any(isnan(PPV_TD)),'NaNs in PPV or its CI') % another specific bug I've encountered, only for DS
                catch
                    if ppv_TD == 0
                        PPV_TD = [0 0 0];
                    else
                        [FALL_TD,PPV_TD] = perfcurve(TDgt,TDfit,1,'XCrit','FPR','YCrit','PPV','Nboot',Nboot,'Alpha',0.05);
                        PPV_TD = PPV_TD(find(PPV_TD(:,1) == ppv_TD,1),:);
                    end
                end
                [FALL_TD,RECALL_TD] = perfcurve(TDgt,TDfit,1,'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_TD(1));
                SPEC_TD = 1-X_TD(find(X_TD(:,1)==opt_TD(1),1),:);
                SPEC_TD = SPEC_TD([1 3 2]); % the CI bounds are swapped, flip them
            end
            
            % sanity checks
            assert(single(ACCU_TD(1))==single(acc_TD),'Accuracies don''t match')
            assert(single(PPV_TD(1))==single(ppv_TD),'PPVs don''t match')
            assert(single(RECALL_TD(1))==single(sens_TD),'sensitivities (recall) don''t match')
            assert(single(SPEC_TD(1))==single(spec_TD),'specificities don''t match')
            
            % extract results
            results(1,1,ift-1,:) = AUC_TD.*100;
            results(1,2,ift-1,:) = ACCU_TD.*100;
            results(1,3,ift-1,:) = PPV_TD.*100;
            results(1,4,ift-1,:) = RECALL_TD.*100;
            results(1,5,ift-1,:) = SPEC_TD.*100;
                            
            fprintf('TD training set accuracy = %1.2f%%, AUC = %1.2f%%, precision = %1.2f%%, recall = %1.2f%%\n',...
                results(1,2,ift-1,1),results(1,1,ift-1,1),results(1,3,ift-1,1),results(1,4,ift-1,1))
            
            %%% AS %%%
            % Test to see if the curve has the correct number of points, 
            % N + 1, where N is the number of data sests
                        
            [X_AS,Y_AS,T_AS,AUC_AS,opt_AS] = perfcurve(ASgt,ASfit,1,...
                'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Alpha',0.05);
              
             % Leave this commented out for now--it seems to do more harm
              % than good if we interpolate curves
%             if size(X_AS,1) ~= sum(ASidx) + 1
%                 warning('Bad number of points! See message below')
%                 fprintf('    Wrong number of test points, interpolating ROC curve\n')
%                 X_AS = interp1(linspace(0,1,length(X_AS)),X_AS,linspace(0,1,sum(ASidx)+1),'linear');
%                 tstvals = X_AS(:,1);
%                 [~,Y_AS,T_AS,AUC_AS,opt_AS] = perfcurve(ASgt,ASfit,1,...
%                     'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Xvals',tstvals,'Alpha',0.05);
%             end
            
            thresh = unique(T_AS(intersect(find(X_AS(:,1) == opt_AS(1)),find(Y_AS(:,1) == opt_AS(2))),1));
            assert(length(thresh)==1,'More than one threshold found')
            TP = sum(ASgt & (ASfit >= thresh));
            FP = sum(~ASgt & (ASfit >= thresh));
            TN = sum(~ASgt & (ASfit < thresh));
            FN = sum(ASgt & (ASfit < thresh));
            sens_AS = TP/(TP+FN);
            spec_AS = TN/(TN+FP);
            ppv_AS = TP/(TP+FP);
            acc_AS = (TP+TN)/(TP+FP+TN+FN);

            % now do perfcurve again to get 95% CIs for accuracy,
            % sensitivity, specificity, etc.
                        
            [FALL_AS,ACCU_AS] = perfcurve(ASgt,ASfit,1,'XCrit','FPR','YCrit','accu','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_AS(1));
            %if we had a wrong accuracy because of the weird operating point
            if opt_AS(1) == 0 
                [FALL_AS,ACCU_AS] = perfcurve(ASgt,ASfit,1,'XCrit','FPR','YCrit','accu','Nboot',Nboot,'Alpha',0.05);
                ACCU_AS = ACCU_AS(find(ACCU_AS(:,1) == acc_AS,1),:);
                [FALL_AS,PPV_AS] = perfcurve(ASgt,ASfit,1,'XCrit','FPR','YCrit','PPV','Nboot',Nboot,'Alpha',0.05);
                PPV_AS = PPV_AS(find(PPV_AS(:,1) == ppv_AS,1),:);
                [FALL_AS,RECALL_AS] = perfcurve(ASgt,ASfit,1,'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Alpha',0.05);
                RECALL_AS = RECALL_AS(find(RECALL_AS(:,1) == sens_AS,1),:);
                SPEC_AS = 1-X_AS;
                SPEC_AS = SPEC_AS(find(SPEC_AS(:,1) == spec_AS,1),:);
                SPEC_AS = SPEC_AS([1 3 2]); % the CI bounds are swapped, flip them
            else
                try % this one will sometimes fail if AUC is very small
                    [FALL_AS,PPV_AS] = perfcurve(ASgt,ASfit,1,'XCrit','FPR','YCrit','PPV','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_AS(1));
                    assert(~any(isnan(PPV_AS)),'NaNs in PPV or its CI') % another specific bug I've encountered, only for DS
                catch
                    if ppv_AS == 0
                        PPV_AS = [0 0 0];
                    else
                        [FALL_AS,PPV_AS] = perfcurve(ASgt,ASfit,1,'XCrit','FPR','YCrit','PPV','Nboot',Nboot,'Alpha',0.05);
                        PPV_AS = PPV_AS(find(PPV_AS(:,1) == ppv_AS,1),:);
                    end
                end
                [FALL_AS,RECALL_AS] = perfcurve(ASgt,ASfit,1,'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_AS(1));
                SPEC_AS = 1-X_AS(find(X_AS(:,1)==opt_AS(1),1),:);
                SPEC_AS = SPEC_AS([1 3 2]); % the CI bounds are swapped, flip them
            end
            
            % sanity checks
            assert(single(ACCU_AS(1))==single(acc_AS),'Accuracies don''t match')
            assert(single(PPV_AS(1))==single(ppv_AS),'PPVs don''t match')
            assert(single(RECALL_AS(1))==single(sens_AS),'sensitivities (recall) don''t match')
            assert(single(SPEC_AS(1))==single(spec_AS),'specificities don''t match')
            
            % extract results
            results(2,1,ift-1,:) = AUC_AS.*100;
            results(2,2,ift-1,:) = ACCU_AS.*100;
            results(2,3,ift-1,:) = PPV_AS.*100;
            results(2,4,ift-1,:) = RECALL_AS.*100;
            results(2,5,ift-1,:) = SPEC_AS.*100;
                            
            fprintf('AS training set accuracy = %1.2f%%, AUC = %1.2f%%, precision = %1.2f%%, recall = %1.2f%%\n',...
                results(2,2,ift-1,1),results(2,1,ift-1,1),results(2,3,ift-1,1),results(2,4,ift-1,1))
            
            %%% Dup15q syndrome %%% 
            % Test to see if the curve has the correct number of points, 
            % N + 1, where N is the number of data sests
           
            [X_DS,Y_DS,T_DS,AUC_DS,opt_DS] = perfcurve(DSgt,DSfit,1,...
                'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Alpha',0.05);
            
              % Leave this commented out for now--it seems to do more harm
              % than good if we interpolate curves
%             if size(X_DS,1) ~= sum(DSidx) + 1
%                 warning('Bad number of points! See message below')
%                 fprintf('    Wrong number of test points, interpolating ROC curve\n')
%                 X_DS = interp1(linspace(0,1,length(X_DS)),X_DS,linspace(0,1,sum(DSidx)+1),'linear');
%                 tstvals = X_DS(:,1);
%                 [~,Y_DS,T_DS,AUC_DS,opt_DS] = perfcurve(DSgt,DSfit,1,...
%                     'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Xvals',tstvals,'Alpha',0.05);
%             end
                        
            thresh = unique(T_DS(intersect(find(X_DS(:,1) == opt_DS(1)),find(Y_DS(:,1) == opt_DS(2)))));
            assert(length(thresh)==1,'More than one threshold found')
            TP = sum(DSgt & (DSfit >= thresh));
            FP = sum(~DSgt & (DSfit >= thresh));
            TN = sum(~DSgt & (DSfit < thresh));
            FN = sum(DSgt & (DSfit < thresh));
            sens_DS = TP/(TP+FN);
            spec_DS = TN/(TN+FP);
            ppv_DS = TP/(TP+FP);
            acc_DS = (TP+TN)/(TP+FP+TN+FN);

            % now do perfcurve again to get 95% CIs for accuracy,
            % sensitivity, specificity, etc.
            
            [FALL_DS,ACCU_DS] = perfcurve(DSgt,DSfit,1,'XCrit','FPR','YCrit','accu','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_DS(1));
            %if we had a wrong accuracy because of the weird operating point
            if opt_DS(1) == 0 
                [FALL_DS,ACCU_DS] = perfcurve(DSgt,DSfit,1,'XCrit','FPR','YCrit','accu','Nboot',Nboot,'Alpha',0.05);
                ACCU_DS = ACCU_DS(find(ACCU_DS(:,1) == acc_DS,1),:);
                [FALL_DS,PPV_DS] = perfcurve(DSgt,DSfit,1,'XCrit','FPR','YCrit','PPV','Nboot',Nboot,'Alpha',0.05);
                PPV_DS = PPV_DS(find(PPV_DS(:,1) == ppv_DS,1),:);
                [FALL_DS,RECALL_DS] = perfcurve(DSgt,DSfit,1,'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Alpha',0.05);
                RECALL_DS = RECALL_DS(find(RECALL_DS(:,1) == sens_DS,1),:);
                SPEC_DS = 1-X_DS;
                SPEC_DS = SPEC_DS(find(SPEC_DS(:,1) == spec_DS,1),:);
                SPEC_DS = SPEC_DS([1 3 2]); % the CI bounds are swapped, flip them
            else
                try % this one will sometimes fail if AUC is very small
                    [FALL_DS,PPV_DS] = perfcurve(DSgt,DSfit,1,'XCrit','FPR','YCrit','PPV','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_DS(1));
                    assert(~any(isnan(PPV_DS)),'NaNs in PPV or its CI') % another specific bug I've encountered, only for DS
                catch
                    if ppv_DS == 0
                        PPV_DS = [0 0 0];
                    else
                        [FALL_DS,PPV_DS] = perfcurve(DSgt,DSfit,1,'XCrit','FPR','YCrit','PPV','Nboot',Nboot,'Alpha',0.05);
                        PPV_DS = PPV_DS(find(PPV_DS(:,1) == ppv_DS,1),:);
                    end
                end
                [FALL_DS,RECALL_DS] = perfcurve(DSgt,DSfit,1,'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_DS(1));
                SPEC_DS = 1-X_DS(find(X_DS(:,1)==opt_DS(1),1),:);
                SPEC_DS = SPEC_DS([1 3 2]); % the CI bounds are swapped, flip them
            end
            
            if any(isnan(PPV_DS)), keyboard, end
            
            % sanity checks
            assert(single(ACCU_DS(1))==single(acc_DS),'Accuracies don''t match')
            assert(single(PPV_DS(1))==single(ppv_DS),'PPVs don''t match')
            assert(single(RECALL_DS(1))==single(sens_DS),'sensitivities (recall) don''t match')
            assert(single(SPEC_DS(1))==single(spec_DS),'specificities don''t match')
            
            % extract results
            results(3,1,ift-1,:) = AUC_DS.*100;
            results(3,2,ift-1,:) = ACCU_DS.*100;
            results(3,3,ift-1,:) = PPV_DS.*100;
            results(3,4,ift-1,:) = RECALL_DS.*100;
            results(3,5,ift-1,:) = SPEC_DS.*100;
                            
            fprintf('DS training set accuracy = %1.2f%%, AUC = %1.2f%%, precision = %1.2f%%, recall = %1.2f%%\n',...
                results(3,2,ift-1,1),results(3,1,ift-1,1),results(3,3,ift-1,1),results(3,4,ift-1,1))
            fprintf('\n')
            

            xTD{ift-1} = X_TD(:,1);
            xAS{ift-1} = X_AS(:,1);
            xDS{ift-1} = X_DS(:,1);
            yTD{ift-1} = Y_TD(:,1);
            yAS{ift-1} = Y_AS(:,1);
            yDS{ift-1} = Y_DS(:,1);
                    
            myfigure2
            plot(X_TD(:,1),Y_TD(:,1),'LineWidth',2)
            plot(X_AS(:,1),Y_AS(:,1),'LineWidth',2)
            plot(X_DS(:,1),Y_DS(:,1),'LineWidth',2)
            plot([-0.02 1],[-0.02 1],'k:')
            legend({sprintf('TD, AUC = %2.1f%%',AUC_TD(1)*100),sprintf('AS, AUC = %2.1f%%',AUC_AS(1)*100),...
                sprintf('Dup15q, AUC = %2.1f%%',AUC_DS(1)*100),'Chance performance'},'location',...
                'northeastoutside','fontsize',12,'autoupdate','off')
            legend box off
            xlabel('False positive rate')
            ylabel('True positive rate')
            title(sprintf('ROC curve, %s with %s feat. selection',feats{ift},method),'fontsize',18)
            axis([0 1 0 1])
            xticks(0:0.2:1)
            yticks(0:0.2:1)
            makefighandsome
            axis square
            pause(0.01)
            print('-dpng',sprintf('./TrainOnTD/Figures/ROC_curve_%s_%s',feats{ift},method))

        end
        
        assert(all(single(results(:,:,:,1)) >= single(results(:,:,:,2)) & ...
            single(results(:,:,:,1)) <= single(results(:,:,:,3)),[1 2 3]),...
            'True values are not contained within 95% CIs!')
        %% Plot ROC curves separately for each group
        
        myfigure
        for i = 1:length(xTD)
            switch method             
                case 'PCA'
                    plot(xTD{i}',yTD{i}','LineWidth',2)
                case 'LMM'
                    plot(xTD{i}',yTD{i}',':','LineWidth',2)
            end
        end
        plot([-0.02 1],[-0.02 1],'k:')
        legend([feats(2:end) {'Chance performance'}],'location',...
            'northeastoutside','fontsize',14,'autoupdate','off')
        legend box off
        xlabel('False positive rate')
        ylabel('True positive rate')
        title(sprintf('TD ROC curves with %s feat. selection',method),'fontsize',18)
        axis([-0.02 1 -0.02 1])
        xticks(0:0.2:1)
        yticks(0:0.2:1)
        makefigpretty
        print('-dpng',sprintf('./TrainOnTD/Figures/TD_ROC_curve_%s.png',method))
        print('-dsvg',sprintf('./TrainOnTD/Figures/TD_ROC_curve_%s.svg',method))
        
        myfigure
        for i = 1:length(xAS)
            switch method             
                case 'PCA'
                    plot(xAS{i}',yAS{i}','LineWidth',2)
                case 'LMM'
                    plot(xAS{i}',yAS{i}',':','LineWidth',2)
            end
        end
        plot([-0.02 1],[-0.02 1],'k:')
        legend([feats(2:end) {'Chance performance'}],'location',...
            'northeastoutside','fontsize',14,'autoupdate','off')
        legend box off
        xlabel('False positive rate')
        ylabel('True positive rate')
        title(sprintf('AS ROC curves with %s feat. selection',method),'fontsize',18)
        axis([-0.02 1 -0.02 1])
        xticks(0:0.2:1)
        yticks(0:0.2:1)
        makefigpretty
        print('-dpng',sprintf('./TrainOnTD/Figures/AS_ROC_curve_%s.png',method))
        print('-dsvg',sprintf('./TrainOnTD/Figures/AS_ROC_curve_%s.svg',method))
        
        myfigure
        for i = 1:length(xDS)
            switch method             
                case 'PCA'
                    plot(xDS{i}',yDS{i}','LineWidth',2)
                case 'LMM'
                    plot(xDS{i}',yDS{i}',':','LineWidth',2)
            end
        end
        plot([-0.02 1],[-0.02 1],'k:')
        legend([feats(2:end) {'Chance performance'}],'location',...
            'northeastoutside','fontsize',14,'autoupdate','off')
        legend box off
        xlabel('False positive rate')
        ylabel('True positive rate')
        title(sprintf('Dup15q ROC curves with %s feat. selection',method),'fontsize',18)
        axis([-0.02 1 -0.02 1])
        xticks(0:0.2:1)
        yticks(0:0.2:1)
        makefigpretty
        print('-dpng',sprintf('./TrainOnTD/Figures/Dup15q_ROC_curve_%s.png',method))
        print('-dsvg',sprintf('./TrainOnTD/Figures/Dup15q_ROC_curve_%s.svg',method))
        
        
        %% Test for differences in AUC between measures
        % feature category names that we will use in the paper
        featnames = feats(2:end);
        
        % initialize vars
        tststr = cell(9,1);
        groupz = cell(9,1);
        e2sp = nan(9,1); % entropy to spectral comparison p-vales
        AUC1 = nan(9,1);
        AUC2 = nan(9,1);
        MWU1 = nan(9,1);
        MWU2 = nan(9,1);
        SampN = nan(9,1);
        winner = cell(9,1);
       
        % Consider the CV AUC separately for TD
        e2sp_cv = nan(3,1); % entropy to spectral comparison p-vales
        AUC1_cv = nan(3,1);
        AUC2_cv = nan(3,1);
        MWU1_cv = nan(3,1);
        MWU2_cv = nan(3,1);
        SampN_cv = nan(3,1);
        winner_cv = cell(3,1);
        groupz_cv = cell(3,1);
        teststr_cv = cell(3,1);

        cnt = 0;
        cntcv = 0;
        g = {'TD';'AS';'DS'};
        for ift = 1:size(results,3)
            for jft = 1:size(results,3)
                for igrp = 1:3
                    if contains(featnames{ift},'Entropy') && contains(featnames{jft},'Spectral') && ...
                            ((contains(featnames{ift},'sc') && contains(featnames{jft},'sc')) || ...
                            (contains(featnames{ift},'fc') && contains(featnames{jft},'fc')))
                        cnt = cnt + 1;
                        groupz{cnt} = g{igrp};
                        tststr{cnt} = sprintf('%s vs %s',featnames{ift},featnames{jft});
                        N = sum(strcmp(cellstr(T.Group),g{igrp}))/2;
                        
                        if strcmp(g{igrp},'TD')
                            % For TD only, also use the CV AUC
                            % Mann-Whitney
                            cntcv = cntcv + 1;
                            groupz_cv{cntcv} = g{igrp};
                            tststr_cv{cntcv} = sprintf('%s vs %s',featnames{ift},featnames{jft});
                            [pVal2,U1,U2] = AUCMannWhitney(AUC1cv(ift),N,N,'both',AUC1cv(jft),N,N);
                            e2sp_cv(cntcv) = pVal2;
                            AUC1_cv(cntcv) = AUC1cv(ift);
                            AUC2_cv(cntcv) = AUC1cv(jft);
                            MWU1_cv(cntcv) = U1;
                            MWU2_cv(cntcv) = U2;
                            SampN_cv(cnt) = N;
                            if U1 > U2
                                winner_cv{cntcv} = 'Entropy';
                            elseif U1 < U2
                                winner_cv{cntcv} = 'Spectral';
                            else
                                winner_cv{cntcv} = 'NA';
                            end
                        end
                        
                        % Also use Mann-Whitney
                        [pVal2,U1,U2] = AUCMannWhitney(results(igrp,1,ift,1)./100,N,N,'both',results(igrp,1,jft,1)./100,N,N);
                        e2sp(cnt) = pVal2;
                        AUC1(cnt) = results(igrp,1,ift)./100;
                        AUC2(cnt) = results(igrp,1,jft)./100;
                        MWU1(cnt) = U1;
                        MWU2(cnt) = U2;
                        SampN(cnt) = N;
                        if U1 > U2
                            winner{cnt} = 'Entropy';
                        elseif U1 < U2
                            winner{cnt} = 'Spectral';
                        else
                            winner{cnt} = 'NA';
                        end
                    end
                end
            end
        end
                
        % Table for main manuscript
        Te2s = table(groupz,SampN,tststr,winner,e2sp,AUC1,AUC2,MWU1,MWU2,'VariableNames',...
            {'Cohort','N','Test','Larger AUC','MWU p-vlue','Entropy AUC','Spectral AUC',...
            'Entropy U','Spectral U'})
        writetable(Te2s,sprintf('./TrainOnTD/EntrVsSpecPvalsTD_%s%i.csv',method,Nboot))
        
        Te2s_cv = table(groupz_cv,repmat(sum(TDgt),3,1),tststr_cv',winner_cv,...
            e2sp_cv,AUC1_cv,AUC2_cv,MWU1_cv,MWU2_cv,'VariableNames',...
            {'Cohort','N','Test','Larger AUC','MWU p-value','Entropy AUC','Spectral AUC',...
            'Entropy U','Spectral U'})
        writetable(Te2s_cv,sprintf('./TrainOnTD/EntrVsSpecPvals_10foldCV_%s%i.csv',method,Nboot))

        
        %% Create table of results
        
        allAUC = reshape(results(:,1,:,1),numel(results)/15,1)./100;
        allACC = reshape(results(:,2,:,1),numel(results)/15,1)./100;
        allPPV = reshape(results(:,3,:,1),numel(results)/15,1)./100; % same as precision
        allSENS = reshape(results(:,4,:,1),numel(results)/15,1)./100; % same as recall
        allSPEC = reshape(results(:,5,:,1),numel(results)/15,1)./100; % specificity
        
        features = [];
        for ift = 1:length(featnames)
            features = [features repmat(featnames(ift),1,3)];
        end
        
        groups = repmat({'NT';'AS';'DS'},size(results,3),1);
        
        cls = repmat({'Training';'Validation';'Validation'},size(results,3),1);
        
        f2 = repmat(unique(feats(~contains(feats,'IV')))',size(unique(feats(~contains(feats,'IV'))),1),3)';
        F = reshape(f2,1,size(f2,1)*size(f2,2))';
        
        LAMBDA2 = repmat(LAMBDA',size(LAMBDA,1),3)';
        L = reshape(LAMBDA2,1,size(LAMBDA2,1)*size(LAMBDA2,2))';
        
        AUC_CIs = squeeze([results(1,1,:,2:3); results(2,1,:,2:3); results(3,1,:,2:3)]);
        AUC_CI1 = reshape(AUC_CIs(:,:,1),15,1);
        AUC_CI2 = reshape(AUC_CIs(:,:,2),15,1);
        
        ACC_CIs = squeeze([results(1,2,:,2:3); results(2,2,:,2:3); results(3,2,:,2:3)]);
        ACC_CI1 = reshape(ACC_CIs(:,:,1),15,1);
        ACC_CI2 = reshape(ACC_CIs(:,:,2),15,1);
        
        PPV_CIs = squeeze([results(1,3,:,2:3); results(2,3,:,2:3); results(3,3,:,2:3)]);
        PPV_CI1 = reshape(PPV_CIs(:,:,1),15,1);
        PPV_CI2 = reshape(PPV_CIs(:,:,2),15,1);
        
        SENS_CIs = squeeze([results(1,4,:,2:3); results(2,4,:,2:3); results(3,4,:,2:3)]);
        SENS_CI1 = reshape(SENS_CIs(:,:,1),15,1);
        SENS_CI2 = reshape(SENS_CIs(:,:,2),15,1);
        
        SPEC_CIs = squeeze([results(1,5,:,2:3); results(2,5,:,2:3); results(3,5,:,2:3)]);
        SPEC_CI1 = reshape(SPEC_CIs(:,:,1),15,1);
        SPEC_CI2 = reshape(SPEC_CIs(:,:,2),15,1);
                
        % Create strings with confidence intervals
        
        AUCstrs = cell(size(allAUC,1),1);
        ACCstrs = cell(size(allACC,1),1);
        PPVstrs = cell(size(allPPV,1),1);
        SENSstrs = cell(size(allSENS,1),1);
        SPECstrs = cell(size(allSPEC,1),1);
        
        % Get p-values for AUCs based on Mann-Whitney U
        MWp = nan(size(allAUC));
        MWU = nan(size(allAUC));
        for irow = 1:size(allAUC,1)
            switch groups{irow}
                case 'NT'
                    [P,U] = AUCMannWhitney(allAUC(irow),sum(TDidx)/2,sum(TDidx)/2,'right');
                    MWp(irow) = P;
                    MWU(irow) = U;
                case 'AS'
                    [P,U] = AUCMannWhitney(allAUC(irow),sum(ASidx)/2,sum(ASidx)/2,'right');
                    MWp(irow) = P;
                    MWU(irow) = U;
                case 'DS'
                    [P,U] = AUCMannWhitney(allAUC(irow),sum(DSidx)/2,sum(DSidx)/2,'right');
                    MWp(irow) = P;
                    MWU(irow) = U;
            end
        end
        
        for irow = 1:size(allAUC,1)
            AUCstrs{irow,1} = sprintf('%2.1f%% (%2.1f%% - %2.1f%%)',allAUC(irow).*100,AUC_CI1(irow),AUC_CI2(irow));
            ACCstrs{irow,1} = sprintf('%2.1f%% (%2.1f%% - %2.1f%%)',allACC(irow).*100,ACC_CI1(irow),ACC_CI2(irow));
            PPVstrs{irow,1} = sprintf('%2.1f%% (%2.1f%% - %2.1f%%)',allPPV(irow).*100,PPV_CI1(irow),PPV_CI2(irow));
            SENSstrs{irow,1} = sprintf('%2.1f%% (%2.1f%% - %2.1f%%)',allSENS(irow).*100,SENS_CI1(irow),SENS_CI2(irow));
            SPECstrs{irow,1} = sprintf('%2.1f%% (%2.1f%% - %2.1f%%)',allSPEC(irow).*100,SPEC_CI1(irow),SPEC_CI2(irow));
        end
            
        
        % create big table (for supplement)
        Ta = table(F,cls,groups,L,MWp,MWU,AUCstrs,ACCstrs,PPVstrs,SENSstrs,SPECstrs,...
            'VariableNames',{'Features','Classifcation','Group','Lambda',...
            'MWU p-values','MW U (test-statistic)','AUC (95% CI)',...
            'Accuracy (95% CI)','Precision (95% CI)','Recall (95% CI)','Specificity (95% CI)'})
        writetable(Ta,sprintf('./TrainOnTD/ClassificationResultsTD%s_%i_resamples.csv',method,Nboot))
               
       % Just to be safe
        clear X_TD Y_TD T_TD AUC_TD opt_TD X_AS Y_AS T_AS AUC_AS opt_AS ...
            X_DS Y_DS T_DS AUC_DS opt_DS
        
        return
        
        %% PCA coefficients
        
        figure('color','w')
        h=subplot(1,1,1);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.5]);
        stem([PC_coefs{2,2}(:,1); PC_coefs{2,3}(:,1); PC_coefs{2,1}(:,1); ...
            PC_coefs{2,4}(:,1); PC_coefs{2,5}(:,1)],'filled','linewidth',2)
        xticks(1:length(varnames))
        xticklabels(varnames(frstvar:end))
        %yticks([-12:4:12])
        %xlim([0 length(varnames)+0.5])
        ylim([-1 1])
        xtickangle(45)
        %xlabel('EEG features')
        ylabel('PC_{1} loadings')
        
        box off
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', 20)
        set(xAX,'color','k')
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', 20)
        set(yAX,'color','k')
        set(gca, 'TickDir', 'out')
        set(gcf,'color','w')
        set(gca,'Layer','top')
        axis normal
        print('-dsvg','./Figures/PCA_coefficeints.svg')

        
        %% Table with multivariate features (all channels)
        
        varnamesx = {'Sleep','Conscious','Group','Subject','Sleep2',...
            'Conscious2','Group2','Age','mMSE','LZc',...
            'CTW','PermEn8','PermEn16','PermEn32','PermEn64',...
            'PermEn128','slow','delta1','delta2',...
            'theta','alpha','beta','slowR','delta1R','delta2R',...
            'thetaR','alphaR','betaR'};
        
%         % Multichannel criticality
%         
%         TDWakeCrit = criticality(TDWakeChaos,alpha(alphaidx));
%         TDSleepCrit = criticality(TDSleepChaos,alpha(alphaidx));
%         ASWakeCrit = criticality(ASWakeChaos,alpha(alphaidx));
%         ASSleepCrit = criticality(ASSleepChaos,alpha(alphaidx));
%         DSWakeCrit = criticality(DSWakeChaos,alpha(alphaidx));
%         DSSleepCrit = criticality(DSSleepChaos,alpha(alphaidx));
        
        % Multichannel permutation entropy

                
        TDPermEnwkTau8MC  = squeeze(TDWakePermEn(:,1,:));
        TDPermEnspTau8MC  = squeeze(TDSleepPermEn(:,1,:));
        DSPermEnwkTau8MC  = squeeze(DSWakePermEn(:,1,:));
        DSPermEnspTau8MC  = squeeze(DSSleepPermEn(:,1,:));
        ASPermEnwkTau8MC  = squeeze(ASWakePermEn(:,1,:));
        ASPermEnspTau8MC  = squeeze(ASSleepPermEn(:,1,:));
         
        TDPermEnwkTau16MC  = squeeze(TDWakePermEn(:,2,:));
        TDPermEnspTau16MC  = squeeze(TDSleepPermEn(:,2,:));
        DSPermEnwkTau16MC  = squeeze(DSWakePermEn(:,2,:));
        DSPermEnspTau16MC  = squeeze(DSSleepPermEn(:,2,:));
        ASPermEnwkTau16MC  = squeeze(ASWakePermEn(:,2,:));
        ASPermEnspTau16MC  = squeeze(ASSleepPermEn(:,2,:));
        
        TDPermEnwkTau32MC  = squeeze(TDWakePermEn(:,3,:));
        TDPermEnspTau32MC  = squeeze(TDSleepPermEn(:,3,:));
        DSPermEnwkTau32MC  = squeeze(DSWakePermEn(:,3,:));
        DSPermEnspTau32MC  = squeeze(DSSleepPermEn(:,3,:));
        ASPermEnwkTau32MC  = squeeze(ASWakePermEn(:,3,:));
        ASPermEnspTau32MC  = squeeze(ASSleepPermEn(:,3,:));
                
        TDPermEnwkTau64MC  = squeeze(TDWakePermEn(:,4,:));
        TDPermEnspTau64MC  = squeeze(TDSleepPermEn(:,4,:));
        DSPermEnwkTau64MC  = squeeze(DSWakePermEn(:,4,:));
        DSPermEnspTau64MC  = squeeze(DSSleepPermEn(:,4,:));
        ASPermEnwkTau64MC  = squeeze(ASWakePermEn(:,4,:));
        ASPermEnspTau64MC  = squeeze(ASSleepPermEn(:,4,:));
        
        TDPermEnwkTau128MC  = squeeze(TDWakePermEn(:,5,:));
        TDPermEnspTau128MC  = squeeze(TDSleepPermEn(:,5,:));
        DSPermEnwkTau128MC  = squeeze(DSWakePermEn(:,5,:));
        DSPermEnspTau128MC  = squeeze(DSSleepPermEn(:,5,:));
        ASPermEnwkTau128MC  = squeeze(ASWakePermEn(:,5,:));
        ASPermEnspTau128MC  = squeeze(ASSleepPermEn(:,5,:));
        
        % Remove criticality and chaoticity from analysis
%         [TDWakeCrit'; ASWakeCrit'; DSWakeCrit'; TDSleepCrit'; ASSleepCrit'; DSSleepCrit'],...
%             [TDWakeChaos'; ASWakeChaos'; DSWakeChaos'; TDSleepChaos'; ASSleepChaos'; DSSleepChaos'],...
%TD15qcn'; AS15qcn'; DS15qcn'; TD15qcn'; AS15qcn'; DS15qcn'],...
        
        Tx = table(issleep,isconscious,group,subjects,sleepBin,groupBin,consciousBin,...
            [TDages'; ASages'; DSages'; TDages'; ASages'; DSages'],...
            [TDWakeMSE'; ASWakeMSE'; DSWakeMSE'; TDSleepMSE'; ASSleepMSE'; DSSleepMSE'],...
            [TDWakeLZc'; ASWakeLZc'; DSWakeLZc'; TDSleepLZc'; ASSleepLZc'; DSSleepLZc'],...
            [TDWakeCTW'; ASWakeCTW'; DSWakeCTW'; TDSleepCTW'; ASSleepCTW'; DSSleepCTW'],...
            [TDPermEnwkTau8MC'; ASPermEnwkTau8MC'; DSPermEnwkTau8MC'; TDPermEnspTau8MC'; ASPermEnspTau8MC'; DSPermEnspTau8MC'],...
            [TDPermEnwkTau16MC'; ASPermEnwkTau16MC'; DSPermEnwkTau16MC'; TDPermEnspTau16MC'; ASPermEnspTau16MC'; DSPermEnspTau16MC'],...
            [TDPermEnwkTau32MC'; ASPermEnwkTau32MC'; DSPermEnwkTau32MC'; TDPermEnspTau32MC'; ASPermEnspTau32MC'; DSPermEnspTau32MC'],...
            [TDPermEnwkTau64MC'; ASPermEnwkTau64MC'; DSPermEnwkTau64MC'; TDPermEnspTau64MC'; ASPermEnspTau64MC'; DSPermEnspTau64MC'],...
            [TDPermEnwkTau128MC'; ASPermEnwkTau128MC'; DSPermEnwkTau128MC'; TDPermEnspTau128MC'; ASPermEnspTau128MC'; DSPermEnspTau128MC'],...
            [slow_wk_TD'; slow_wk_AS'; slow_wk_DS'; slow_wk_TD'; slow_sp_AS'; slow_sp_DS'],...
            [delta1_wk_TD'; delta1_wk_AS'; delta1_wk_DS'; delta1_sp_TD'; delta1_sp_AS'; delta1_sp_DS'],...
            [delta2_wk_TD'; delta2_wk_AS'; delta2_wk_DS'; delta2_sp_TD'; delta2_sp_AS'; delta2_sp_DS'],...
            [theta_wk_TD'; theta_wk_AS'; theta_wk_DS'; theta_sp_TD'; theta_sp_AS'; theta_sp_DS'],...
            [alpha_wk_TD'; alpha_wk_AS'; alpha_wk_DS'; alpha_sp_TD'; alpha_sp_AS'; alpha_sp_DS'],...
            [beta_wk_TD'; beta_wk_AS'; beta_wk_DS'; beta_sp_TD'; beta_sp_AS'; beta_sp_DS'],...
            [slow_rel_wk_TD'; slow_rel_wk_AS'; slow_rel_wk_DS'; slow_rel_sp_TD';slow_rel_sp_AS'; slow_rel_sp_DS'],...
            [delta1_rel_wk_TD'; delta1_rel_wk_AS'; delta1_rel_wk_DS'; delta1_rel_sp_TD';delta1_rel_sp_AS'; delta1_rel_sp_DS'],...
            [delta2_rel_wk_TD'; delta2_rel_wk_AS'; delta2_rel_wk_DS'; delta2_rel_sp_TD';delta2_rel_sp_AS'; delta2_rel_sp_DS'],...
            [theta_rel_wk_TD'; theta_rel_wk_AS'; theta_rel_wk_DS'; theta_rel_sp_TD';theta_rel_sp_AS'; theta_rel_sp_DS'],...
            [alpha_rel_wk_TD'; alpha_rel_wk_AS'; alpha_rel_wk_DS'; alpha_rel_sp_TD';alpha_rel_sp_AS'; alpha_rel_sp_DS'],...
            [beta_rel_wk_TD'; beta_rel_wk_AS'; beta_rel_wk_DS'; beta_rel_sp_TD';beta_rel_sp_AS'; beta_rel_sp_DS'],'VariableNames',varnamesx);
       
            % remove rows with NaNs
            Tx(any(isnan(Tx{:,max(find(strcmp(vartype,'IV')))+1:end}')),:) = []
        
            writetable(Tx,'./TrainOnTD/MultichannelEEGFeaturesTD.csv')
            
             
%         %% Look at AUC topography for vars
%         
%         varnamesy = Tx.Properties.VariableNames(max(find(strcmp(vartype,'IV')))+1:end);
%         
%         AUC = nan(3,length(varnamesy),19);
%                
%         for ivar = max(find(strcmp(vartype,'IV')))+1:length(varnamesx)
%             % Choose regularization parameter using 10-fold cross validation on all
%             % training data (TD)
%             
%             eval(sprintf('xvar = Tx.%s;',varnamesx{ivar}))
%             fprintf('Univariate prediction using %s\n',varnamesx{ivar})
%             for ich = 1:size(xvar,2) % for each channel
%             
%                 % Choose regularization parameter using 10-fold cross validation on all
%                 % training data (TD)
%                 [B0,FitInfo] = mylassoglm(xvar(TDidx,ich),Tx.Sleep(TDidx),'binomial','Link','logit','CV',kfold,'Balanced',cat(3,mytraining,mytest));
% 
%                 % Examine the cross-validation plot to see the effect of the Lambda regularization parameter.
%                 % lassoPlot(B0,FitInfo,'plottype','CV');
%                 % legend('show') % Show legend
% 
%                 [~,useme] = min(FitInfo.Deviance); % find optimal regularization
%                 lambda = FitInfo.Lambda(useme) % regularization parameter
% 
%                 % Do logistic regression just on TD data using the chosen lambda value
%                 [B,LassoFit] = lassoglm(xvar(TDidx,ich),Tx.Sleep(TDidx),'binomial','Link','logit','Lambda',lambda);
% 
%                 coef = [LassoFit.Intercept; B];
% 
%                 % apply model fit with optimal regularization to validation sets
%                 TDfit = glmval(coef,xvar(TDidx,ich),'logit','Constant','on');
%                 ASfit = glmval(coef,xvar(ASidx,ich),'logit','Constant','on');
%                 DSfit = glmval(coef,xvar(DSidx,ich),'logit','Constant','on');
% 
% 
%                 fprintf('Channel #%i: %s\n',ich,lay.label{ich})
%                 [X_TD,Y_TD,T_TD,AUC_TD,opt_TD] = perfcurve(TDgt,TDfit,1);
%                 thresh = T_TD(intersect(find(X_TD == opt_TD(1)),find(Y_TD == opt_TD(2))));
%                 acc_TD = sum(TDgt == (TDfit >= thresh))/length(TDfit)*100;
%                 fprintf('     TD training set accuracy = %1.2f%%, AUC = %1.2f%%\n',acc_TD,AUC_TD*100)
% 
%                 [X_AS,Y_AS,T_AS,AUC_AS,opt_AS] = perfcurve(ASgt,ASfit,1);
%                 thresh = T_AS(intersect(find(X_AS == opt_AS(1)),find(Y_AS == opt_AS(2))));
%                 acc_AS = sum(ASgt == (ASfit >= thresh))/length(ASfit)*100;
%                 fprintf('    AS validation set accuracy = %1.2f%%, AUC = %1.2f%%\n',acc_AS,AUC_AS*100)
% 
%                 [X_DS,Y_DS,T_DS,AUC_DS,opt_DS] = perfcurve(DSgt,DSfit,1);
%                 thresh = T_DS(intersect(find(X_DS == opt_DS(1)),find(Y_DS == opt_DS(2))));
%                 acc_DS = sum(DSgt == (DSfit >= thresh))/length(DSfit)*100;
%                 fprintf('    Dup15q validation set accuracy = %1.2f%%, AUC = %1.2f%%\n',acc_DS,AUC_DS*100)
%                 fprintf('\n')
%                 
%                 AUC(1,ivar,ich) = AUC_TD;
%                 AUC(2,ivar,ich) = AUC_AS;
%                 AUC(3,ivar,ich) = AUC_DS;
%             end
%             myfigure
%             plot_topo_AS(squeeze(AUC(1,ivar,:)))
%             caxis([0 1])
%             colormap jet
%             mycolorbar
%             title(sprintf('TD AUC %s',varnamesx{ivar}),'fontsize',18)
%             print('-dpng',sprintf('./TrainOnTD/Figures/TD_AUC_%s.png',varnamesx{ivar}))
%             
%             myfigure
%             plot_topo_AS(squeeze(AUC(2,ivar,:)))
%             caxis([0 1])
%             colormap jet
%             mycolorbar
%             title(sprintf('AS AUC %s',varnamesx{ivar}),'fontsize',18)
%             print('-dpng',sprintf('./TrainOnTD/Figures/AS_AUC_%s.png',varnamesx{ivar}))   
%             
%             myfigure
%             plot_topo_AS(squeeze(AUC(3,ivar,:)))
%             caxis([0 1])
%             colormap jet
%             mycolorbar
%             title(sprintf('Dup15q AUC %s',varnamesx{ivar}),'fontsize',18)
%             print('-dpng',sprintf('./TrainOnTD/Figures/Dup15q_AUC_%s.png',varnamesx{ivar}))
%             
%             pause(0.01)
%             
%         end
%         
%         % Just to be safe
%         clear X_TD Y_TD T_TD AUC_TD opt_TD X_AS Y_AS T_AS AUC_AS opt_AS ...
%             X_DS Y_DS T_DS AUC_DS opt_DS
%         
%         
%         %% Look at AUC topography for PCs
%         
%         AUC_PC = nan(3,3,19);
%         f = {'scEntropy','scSpectralA','scSpectralR'};
%         % Readload table
%         Tx = readtable('MultichannelEEGFeatures.csv');
%         for ift = 1:3
%             XSLEEP = nan(19,size(Tx,1)/2);
%             XWAKE  = nan(19,size(Tx,1)/2);
%             
%             fprintf('PC1 for %s\n',f{ift})
%             for ich = 1:19
%                 vrs = varnames(contains(vartype,f{ift}));
%                 if strcmp(f{ift},'scSpectralA')
%                     % Make sure relative power is not included (exclude
%                     % "R")
%                     NDX = endsWith(Tx.Properties.VariableNames,sprintf('_%i',ich)) & ...
%                         contains(Tx.Properties.VariableNames,vrs) & ~contains(Tx.Properties.VariableNames,'R');
%                 else
%                     NDX = endsWith(Tx.Properties.VariableNames,sprintf('_%i',ich)) & contains(Tx.Properties.VariableNames,vrs);
%                 end
%                 Tx.Properties.VariableNames(NDX);
%                 spvars = Tx{logical(Tx.Sleep2),NDX};
%                 wkvars = Tx{~logical(Tx.Sleep2),NDX};
%                 XSLEEP(ich,:)  = spvars*PC_coefs{2,contains(PC_coefs(1,:),f{ift})}(:,1);
%                 XWAKE(ich,:)  = wkvars*PC_coefs{2,contains(PC_coefs(1,:),f{ift})}(:,1);
%                 
%                 % Check this assumption below, necessary for index 
%                 assert(all(TDidx(1:end/2) == TDidx(end/2+1:end)),...
%                     'TD subject indices not symmatrical across wake/sleep')
%                 assert(all(ASidx(1:end/2) == ASidx(end/2+1:end)),...
%                     'AS subject indices not symmatrical across wake/sleep')
%                 assert(all(DSidx(1:end/2) == DSidx(end/2+1:end)),...
%                     'Dup15q subject indices not symmatrical across wake/sleep')
%                 
%                 % Choose regularization parameter using 10-fold cross validation on all
%                 % training data (TD)
%                 [B0,FitInfo] = mylassoglm([XWAKE(ich,TDidx(1:end/2))  XSLEEP(ich,TDidx(1:end/2))]',Tx.Sleep2(TDidx),'binomial','Link','logit','CV',kfold,'Balanced',cat(3,mytraining,mytest));
% 
%                 [~,useme] = min(FitInfo.Deviance); % find optimal regularization
%                 lambda = FitInfo.Lambda(useme); % regularization parameter
% 
%                 % Do logistic regression just on TD data using the chosen lambda value
%                 [B,LassoFit] = lassoglm([XWAKE(ich,TDidx(1:end/2))  XSLEEP(ich,TDidx(1:end/2))]',Tx.Sleep2(TDidx),'binomial','Link','logit','Lambda',lambda);
% 
%                 coef = [LassoFit.Intercept; B];
% 
%                 % apply model fit with optimal regularization to validation sets
%                 TDfit = glmval(coef,[XWAKE(ich,TDidx(1:end/2))  XSLEEP(ich,TDidx(1:end/2))]','logit','Constant','on');
%                 ASfit = glmval(coef,[XWAKE(ich,ASidx(1:end/2))  XSLEEP(ich,ASidx(1:end/2))]','logit','Constant','on');
%                 DSfit = glmval(coef,[XWAKE(ich,DSidx(1:end/2))  XSLEEP(ich,DSidx(1:end/2))]','logit','Constant','on');
% 
%                 fprintf('Channel #%i: %s\n',ich,lay.label{ich})
%                 [X_TD,Y_TD,T_TD,AUC_TD,opt_TD] = perfcurve(TDgt,TDfit,1);
%                 thresh = T_TD(intersect(find(X_TD == opt_TD(1)),find(Y_TD == opt_TD(2))));
%                 acc_TD = sum(TDgt == (TDfit >= thresh))/length(TDfit)*100;
%                 fprintf('     TD training set accuracy = %1.2f%%, AUC = %1.2f%%\n',acc_TD,AUC_TD*100)
% 
%                 [X_AS,Y_AS,T_AS,AUC_AS,opt_AS] = perfcurve(ASgt,ASfit,1);
%                 thresh = T_AS(intersect(find(X_AS == opt_AS(1)),find(Y_AS == opt_AS(2))));
%                 acc_AS = sum(ASgt == (ASfit >= thresh))/length(ASfit)*100;
%                 fprintf('    AS validation set accuracy = %1.2f%%, AUC = %1.2f%%\n',acc_AS,AUC_AS*100)
% 
%                 [X_DS,Y_DS,T_DS,AUC_DS,opt_DS] = perfcurve(DSgt,DSfit,1);
%                 thresh = T_DS(intersect(find(X_DS == opt_DS(1)),find(Y_DS == opt_DS(2))));
%                 acc_DS = sum(DSgt == (DSfit >= thresh))/length(DSfit)*100;
%                 fprintf('    Dup15q validation set accuracy = %1.2f%%, AUC = %1.2f%%\n',acc_DS,AUC_DS*100)
%                 fprintf('\n')
%                 
%                 AUC_PC(1,ift,ich) = AUC_TD;
%                 AUC_PC(2,ift,ich) = AUC_AS;
%                 AUC_PC(3,ift,ich) = AUC_DS;
%             end
%             % Topoplots
%             
%             % TD
%             myfigure
%             plot_topo_AS(squeeze(AUC_PC(1,ift,:)))
%             caxis([0 1])
%             colormap jet
%             mycolorbar
%             title(sprintf('TD_AUC_PC1_%s',f{ift}),'fontsize',18)
%             print('-dpng',sprintf('./TrainOnTD/Figures/TD_AUC_PC1_%s.png',f{ift}))
%             
%             % AS
%             myfigure
%             plot_topo_AS(squeeze(AUC_PC(2,ift,:)))
%             caxis([0 1])
%             colormap jet
%             mycolorbar
%             title(sprintf('AS_AUC_PC1_%s',f{ift}),'fontsize',18)
%             print('-dpng',sprintf('./TrainOnTD/Figures/AS_AUC_PC1_%s.png',f{ift}))
%             
%             % Dup15q
%             myfigure
%             plot_topo_AS(squeeze(AUC_PC(3,ift,:)))
%             caxis([0 1])
%             colormap jet
%             mycolorbar
%             title(sprintf('Dup15q_AUC_PC1_%s',f{ift}),'fontsize',18)
%             print('-dpng',sprintf('./TrainOnTD/Figures/Dup15q_AUC_PC1_%s.png',f{ift}))
%             
%         end
%                 
%         
%         %% Look at AUC for connectivity principal components
%         
%         idx = strcmp(PC_coefs(1,:),'fcEntropy');
%         % Average short-range and long-range coefficients
%         wSMI_PC_COEF = mean([PC_coefs{2,idx}(1:size(PC_coefs{2,idx},1)/2,1) ...
%             PC_coefs{2,idx}(size(PC_coefs{2,idx},1)/2+1:end,1)],2);
%         r = corr([PC_coefs{2,idx}(1:size(PC_coefs{2,idx},1)/2,1) PC_coefs{2,idx}(size(PC_coefs{2,idx},1)/2+1:end,1)]);
%         fprintf('Correlation between short-range and long-range wSMI coefficients is %1.2f\n',r(2,1))
%         
%         idx = strcmp(PC_coefs(1,:),'fcSpectral');
%         % Average short-range and long-range coefficients
%         dwPLI_PC_COEF = mean([PC_coefs{2,idx}(1:size(PC_coefs{2,idx},1)/2,1) ...
%             PC_coefs{2,idx}(size(PC_coefs{2,idx},1)/2+1:end,1)],2);
%         r = corr([PC_coefs{2,idx}(1:size(PC_coefs{2,idx},1)/2,1) PC_coefs{2,idx}(size(PC_coefs{2,idx},1)/2+1:end,1)]);
%         fprintf('Correlation between short-range and long-range dwPLI coefficients is %1.2f\n',r(2,1))
%         
%         % channel x channel x subject
%         
%         % Detect and remove nan subjects
%              
%         assert(size(TDSleepdwPLI,4)==size(TDWakewSMI,4),'Unequal number of subjects')
%         assert(size(ASSleepdwPLI,4)==size(ASWakewSMI,4),'Unequal number of subjects')
%         assert(size(DSSleepdwPLI,4)==size(DSWakewSMI,4),'Unequal number of subjects')
% 
%         rmv = find(isnan(squeeze(TDSleepwSMI(2,1,3,:))) | isnan(squeeze(TDWakewSMI(2,1,3,:)))...
%             | isnan(squeeze(TDSleepdwPLI(2,1,3,:))) | isnan(squeeze(TDWakedwPLI(2,1,3,:))));
%         fprintf('Removing %i matrix layers with all NaNs (TD)\n',length(rmv))
%         TDSleepwSMI(:,:,:,rmv) = [];
%         TDWakewSMI(:,:,:,rmv) = [];
%         TDSleepdwPLI(:,:,:,rmv) = [];
%         TDWakedwPLI(:,:,:,rmv) = [];
%         
%         rmv = find(isnan(squeeze(ASSleepwSMI(2,1,3,:))) | isnan(squeeze(ASWakewSMI(2,1,3,:)))...
%             | isnan(squeeze(ASSleepdwPLI(2,1,3,:))) | isnan(squeeze(ASWakedwPLI(2,1,3,:))));
%         fprintf('Removing %i matrix layers with all NaNs (AS)\n',length(rmv))
%         ASSleepwSMI(:,:,:,rmv) = [];
%         ASWakewSMI(:,:,:,rmv) = [];
%         ASSleepdwPLI(:,:,:,rmv) = [];
%         ASWakedwPLI(:,:,:,rmv) = [];
%         
%         rmv = find(isnan(squeeze(DSSleepwSMI(2,1,3,:))) | isnan(squeeze(DSWakewSMI(2,1,3,:)))...
%             | isnan(squeeze(DSSleepdwPLI(2,1,3,:))) | isnan(squeeze(DSWakedwPLI(2,1,3,:))));
%         fprintf('Removing %i matrix layers with all NaNs (Dup15q)\n',length(rmv))
%         DSSleepwSMI(:,:,:,rmv) = [];
%         DSWakewSMI(:,:,:,rmv) = [];
%         DSSleepdwPLI(:,:,:,rmv) = [];
%         DSWakedwPLI(:,:,:,rmv) = [];
%         
%         TD_SLEEP_dwPLI = nan(size(TDSleepdwPLI,[1 2 4]));
%         TD_WAKE_dwPLI = nan(size(TDSleepdwPLI,[1 2 4]));
%         TD_SLEEP_wSMI = nan(size(TDSleepwSMI,[1 2 4]));
%         TD_WAKE_wSMI = nan(size(TDSleepwSMI,[1 2 4]));
%         TD_dwPLI_AUC = nan(19,19);
%         TD_wSMI_AUC = nan(19,19);
%         
%         
%         AS_SLEEP_dwPLI = nan(size(ASSleepdwPLI,[1 2 4]));
%         AS_WAKE_dwPLI = nan(size(ASSleepdwPLI,[1 2 4]));
%         AS_SLEEP_wSMI = nan(size(ASSleepwSMI,[1 2 4]));
%         AS_WAKE_wSMI = nan(size(ASSleepwSMI,[1 2 4]));
%         AS_dwPLI_AUC = nan(19,19);
%         AS_wSMI_AUC = nan(19,19);
%         
%         DS_SLEEP_dwPLI = nan(size(DSSleepdwPLI,[1 2 4]));
%         DS_WAKE_dwPLI = nan(size(DSSleepdwPLI,[1 2 4]));
%         DS_SLEEP_wSMI = nan(size(DSSleepwSMI,[1 2 4]));
%         DS_WAKE_wSMI = nan(size(DSSleepwSMI,[1 2 4]));
%         DS_dwPLI_AUC = nan(19,19);
%         DS_wSMI_AUC = nan(19,19);
%         
%         % compute PC for each channel combo
%         chcnt = 0; % counter for reporting progress
%         for ich = 1:19          
%             for jch = 1:19               
%                 % if they're not neighbors; only do one half triangle
%                 if ich > jch && ~neighbors(ich,jch)
%                     chcnt = chcnt + 1;
%                     fprintf('%i%% complete\n',round(chcnt/((19^2-19)/2)*100))
%                     % TD
%                     TD_SLEEP_wSMI(ich,jch,:)  = wSMI_PC_COEF'*squeeze(TDSleepwSMI(ich,jch,:,:));
%                     TD_WAKE_wSMI(ich,jch,:)   = wSMI_PC_COEF'*squeeze(TDWakewSMI(ich,jch,:,:));
%                     TD_SLEEP_dwPLI(ich,jch,:) = dwPLI_PC_COEF'*squeeze(TDSleepdwPLI(ich,jch,:,:));
%                     TD_WAKE_dwPLI(ich,jch,:)  = dwPLI_PC_COEF'*squeeze(TDWakedwPLI(ich,jch,:,:));
%                     % Define the IV (sleep)
%                     s = [repmat(categorical({'Wake'}),size(TD_WAKE_wSMI,3),1);...
%                         repmat(categorical({'Sleep'}),size(TD_SLEEP_wSMI,3),1)];
%                     TDgt2 = s == categorical({'Sleep'});
%                     
%                     %%%% Find AUC for wSMI %%%%
%                     tmp_wSMI = [squeeze(TD_WAKE_wSMI(ich,jch,:)); squeeze(TD_SLEEP_wSMI(ich,jch,:))];
%                     [B0,FitInfo] = mylassoglm(tmp_wSMI,s,'binomial',...
%                         'Link','logit','CV',kfold,'Balanced',cat(3,mytraining,mytest));
%                                         
%                     [~,useme] = min(FitInfo.Deviance); % find optimal regularization
%                     lambda = FitInfo.Lambda(useme);% regularization parameter
%                     
%                     % Do logistic regression just on TD data using the chosen lambda value
%                     [B,LassoFit] = lassoglm(tmp_wSMI,s,'binomial','Link','logit','Lambda',lambda);
%                     
%                     coef = [LassoFit.Intercept; B];
%                     
%                     % apply model fit with optimal regularization to validation sets
%                     TDfit = glmval(coef,tmp_wSMI,'logit','Constant','on');
%                     
%                     [~,~,~,AUC_TD] = perfcurve(TDgt2,TDfit,1);
%                     TD_wSMI_AUC(ich,jch) = AUC_TD;
%                     
%                     %%%% Find AUC for dwPLI %%%%
%                     tmp_dwPLI = [squeeze(TD_WAKE_dwPLI(ich,jch,:)); squeeze(TD_SLEEP_dwPLI(ich,jch,:))];
%                     [B0,FitInfo] = mylassoglm(tmp_dwPLI,s,'binomial',...
%                         'Link','logit','CV',kfold,'Balanced',cat(3,mytraining,mytest));
%                                         
%                     [~,useme] = min(FitInfo.Deviance); % find optimal regularization
%                     lambda = FitInfo.Lambda(useme);% regularization parameter
%                     
%                     % Do logistic regression just on TD data using the chosen lambda value
%                     [B,LassoFit] = lassoglm(tmp_dwPLI,s,'binomial','Link','logit','Lambda',lambda);
%                     
%                     coef = [LassoFit.Intercept; B];
%                     
%                     % apply model fit with optimal regularization to validation sets
%                     TDfit = glmval(coef,tmp_dwPLI,'logit','Constant','on');
%                     
%                     [~,~,~,AUC_TD] = perfcurve(TDgt2,TDfit,1);
%                     TD_dwPLI_AUC(ich,jch) = AUC_TD;
%                                        
%                     % AS
%                     AS_SLEEP_wSMI(ich,jch,:)  = wSMI_PC_COEF'*squeeze(ASSleepwSMI(ich,jch,:,:));
%                     AS_WAKE_wSMI(ich,jch,:)   = wSMI_PC_COEF'*squeeze(ASWakewSMI(ich,jch,:,:));
%                     AS_SLEEP_dwPLI(ich,jch,:) = dwPLI_PC_COEF'*squeeze(ASSleepdwPLI(ich,jch,:,:));
%                     AS_WAKE_dwPLI(ich,jch,:)  = dwPLI_PC_COEF'*squeeze(ASWakedwPLI(ich,jch,:,:));
%                     % Define the IV (sleep)
%                     s = [repmat(categorical({'Wake'}),size(AS_WAKE_wSMI,3),1);...
%                         repmat(categorical({'Sleep'}),size(AS_SLEEP_wSMI,3),1)];
%                     ASgt2 = s == categorical({'Sleep'});
%                                         
%                     % apply model fit with optimal regularization to validation sets
%                     tmp_wSMI = [squeeze(AS_WAKE_wSMI(ich,jch,:)); squeeze(AS_SLEEP_wSMI(ich,jch,:))];
%                     ASfit = glmval(coef,tmp_wSMI,'logit','Constant','on');
%                     
%                     [~,~,~,AUC_AS] = perfcurve(ASgt2,ASfit,1);
%                     AS_wSMI_AUC(ich,jch) = AUC_AS;
%                     
%                     %%%% Find AUC for dwPLI %%%%
%                     
%                     % apply model fit with optimal regularization to validation sets
%                     tmp_dwPLI = [squeeze(AS_WAKE_dwPLI(ich,jch,:)); squeeze(AS_SLEEP_dwPLI(ich,jch,:))];
%                     ASfit = glmval(coef,tmp_dwPLI,'logit','Constant','on');
%                     
%                     [~,~,~,AUC_AS] = perfcurve(ASgt2,ASfit,1);
%                     AS_dwPLI_AUC(ich,jch) = AUC_AS;
%                     
%                     % Dup15q
%                     DS_SLEEP_wSMI(ich,jch,:)  = wSMI_PC_COEF'*squeeze(DSSleepwSMI(ich,jch,:,:));
%                     DS_WAKE_wSMI(ich,jch,:)   = wSMI_PC_COEF'*squeeze(DSWakewSMI(ich,jch,:,:));
%                     DS_SLEEP_dwPLI(ich,jch,:) = dwPLI_PC_COEF'*squeeze(DSSleepdwPLI(ich,jch,:,:));
%                     DS_WAKE_dwPLI(ich,jch,:)  = dwPLI_PC_COEF'*squeeze(DSWakedwPLI(ich,jch,:,:));
%                     % Define the IV (sleep)
%                     s = [repmat(categorical({'Wake'}),size(DS_WAKE_wSMI,3),1);...
%                         repmat(categorical({'Sleep'}),size(DS_SLEEP_wSMI,3),1)];
%                     DSgt2 = s == categorical({'Sleep'});
%                     
%                     %%%% Find AUC for wSMI %%%%
%                                         
%                     % apply model fit with optimal regularization to validation sets
%                     tmp_wSMI = [squeeze(DS_WAKE_wSMI(ich,jch,:)); squeeze(DS_SLEEP_wSMI(ich,jch,:))];
%                     DSfit = glmval(coef,tmp_wSMI,'logit','Constant','on');
%                     
%                     [~,~,~,AUC_DS] = perfcurve(DSgt2,DSfit,1);
%                     DS_wSMI_AUC(ich,jch) = AUC_DS;
%                     
%                     %%%% Find AUC for dwPLI %%%%
%                                         
%                     % apply model fit with optimal regularization to validation sets
%                     tmp_dwPLI = [squeeze(DS_WAKE_dwPLI(ich,jch,:)); squeeze(DS_SLEEP_dwPLI(ich,jch,:))];
%                     DSfit = glmval(coef,tmp_dwPLI,'logit','Constant','on');
%                     
%                     [~,~,~,AUC_DS] = perfcurve(DSgt2,DSfit,1);
%                     DS_dwPLI_AUC(ich,jch) = AUC_DS;
%                 end
%             end
%         end  
%         
%         
%         %% Plot figures for connecitivity AUCs
%         
%         TD_conn_AUC = nanmean(cat(3,TD_wSMI_AUC,rot90(TD_dwPLI_AUC,2)),3);
%         myfigure2
%         pcolor(TD_conn_AUC), colormap parula, caxis([0 1])
%         title('TD PC_{1} AUC (top triangle: wSMI; bottom triangle: dwPLI)','fontsize',18)
%         xticks([1:nchan])
%         xticklabels(lay.label)
%         yticks([1:nchan])
%         yticklabels(lay.label)
%         mycolorbar
%         makefighandsome, axis square
%         print('-dpng','./Figures/TD_PC1_AUC_connectivity.png')
%         print('-dsvg','./Figures/TD_PC1_AUC_connectivity.svg')
%         
%         AS_conn_AUC = nanmean(cat(3,AS_wSMI_AUC,rot90(AS_dwPLI_AUC,2)),3);
%         myfigure2
%         pcolor(AS_conn_AUC), colormap parula, caxis([0 1])
%         title('AS PC_{1} AUC (top triangle: wSMI; bottom triangle: dwPLI)','fontsize',18)
%         xticks([1:nchan])
%         xticklabels(lay.label)
%         yticks([1:nchan])
%         yticklabels(lay.label)
%         mycolorbar
%         makefighandsome, axis square
%         print('-dpng','./Figures/AS_PC1_AUC_connectivity.png')
%         print('-dsvg','./Figures/AS_PC1_AUC_connectivity.svg')
%         
%         DS_conn_AUC = nanmean(cat(3,DS_wSMI_AUC,rot90(DS_dwPLI_AUC,2)),3);
%         myfigure2
%         pcolor(DS_conn_AUC), colormap parula, caxis([0 1])
%         title('DS PC_{1} AUC (top triangle: wSMI; bottom triangle: dwPLI)','fontsize',18)
%         xticks([1:nchan])
%         xticklabels(lay.label)
%         yticks([1:nchan])
%         yticklabels(lay.label)
%         mycolorbar
%         makefighandsome, axis square
%         print('-dpng','./Figures/DS_PC1_AUC_connectivity.png')
%         print('-dsvg','./Figures/DS_PC1_AUC_connectivity.svg')        



        %% Look at univariate AUCs to produce figure like that from Sitt's group
        npc = 1;
        varnames2 = T.Properties.VariableNames(max(find(strcmp(vartype,'IV')))+1:end);
        vartype2 = [vartype T.Properties.VariableNames(contains(T.Properties.VariableNames,'PC'))];
        AUC = nan(3,length(varnames2));
        ACC = nan(3,length(varnames2));
        SENS = nan(3,length(varnames2));
        SPEC = nan(3,length(varnames2));
        PPV = nan(3,length(varnames2));
        
               
        for ivar = frstvar:length(T.Properties.VariableNames)
            ivar
            % Choose regularization parameter using 10-fold cross validation on all
            % training data (TD)
            
            eval(sprintf('xvar = T.%s;',T.Properties.VariableNames{ivar}))
            
            % use UNREGULARIZED logistic regression
            glm = fitglm(xvar(ASidx),T.Conscious(ASidx),'Distribution','binomial','Link','logit')
            coef = glm.Coefficients{:,1};
            
            % apply model fit with optimal regularization to validation sets
            TDfit = glmval(coef,xvar(TDidx),'logit','Constant','on');
            ASfit = glmval(coef,xvar(ASidx),'logit','Constant','on');
            DSfit = glmval(coef,xvar(DSidx),'logit','Constant','on');
            
            
            [X_TD,Y_TD,T_TD,AUC_TD,opt_TD] = perfcurve(TDgt,TDfit,1);
            thresh = T_TD(intersect(find(X_TD == opt_TD(1)),find(Y_TD == opt_TD(2))));
            acc_TD = sum(TDgt == (TDfit >= thresh))/length(TDfit)*100;
            fprintf('     TD training set accuracy = %1.2f%%, AUC = %1.2f%%\n',acc_TD,AUC_TD*100)
            TP = sum(TDgt & (TDfit >= thresh)); % true positives
            FP = sum(~TDgt & (TDfit >= thresh)); % false positives
            TN = sum(~TDgt & (TDfit < thresh)); % true negatives
            FN = sum(TDgt & (TDfit < thresh)); % false negatives
            assert(acc_TD==((TP+TN)/(TP+FP+TN+FN))*100,'Accuracy is wrong')
            
            if TP == 0 && FN == 0
                sens_TD = 0;
            else
                sens_TD = TP/(TP+FN);
            end
            
            if TN == 0 && FP == 0
                spec_TD = 0;
            else
                spec_TD = TN/(TN+FP);
            end
            
            if TP == 0 && FP == 0
                ppv_TD = 0;
            else
                ppv_TD = TP/(TP+FP);
            end
            
            [X_AS,Y_AS,T_AS,AUC_AS,opt_AS] = perfcurve(ASgt,ASfit,1);
            thresh = T_AS(intersect(find(X_AS == opt_AS(1)),find(Y_AS == opt_AS(2))));
            acc_AS = sum(ASgt == (ASfit >= thresh))/length(ASfit)*100;
            fprintf('    AS validation set accuracy = %1.2f%%, AUC = %1.2f%%\n',acc_AS,AUC_AS*100)
            TP = sum(ASgt & (ASfit >= thresh)); % true positives
            FP = sum(~ASgt & (ASfit >= thresh)); % false positives
            TN = sum(~ASgt & (ASfit < thresh)); % true negatives
            FN = sum(ASgt & (ASfit < thresh)); % false negatives
            assert(acc_AS==((TP+TN)/(TP+FP+TN+FN))*100,'Accuracy is wrong')
            
            if TP == 0 && FN == 0
                sens_AS = 0;
            else
                sens_AS = TP/(TP+FN);
            end
            
            if TN == 0 && FP == 0
                spec_AS = 0;
            else
                spec_AS = TN/(TN+FP);
            end
            
            if TP == 0 && FP == 0
                ppv_AS = 0;
            else
                ppv_AS = TP/(TP+FP);
            end
            
            [X_DS,Y_DS,T_DS,AUC_DS,opt_DS] = perfcurve(DSgt,DSfit,1);
            thresh = T_DS(intersect(find(X_DS == opt_DS(1)),find(Y_DS == opt_DS(2))));
            acc_DS = sum(DSgt == (DSfit >= thresh))/length(DSfit)*100;
            fprintf('    Dup15q validation set accuracy = %1.2f%%, AUC = %1.2f%%\n',acc_DS,AUC_DS*100)
            fprintf('\n')
            TP = sum(DSgt & (DSfit >= thresh)); % true positives
            FP = sum(~DSgt & (DSfit >= thresh)); % false positives
            TN = sum(~DSgt & (DSfit < thresh)); % true negatives
            FN = sum(DSgt & (DSfit < thresh)); % false negatives
            assert(acc_DS==((TP+TN)/(TP+FP+TN+FN))*100,'Accuracy is wrong')
            
            if TP == 0 && FN == 0
                sens_DS = 0;
            else
                sens_DS = TP/(TP+FN);
            end
            
            if TN == 0 && FP == 0
                spec_DS = 0;
            else
                spec_DS = TN/(TN+FP);
            end
            
            if TP == 0 && FP == 0
                ppv_DS = 0;
            else
                ppv_DS = TP/(TP+FP);
            end
            
            % Area under the ROC curve for each group
            AUC(1,ivar) = AUC_TD;
            AUC(2,ivar) = AUC_AS;
            AUC(3,ivar) = AUC_DS;
            
            ACC(1,ivar) = acc_TD/100;
            ACC(2,ivar) = acc_AS/100;
            ACC(3,ivar) = acc_DS/100;
            
            % Sensitivity for each group
            SENS(1,ivar) = sens_TD;
            SENS(2,ivar) = sens_AS;
            SENS(3,ivar) = sens_DS;
            
            % Specificity for each group
            SPEC(1,ivar) = spec_TD;
            SPEC(2,ivar) = spec_AS;
            SPEC(3,ivar) = spec_DS;
            
            % Precision for each group
            PPV(1,ivar) = ppv_TD;
            PPV(2,ivar) = ppv_AS;
            PPV(3,ivar) = ppv_DS; 
        end
        

        Tall = table(T.Properties.VariableNames(frstvar:end)',AUC(1,frstvar:end)',AUC(2,frstvar:end)',AUC(3,frstvar:end)',...
            ACC(1,frstvar:end)',ACC(2,frstvar:end)',ACC(3,frstvar:end)',...
            PPV(1,frstvar:end)',PPV(2,frstvar:end)',PPV(3,frstvar:end)',...
            SENS(1,frstvar:end)',SENS(2,frstvar:end)',SENS(3,frstvar:end)',...
            SPEC(1,frstvar:end)',SPEC(2,frstvar:end)',SPEC(3,frstvar:end)','VariableNames',...
            {'Vars','NT AUC','AS AUC','DS AUC',...
            'NT ACC','AS ACC','DS ACC',...
            'NT PPV','AS PPV','DS PPV',...
            'NT TPR','AS TPR','DS TPR',...
            'NT TNR','AS TNR','DS TNR'})
        
        writetable(Tall,'./TrainOnTD/AllVariableScoresTD.csv')
        
       
        
        %% Show classification stats for all variables
        tmp = properties(Tall);
        vartype3 = vartype2(frstvar:end);
        reorder = [find(contains(vartype3,'scEntropy')) find(contains(vartype3,'fcEntropy')) ...
            find(contains(vartype3,'scSpectralA')) find(contains(vartype3,'scSpectralR')) ...
            find(contains(vartype3,'fcSpectral'))];
        tmp2 = Tall.Vars(reorder);
        varlabels = replace(tmp2,'_',' ');
        varlabels = replace(varlabels,'slow','s');
        varlabels = replace(varlabels,'delta','?');
        varlabels = replace(varlabels,'theta','?');
        varlabels = replace(varlabels,'alpha','?-?');
        varlabels = replace(varlabels,'beta','?');
        
        
        myfigure2
        % do it this way because pcolor function needs buffer on edges
        plttbl = ones(size(Tall,2),size(Tall,1)+1);
        plttbl(1:size(Tall,2)-1,1:size(Tall,1)) = flipud(Tall{reorder,2:end}');
        pcolor(plttbl)
        colormap parula
        caxis([0 1])
        axis([1 size(Tall,1)+1 1 size(Tall,2)])
        xticks(1.5:1:size(Tall,1)+1)
        xticklabels(varlabels)
        xtickangle(45)
        yticks(1.5:1:size(Tall,2))
        yticklabels(flipud(tmp(2:size(Tall,2))))
        box off
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', 16)
        set(xAX,'color','k')
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', 20)
        set(yAX,'color','k')
        set(gca, 'TickDir', 'out')
        set(gcf,'color','w')
        set(gca,'Layer','top')
        axis normal
        print('-dpng','./TrainOnTD/Figures/AllVariableScores.png')
        print('-dsvg','./TrainOnTD/Figures/AllVariableScores.svg')
        
         %% Plot precision vs recall (specificity)
        
        % Include indices for the IVs below 
        
        PCidx = find(contains(vartype2,'PC'));
%         
%         rmvidx = contains(vartype2,'IV');
%         vartype2(rmvidx) = [];
       
                
        g = {'TD','AS','Dup15q'};
        
        idxLMM = abs(TDsleepES(1,1:end)) >= betathresh; % features selected using LMMs
        idxPC = [false(1,length(idxLMM)) true(1,5)]; % indices of PCs
        idxLMM = [idxLMM false(1,5)]; % pad with false (for PCs)
        size(idxLMM)
        TDclr = [0, 0.4470, 0.7410];
        CXclr = [0.8500, 0.3250, 0.0980];
        MIclr = [0.9290, 0.6940, 0.1250];
        RSclr = [0.4940, 0.1840, 0.5560];
        PIclr = [0.4660, 0.6740, 0.1880];
                
        for igrp = 1:3
            myfigure2
            scatter(PPV(igrp,contains(vartype2,'scEntropy')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'scEntropy')&~idxLMM&~idxPC),135,CXclr,'fill')
            scatter(PPV(igrp,contains(vartype2,'fcEntropy')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'fcEntropy')&~idxLMM&~idxPC),135,MIclr,'fill')
            scatter(PPV(igrp,contains(vartype2,'fcSpectral')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'fcSpectral')&~idxLMM&~idxPC),135,PIclr,'fill')
            scatter(PPV(igrp,contains(vartype2,'scSpectralA')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'scSpectralA')&~idxLMM&~idxPC),135,TDclr,'fill')
            scatter(PPV(igrp,contains(vartype2,'scSpectralR')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'scSpectralR')&~idxLMM&~idxPC),135,RSclr,'fill')
            
            buf = -0.01;
            plot([buf 1],[buf 1],'k--')

            legend({'scEntropy','fcEntropy','fcSpectral','scSpectralA',...
                'scSpectralR','Chance performance'},'autoupdate','off',...
                'location','northeastoutside','fontsize',16)

            % features selected by LMMs
            scatter(PPV(igrp,contains(vartype2,'scEntropy')&idxLMM),SENS(igrp,contains(vartype2,'scEntropy')&idxLMM),150,CXclr,'d','fill')
            scatter(PPV(igrp,contains(vartype2,'fcEntropy')&idxLMM),SENS(igrp,contains(vartype2,'fcEntropy')&idxLMM),150,MIclr,'d','fill')
            scatter(PPV(igrp,contains(vartype2,'fcSpectral')&idxLMM),SENS(igrp,contains(vartype2,'fcSpectral')&idxLMM),150,PIclr,'d','fill')
            scatter(PPV(igrp,contains(vartype2,'scSpectralA')&idxLMM),SENS(igrp,contains(vartype2,'scSpectralA')&idxLMM),150,TDclr,'d','fill')
            scatter(PPV(igrp,contains(vartype2,'scSpectralR')&idxLMM),SENS(igrp,contains(vartype2,'scSpectralR')&idxLMM),150,RSclr,'d','fill')
            
            % First principal components
            scatter(PPV(igrp,contains(vartype2,'scEntropy')&idxPC),SENS(igrp,contains(vartype2,'scEntropy')&idxPC),200,CXclr,'p','fill')
            scatter(PPV(igrp,contains(vartype2,'fcEntropy')&idxPC),SENS(igrp,contains(vartype2,'fcEntropy')&idxPC),200,MIclr,'p','fill')
            scatter(PPV(igrp,contains(vartype2,'fcSpectral')&idxPC),SENS(igrp,contains(vartype2,'fcSpectral')&idxPC),200,PIclr,'p','fill')
            scatter(PPV(igrp,contains(vartype2,'scSpectralA')&idxPC),SENS(igrp,contains(vartype2,'scSpectralA')&idxPC),200,TDclr,'p','fill')
            scatter(PPV(igrp,contains(vartype2,'scSpectralR')&idxPC),SENS(igrp,contains(vartype2,'scSpectralR')&idxPC),200,RSclr,'p','fill')
            
            legend boxoff
            
            title(g{igrp},'fontsize',18)
            xlabel('Precision (%)')
            ylabel('Recall (%)')
            xticks(0:0.2:1)
            xticklabels(0:20:100)
            yticks(0:0.2:1)
            yticklabels(0:20:100)
            axis([buf 1 buf 1])
            makefigpretty
            
            DT = delaunayTriangulation(PPV(igrp,contains(vartype2,'scEntropy'))',...
                SENS(igrp,contains(vartype2,'scEntropy'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),CXclr)
            
            DT = delaunayTriangulation(PPV(igrp,contains(vartype2,'fcEntropy'))',...
                SENS(igrp,contains(vartype2,'fcEntropy'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),MIclr)
            
            DT = delaunayTriangulation(PPV(igrp,contains(vartype2,'scSpectralA'))',...
                SENS(igrp,contains(vartype2,'scSpectralA'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),TDclr)
            
            DT = delaunayTriangulation(PPV(igrp,contains(vartype2,'scSpectralR'))',...
                SENS(igrp,contains(vartype2,'scSpectralR'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),RSclr)
            
            DT = delaunayTriangulation(PPV(igrp,contains(vartype2,'fcSpectral'))',...
                SENS(igrp,contains(vartype2,'fcSpectral'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),PIclr)
                      
            axis([0 1 0 1])
            axis square
            
            print('-dsvg',sprintf('./TrainOnTD/Figures/ppv_vs_recall_scatter_plot_%s.svg',g{igrp}))
            print('-dpng',sprintf('./TrainOnTD/Figures/ppv_vs_recall_scatter_plot_%s.png',g{igrp}))
        end
        %% Plot precision vs recall upper quadrant
                
        % Make names better suited for these labels
        varnames3 = {'mMSE','LZc','CTWc','PE8','PE16','PE32','PE64','PE128',...
            'SR','SR16','SR32','SR64','SR128','LR8','LR16','LR32','LR64','LR128',...
            'sA','\delta1A','\delta2A','\thetaA','\alpha-\sigmaA','\betaA','sR','\delta1R',...
            '\delta2R','\thetaR','\alpha-\sigmaR','\betaR',...
            'SRs','SR\delta1','SR\delta2','SR\theta','SR\alpha-\sigma',...
            'SR\beta','LRs','LR\delta1','LR\delta2','LR\theta','LR\alpha-LR\sigma','LR\beta'};
        
        TDclr = [0, 0.4470, 0.7410];
        CXclr = [0.8500, 0.3250, 0.0980];
        MIclr = [0.9290, 0.6940, 0.1250];
        RSclr = [0.4940, 0.1840, 0.5560];
        PIclr = [0.4660, 0.6740, 0.1880];
        
        for igrp = 1:3
            myfigure2
            scatter(PPV(igrp,contains(vartype2,'scEntropy')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'scEntropy')&~idxLMM&~idxPC),135,CXclr,'fill')
            scatter(PPV(igrp,contains(vartype2,'fcEntropy')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'fcEntropy')&~idxLMM&~idxPC),135,MIclr,'fill')
            scatter(PPV(igrp,contains(vartype2,'fcSpectral')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'fcSpectral')&~idxLMM&~idxPC),135,PIclr,'fill')
            scatter(PPV(igrp,contains(vartype2,'scSpectralA')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'scSpectralA')&~idxLMM&~idxPC),135,TDclr,'fill')
            scatter(PPV(igrp,contains(vartype2,'scSpectralR')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'scSpectralR')&~idxLMM&~idxPC),135,RSclr,'fill')
            
            buf = -0.01;
            plot([buf 1],[buf 1],'k--')

            legend({'scEntropy','fcEntropy','fcSpectral','scSpectralA',...
                'scSpectralR','Chance performance'},'autoupdate','off',...
                'location','northeastoutside','fontsize',16)

            % features selected by LMMs
            scatter(PPV(igrp,contains(vartype2,'scEntropy')&idxLMM),SENS(igrp,contains(vartype2,'scEntropy')&idxLMM),150,CXclr,'d','fill')
            scatter(PPV(igrp,contains(vartype2,'fcEntropy')&idxLMM),SENS(igrp,contains(vartype2,'fcEntropy')&idxLMM),150,MIclr,'d','fill')
            scatter(PPV(igrp,contains(vartype2,'fcSpectral')&idxLMM),SENS(igrp,contains(vartype2,'fcSpectral')&idxLMM),150,PIclr,'d','fill')
            scatter(PPV(igrp,contains(vartype2,'scSpectralA')&idxLMM),SENS(igrp,contains(vartype2,'scSpectralA')&idxLMM),150,TDclr,'d','fill')
            scatter(PPV(igrp,contains(vartype2,'scSpectralR')&idxLMM),SENS(igrp,contains(vartype2,'scSpectralR')&idxLMM),150,RSclr,'d','fill')
            
            % First principal components
            scatter(PPV(igrp,contains(vartype2,'scEntropy')&idxPC),SENS(igrp,contains(vartype2,'scEntropy')&idxPC),200,CXclr,'p','fill')
            scatter(PPV(igrp,contains(vartype2,'fcEntropy')&idxPC),SENS(igrp,contains(vartype2,'fcEntropy')&idxPC),200,MIclr,'p','fill')
            scatter(PPV(igrp,contains(vartype2,'fcSpectral')&idxPC),SENS(igrp,contains(vartype2,'fcSpectral')&idxPC),200,PIclr,'p','fill')
            scatter(PPV(igrp,contains(vartype2,'scSpectralA')&idxPC),SENS(igrp,contains(vartype2,'scSpectralA')&idxPC),200,TDclr,'p','fill')
            scatter(PPV(igrp,contains(vartype2,'scSpectralR')&idxPC),SENS(igrp,contains(vartype2,'scSpectralR')&idxPC),200,RSclr,'p','fill')
            
            legend boxoff
            title(g{igrp},'fontsize',18)
            xlabel('Precision (%)')
            ylabel('Recall (%)')
            xticks(0:0.1:1)
            xticklabels(0:10:100)
            yticks(0:0.1:1)
            yticklabels(0:10:100)
            axis([0.5 1 0.5 1])
            makefigpretty
            
            %jitter = randn(2,length(T.Properties.VariableNames))./100;
            jitter = zeros(2,length(T.Properties.VariableNames));
            for ivar = frstvar:length(T.Properties.VariableNames)
                idx = ivar-frstvar+1;
                switch vartype2{ivar}
                    case 'scEntropy'
                        text(PPV(igrp,ivar)+jitter(1,ivar),SENS(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',CXclr,'FontSize',14)
                    case 'fcEntropy'
                        text(PPV(igrp,ivar)+jitter(1,ivar),SENS(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',MIclr,'FontSize',14)
                    case 'scSpectralA'
                        text(PPV(igrp,ivar)+jitter(1,ivar),SENS(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',TDclr,'FontSize',14)
                    case 'scSpectralR'
                        text(PPV(igrp,ivar)+jitter(1,ivar),SENS(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',RSclr,'FontSize',14)
                    case 'fcSpectral'
                        text(PPV(igrp,ivar)+jitter(1,ivar),SENS(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',PIclr,'FontSize',14)
                end
            end
                                
            % vertical and horizontal lines for AUC thresholds
            tcks = 0:0.1:1;

            for i = 1:length(tcks)      
                plot(ones(1,100).*tcks(i),linspace(buf,1,100),'k:')
                plot(linspace(buf,1,100),ones(1,100).*tcks(i),'k:')
            end
            
            %plot2svg(sprintf('./TrainOnTD/Figures/sens_vs_spec_scatter_plot_upper_quadrant_%s.svg',g{igrp}));
            print('-dsvg',sprintf('./TrainOnTD/Figures/ppv_vs_recall_scatter_plot_upper_quadrant_%s.svg',g{igrp}))
            print('-dpng',sprintf('./TrainOnTD/Figures/ppv_vs_recall_scatter_plot_upper_quadrant_%s.png',g{igrp}))
        end

        
        %% Plot sensitivity vs specificity
        
        % Include indices for the IVs below 
        
        PCidx = find(contains(vartype2,'PC'));
        
%         rmvidx = contains(vartype2,'IV');
%         vartype2(rmvidx) = [];
       
                
        g = {'TD','AS','Dup15q'};
        
        idxLMM = abs(TDsleepES(1,1:end)) >= 0.5; % features selected using LMMs
        idxPC = [false(1,length(idxLMM)) true(1,5)]; % indices of PCs
        idxLMM = [idxLMM false(1,5)]; % pad with false (for PCs)
        size(idxLMM)
        TDclr = [0, 0.4470, 0.7410];
        CXclr = [0.8500, 0.3250, 0.0980];
        MIclr = [0.9290, 0.6940, 0.1250];
        RSclr = [0.4940, 0.1840, 0.5560];
        PIclr = [0.4660, 0.6740, 0.1880];
                
        for igrp = 1:3
            myfigure2
            scatter(SENS(igrp,contains(vartype2,'scEntropy')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'scEntropy')&~idxLMM&~idxPC),135,CXclr,'fill')
            scatter(SENS(igrp,contains(vartype2,'fcEntropy')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'fcEntropy')&~idxLMM&~idxPC),135,MIclr,'fill')
            scatter(SENS(igrp,contains(vartype2,'fcSpectral')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'fcSpectral')&~idxLMM&~idxPC),135,PIclr,'fill')
            scatter(SENS(igrp,contains(vartype2,'scSpectralA')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'scSpectralA')&~idxLMM&~idxPC),135,TDclr,'fill')
            scatter(SENS(igrp,contains(vartype2,'scSpectralR')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'scSpectralR')&~idxLMM&~idxPC),135,RSclr,'fill')
            
            buf = -0.01;
            plot([buf 1],[buf 1],'k--')

            legend({'scEntropy','fcEntropy','fcSpectral','scSpectralA',...
                'scSpectralR','Chance performance'},'autoupdate','off',...
                'location','northeastoutside','fontsize',16)

            % features selected by LMMs
            scatter(SENS(igrp,contains(vartype2,'scEntropy')&idxLMM),SPEC(igrp,contains(vartype2,'scEntropy')&idxLMM),150,CXclr,'d','fill')
            scatter(SENS(igrp,contains(vartype2,'fcEntropy')&idxLMM),SPEC(igrp,contains(vartype2,'fcEntropy')&idxLMM),150,MIclr,'d','fill')
            scatter(SENS(igrp,contains(vartype2,'fcSpectral')&idxLMM),SPEC(igrp,contains(vartype2,'fcSpectral')&idxLMM),150,PIclr,'d','fill')
            scatter(SENS(igrp,contains(vartype2,'scSpectralA')&idxLMM),SPEC(igrp,contains(vartype2,'scSpectralA')&idxLMM),150,TDclr,'d','fill')
            scatter(SENS(igrp,contains(vartype2,'scSpectralR')&idxLMM),SPEC(igrp,contains(vartype2,'scSpectralR')&idxLMM),150,RSclr,'d','fill')
            
            % First principal components
            scatter(SENS(igrp,contains(vartype2,'scEntropy')&idxPC),SPEC(igrp,contains(vartype2,'scEntropy')&idxPC),200,CXclr,'p','fill')
            scatter(SENS(igrp,contains(vartype2,'fcEntropy')&idxPC),SPEC(igrp,contains(vartype2,'fcEntropy')&idxPC),200,MIclr,'p','fill')
            scatter(SENS(igrp,contains(vartype2,'fcSpectral')&idxPC),SPEC(igrp,contains(vartype2,'fcSpectral')&idxPC),200,PIclr,'p','fill')
            scatter(SENS(igrp,contains(vartype2,'scSpectralA')&idxPC),SPEC(igrp,contains(vartype2,'scSpectralA')&idxPC),200,TDclr,'p','fill')
            scatter(SENS(igrp,contains(vartype2,'scSpectralR')&idxPC),SPEC(igrp,contains(vartype2,'scSpectralR')&idxPC),200,RSclr,'p','fill')
            
            legend boxoff
            
            title(g{igrp},'fontsize',18)
            xlabel('Sensitivity (%)')
            ylabel('Specificity (%)')
            xticks(0:0.2:1)
            xticklabels(0:20:100)
            yticks(0:0.2:1)
            yticklabels(0:20:100)
            axis([buf 1 buf 1])
            makefigpretty
            
            DT = delaunayTriangulation(SENS(igrp,contains(vartype2,'scEntropy'))',...
                SPEC(igrp,contains(vartype2,'scEntropy'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),CXclr)
            
            DT = delaunayTriangulation(SENS(igrp,contains(vartype2,'fcEntropy'))',...
                SPEC(igrp,contains(vartype2,'fcEntropy'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),MIclr)
            
            DT = delaunayTriangulation(SENS(igrp,contains(vartype2,'scSpectralA'))',...
                SPEC(igrp,contains(vartype2,'scSpectralA'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),TDclr)
            
            DT = delaunayTriangulation(SENS(igrp,contains(vartype2,'scSpectralR'))',...
                SPEC(igrp,contains(vartype2,'scSpectralR'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),RSclr)
            
            DT = delaunayTriangulation(SENS(igrp,contains(vartype2,'fcSpectral'))',...
                SPEC(igrp,contains(vartype2,'fcSpectral'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),PIclr)
                      
            axis([0 1 0 1])
            axis square
            
            print('-dsvg',sprintf('./TrainOnTD/Figures/sens_vs_spec_scatter_plot_%s.svg',g{igrp}))
            print('-dpng',sprintf('./TrainOnTD/Figures/sens_vs_spec_scatter_plot_%s.png',g{igrp}))
        end
        %% Plot sensitivity vs specificity upper quadrant
                
        % Make names better suited for these labels
        varnames3 = {'mMSE','LZc','CTWc','PE8','PE16','PE32','PE64','PE128',...
            'SR','SR16','SR32','SR64','SR128','LR8','LR16','LR32','LR64','LR128',...
            's','\delta1','\delta2','\theta','\alpha-\sigma','\beta','s','\delta1',...
            '\delta2','\theta','\alpha-\sigma','\beta',...
            'SRs','SR\delta1','SR\delta2','SR\theta','SR\alpha-\sigma',...
            'SR\beta','LRs','LR\delta1','LR\delta2','LR\theta','LR\alpha-LR\sigma','LR\beta'};
        
        TDclr = [0, 0.4470, 0.7410];
        CXclr = [0.8500, 0.3250, 0.0980];
        MIclr = [0.9290, 0.6940, 0.1250];
        RSclr = [0.4940, 0.1840, 0.5560];
        PIclr = [0.4660, 0.6740, 0.1880];
        
        for igrp = 1:3
            myfigure2
            scatter(SENS(igrp,contains(vartype2,'scEntropy')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'scEntropy')&~idxLMM&~idxPC),135,CXclr,'fill')
            scatter(SENS(igrp,contains(vartype2,'fcEntropy')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'fcEntropy')&~idxLMM&~idxPC),135,MIclr,'fill')
            scatter(SENS(igrp,contains(vartype2,'fcSpectral')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'fcSpectral')&~idxLMM&~idxPC),135,PIclr,'fill')
            scatter(SENS(igrp,contains(vartype2,'scSpectralA')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'scSpectralA')&~idxLMM&~idxPC),135,TDclr,'fill')
            scatter(SENS(igrp,contains(vartype2,'scSpectralR')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'scSpectralR')&~idxLMM&~idxPC),135,RSclr,'fill')
            
            buf = -0.01;
            plot([buf 1],[buf 1],'k--')

            legend({'scEntropy','fcEntropy','fcSpectral','scSpectralA',...
                'scSpectralR','Chance performance'},'autoupdate','off',...
                'location','northeastoutside','fontsize',16)

            % features selected by LMMs
            scatter(SENS(igrp,contains(vartype2,'scEntropy')&idxLMM),SPEC(igrp,contains(vartype2,'scEntropy')&idxLMM),150,CXclr,'d','fill')
            scatter(SENS(igrp,contains(vartype2,'fcEntropy')&idxLMM),SPEC(igrp,contains(vartype2,'fcEntropy')&idxLMM),150,MIclr,'d','fill')
            scatter(SENS(igrp,contains(vartype2,'fcSpectral')&idxLMM),SPEC(igrp,contains(vartype2,'fcSpectral')&idxLMM),150,PIclr,'d','fill')
            scatter(SENS(igrp,contains(vartype2,'scSpectralA')&idxLMM),SPEC(igrp,contains(vartype2,'scSpectralA')&idxLMM),150,TDclr,'d','fill')
            scatter(SENS(igrp,contains(vartype2,'scSpectralR')&idxLMM),SPEC(igrp,contains(vartype2,'scSpectralR')&idxLMM),150,RSclr,'d','fill')
            
            % First principal components
            scatter(SENS(igrp,contains(vartype2,'scEntropy')&idxPC),SPEC(igrp,contains(vartype2,'scEntropy')&idxPC),200,CXclr,'p','fill')
            scatter(SENS(igrp,contains(vartype2,'fcEntropy')&idxPC),SPEC(igrp,contains(vartype2,'fcEntropy')&idxPC),200,MIclr,'p','fill')
            scatter(SENS(igrp,contains(vartype2,'fcSpectral')&idxPC),SPEC(igrp,contains(vartype2,'fcSpectral')&idxPC),200,PIclr,'p','fill')
            scatter(SENS(igrp,contains(vartype2,'scSpectralA')&idxPC),SPEC(igrp,contains(vartype2,'scSpectralA')&idxPC),200,TDclr,'p','fill')
            scatter(SENS(igrp,contains(vartype2,'scSpectralR')&idxPC),SPEC(igrp,contains(vartype2,'scSpectralR')&idxPC),200,RSclr,'p','fill')
            
            legend boxoff
            title(g{igrp},'fontsize',18)
            xlabel('Sensitivity (%)')
            ylabel('Specificity (%)')
            xticks(0:0.1:1)
            xticklabels(0:10:100)
            yticks(0:0.1:1)
            yticklabels(0:10:100)
            axis([0.5 1 0.5 1])
            makefigpretty
            
            %jitter = randn(2,length(T.Properties.VariableNames))./100;
            jitter = zeros(2,length(T.Properties.VariableNames));
            for ivar = frstvar:length(T.Properties.VariableNames)
                idx = ivar-frstvar+1;
                switch vartype2{ivar}
                    case 'scEntropy'
                        text(SENS(igrp,ivar)+jitter(1,ivar),SPEC(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',CXclr,'FontSize',14)
                    case 'fcEntropy'
                        text(SENS(igrp,ivar)+jitter(1,ivar),SPEC(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',MIclr,'FontSize',14)
                    case 'scSpectralA'
                        text(SENS(igrp,ivar)+jitter(1,ivar),SPEC(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',TDclr,'FontSize',14)
                    case 'scSpectralR'
                        text(SENS(igrp,ivar)+jitter(1,ivar),SPEC(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',RSclr,'FontSize',14)
                    case 'fcSpectral'
                        text(SENS(igrp,ivar)+jitter(1,ivar),SPEC(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',PIclr,'FontSize',14)
                end
            end
                                
            % vertical and horizontal lines for AUC thresholds
            tcks = 0:0.1:1;

            for i = 1:length(tcks)      
                plot(ones(1,100).*tcks(i),linspace(buf,1,100),'k:')
                plot(linspace(buf,1,100),ones(1,100).*tcks(i),'k:')
            end
            
            %plot2svg(sprintf('./TrainOnTD/Figures/sens_vs_spec_scatter_plot_upper_quadrant_%s.svg',g{igrp}));
            print('-dsvg',sprintf('./TrainOnTD/Figures/sens_vs_spec_scatter_plot_upper_quadrant_%s.svg',g{igrp}))
            print('-dpng',sprintf('./TrainOnTD/Figures/sens_vs_spec_scatter_plot_upper_quadrant_%s.png',g{igrp}))
        end      
end
        
        
        