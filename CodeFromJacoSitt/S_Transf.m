function [sym count] = S_transf(X,cfg)
% ---------------------------------------------------------------------
% (c) Jacobo Sitt
% ---------------------------------------------------------------------
tic

chan_sel  = cfg.chan_sel;
data_sel  = cfg.data_sel;
taus   = cfg.taus;
kernel = cfg.kernel;

for tau = 1:size(taus,2)  %% 1: or 3:
    
    %%% filter
    
    filter_freq = cfg.sf/cfg.kernel/taus(tau);
    %filter_freq = cfg.sf/taus(tau)/cfg.kernel;
    
    %disp(['Filtering @ ' num2str(filter_freq) ' Hz'])
    
    nTrials = size(X,3);
    %Xf      = zeros(size(X));
    
    %%% select data
    % JF edit #1 06.10.21 -- try catch block to allow 2D data
    % JF edit #2 09.02.21 -- prune data first so that filtering goes faster
    try
        X = X(chan_sel,data_sel,:);
    catch
        X = X(chan_sel,data_sel);
    end
    Xf      = nan(size(X));
   
    for nT =1:nTrials
        Xf(:,:,nT) = ft_preproc_lowpassfilter(squeeze(X(:,:,nT)),cfg.sf,filter_freq );
    end
    
    
%     %%% select data
%     % JF edit 06.10.21 -- try catch block to allow 2D data
%     try
%         Xf = Xf(chan_sel,data_sel,:);
%     catch
%         Xf = Xf(chan_sel,data_sel);
%     end
    
    %%%%
    %%% PE calculation
    
    %disp('PE calculation')
    
    %%% initilize variables
    sym{tau} = zeros(size(Xf,1),size(Xf,2)-(kernel-1)*taus(tau),size(Xf,3));
    count{tau} = zeros(size(Xf,1),factorial(kernel),size(Xf,3));
    
    
    for trial = 1:size(Xf,3)
        
        %disp(['trial ' num2str(trial) ' of ' num2str(size(Xf,3))])
        
        [sym{tau}(:,:,trial) count{tau}(:,:,trial) c] = PE_paralel(squeeze(Xf(:,:,trial))',kernel,taus(tau));
        
        count{tau}(:,:,trial) = count{tau}(:,:,trial)/size(sym{tau},2); %%% probabilistic format
        
    end
end


end




