clearvars
rng(45783) % seed added 02.26.19
OS = computer; % detect opperating system
dbstop if error
new_fsample = 125; % see King et al. 2013 supplemental results
ft_defaults

switch OS
    case 'GLNXA64'
        pth = './2021_analysis/postprocessed/Dup15q/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
    case 'MACI64'
        pth = './2021_analysis/postprocessed/Dup15q/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
    case 'PCWIN64'
        pth = './2021_analysis/postprocessed/Dup15q/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
end

% these are the files we already know we want to use
% load('select_files','select_files')

sleep = true; % asleep or awake?
surrogate = false; % scramble FFT phases?
compute_wSMI(sleep,surrogate,files,pth,new_fsample)

sleep = false; % asleep or awake?
surrogate = false; % scramble FFT phases?
compute_wSMI(sleep,surrogate,files,pth,new_fsample)
% 
% sleep = true; % asleep or awake?
% surrogate = true; % scramble FFT phases?
% compute_wSMI(sleep,surrogate,files,pth,match,new_fsample)
% 
% sleep = false; % asleep or awake?
% surrogate = true; % scramble FFT phases?
% compute_wSMI(sleep,surrogate,files,pth,match,new_fsample)
%%

function[] = compute_wSMI(sleep,surrogate,files,pth,new_fsample)

switch sleep
    case true
        DIRRESULT = './2021_analysis/Dup15q_output/wSMI/sleep/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end    
        ntrl = 2;
    case false       
        DIRRESULT = './2021_analysis/Dup15q_output/wSMI/wake/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end     
        ntrl = 1;
end

tic
for ifile = 1:length(files)
    %if logical(sum(ismember(select_files,files(ifile).name(1:end-4))))
    load(sprintf('%s%s',pth,files(ifile).name),'data')
    fprintf('Now loading %s \n',files(ifile).name)
    
    % check if this file had useable sleep
    
    % The second trial is always sleep, so this might not be necessary
    if ~isfield(data.cfg.dattype,'SLEEP_MONTI') || isempty(data.cfg.dattype.SLEEP_MONTI)
        fprintf('Skipping this one, no useable sleep ... \n')
        continue
    end
    
    %Get datamat var for freq transform
    datamat = data.trial{ntrl};
    
    % Do some santiy checks
    assert(size(datamat,1) == 19,'channel(s) missing!')
    % allow for small number of NaNs at beginning and end of file
    assert(sum(isnan(datamat(1:end)))/numel(datamat) < 0.01, 'Many NaNs in data!')
        
    % Set bad data to NaNs for purpose of taking freq transform
    badart = [data.cfg.dattype.bad; data.cfg.dattype.flash ];
    
    % add muscle to artifacts if it is a field
    if isfield(data.cfg.dattype,'muscle')
        badart = [badart; data.cfg.dattype.muscle];
    end
    
    badart(badart==0) = 1; % do this so no 0 indicies that give bugs
    
    for iart = 1:size(badart,1)
        idpt = badart(iart,1) - data.sampleinfo(ntrl,1) + 1;
        if idpt > 0 && idpt+diff(badart(iart,:)) <= size(datamat,2)
            datamat(:,idpt:idpt+diff(badart(iart,:))) = nan; 
        end
    end
    
    % determine amount of good data
    c = sum(~isnan(datamat(1,:)));
    good = c/data.fsample; % seconds of good data
    data.good = good;
    if good < 30
        fprintf('Skipping dataset, less than 30 seconds of good data ...\n')
        continue
    else
        display(strcat('Data has this many seconds of good data: ',num2str(round(good))))
    end
    
    percent = c/length(datamat)*100;
    durgood = c/data.fsample;
    
    % downsample to 125Hz prior to wSMI 
    datamat = single(resample(double(datamat)',new_fsample,data.fsample))';
    
    %% partion into n sec windows and compute wSMI
    
    nwin = 5;
    win_len = new_fsample*nwin; % n sec windows
    
    % BAND (tau): 32-80 Hz (? = 4 ms), 16-40 Hz (? = 8 ms), 8-20 Hz 
    % (? = 16 ms), 4-10 Hz (? = 32 ms), 2-5 Hz (? =  64 ms) and 1-2.5 Hz 
    %(? = 128 303ms) [from   Bourdillon et al. 2019 BioRxiv) 
    
    cfg = [];
    cfg.sf = new_fsample;
    cfg.verbose = false;
    cfg.kernel = 3;
    cfg.taus = round([8 16 32 64 128].*(new_fsample/1000)); % covert tau from ms to samples (e.g., 32 ms becomes 8 samples at fs = 250 Hz)
    cfg.taus_ms = [8 16 32 64 128];
    cfg.chan_sel = 1:size(datamat,1); % use all channels
    cfg.data_sel = 1:win_len; % use all samples

    gcount = 0;
    bcount = 0;
 
    % windowing parameters
    % 1.0 --> no overlap; 0.5 --> 50 percent
    window_shift=0.5; % do 50% overlap 
    n_shift = round(win_len*window_shift);
    for isection = 1:n_shift:size(datamat,2)-win_len+1
        section = double(datamat(:,isection:isection+win_len-1));
        if ~any(isnan(section(:))) > 0 && surrogate
            surdat = nan(size(section,1),sum(~isnan(mean(section)))); % allocation
            for ich = 1:size(section,1)
                surdat(ich,:) = surrogate_data(section(ich,:),new_fsample);
            end
            % overwrite data section with the surrogate data
            section = surdat;
        else
            %%% DO NOTHING if data are just nans or if we don't need to
            %%% make a surrogate anyway
        end
        
        if ~any(isnan(section(:))) % IF there are no NaNs in the window
            [cfg_out] = ro_wSMI(section,cfg);
            for itau = 1:length(cfg.taus)
                % create variable if it doesn't already exist
                %%% wSMI %%%
                if ~exist(sprintf('wSMI_%i',cfg.taus(itau)),'var')
                    eval(sprintf('wSMI_%i = [];',cfg.taus(itau)));
                end
                eval(sprintf('wSMI_%i = cat(3,cfg_out.wSMI{itau},wSMI_%i);',...
                    cfg.taus(itau),cfg.taus(itau)))
                %%% SMI (unweighted) %%%
                if ~exist(sprintf('SMI_%i',cfg.taus(itau)),'var')
                    eval(sprintf('SMI_%i = [];',cfg.taus(itau)));
                end
                eval(sprintf('SMI_%i = cat(3,cfg_out.SMI{itau},SMI_%i);',...
                    cfg.taus(itau),cfg.taus(itau)))
                %%% PE %%%
                if ~exist(sprintf('PE_%i',cfg.taus(itau)),'var')
                    eval(sprintf('PE_%i = [];',cfg.taus(itau)));
                end
                eval(sprintf('PE_%i = cat(2,cfg_out.PE{itau},PE_%i);',...
                    cfg.taus(itau),cfg.taus(itau)))
            end
            gcount = gcount + 1;
            %fprintf('%2.2f percent complete \n',(isection/(size(datamat,2)-win_len+1))*100)
            %toc
        else
            bcount = bcount + 1;
            %fprintf('Window contains NaNs, continuing to next window ...\n')
        end
    end
    
    wSMIout.cfg = cfg_out;
    wSMIout.window_shift = window_shift; % window overlap
    
    for itau = 1:length(cfg.taus)
        eval(sprintf('wSMIout.wSMI_%i = wSMI_%i;',cfg.taus(itau),cfg.taus(itau)));
        eval(sprintf('wSMIout.SMI_%i = SMI_%i;',cfg.taus(itau),cfg.taus(itau)));
        eval(sprintf('wSMIout.PE_%i = PE_%i;',cfg.taus(itau),cfg.taus(itau)));
    end
    
    wSMIout.prc_used_raw = percent; % percent of data used
    wSMIout.dur_used_raw = durgood; % seconds of data used
    wSMIout.prc_win_used = gcount/(gcount+bcount)*100;
    wSMIout.win_len_used = gcount;
    wSMIout.win_len      = nwin;
    wSMIout.fsample      = cfg.sf;
    wSMIout.surrogate    = surrogate; % Boolean, is this surrogate data?
    
    outstr = data.fstr; % cutoff '.mat' suffix
    
    if sleep && surrogate
        outname = sprintf('%s%s%s%s',DIRRESULT,outstr,'_sleep_wSMI_SURROGATE','.mat');
    elseif sleep && ~surrogate
        outname = sprintf('%s%s%s%s',DIRRESULT,outstr,'_sleep_wSMI','.mat');
    elseif ~sleep && surrogate
        outname = sprintf('%s%s%s%s',DIRRESULT,outstr,'_wake_wSMI_SURROGATE','.mat');
    elseif ~sleep && ~surrogate
        outname = sprintf('%s%s%s%s',DIRRESULT,outstr,'_wake_wSMI','.mat');
    end
    
    save(outname,'wSMIout','-v7.3');
    clear data wSMIout
    % end
end


end
