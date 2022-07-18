% run function once for sleep data, once for awake data
clearvars
rng(45783) % seed added 02.26.19
OS = computer; % detect opperating system
dbstop if error

switch OS
    case 'GLNXA64'
        pth = './2021_analysis/postprocessed/TD/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
    case 'MACI64'
        pth = './2021_analysis/postprocessed/TD/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
    case 'PCWIN64'
        pth = './2021_analysis/postprocessed/TD/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
end

sleep = true; % asleep or awake?
surrogate = false; % scramble FFT phases?
compute_criticality(sleep,surrogate,files,pth)

sleep = false; % asleep or awake?
surrogate = false; % scramble FFT phases?
compute_criticality(sleep,surrogate,files,pth)

sleep = true; % asleep or awake?
surrogate = true; % scramble FFT phases?
compute_criticality(sleep,surrogate,files,pth)

sleep = false; % asleep or awake?
surrogate = true; % scramble FFT phases?
compute_criticality(sleep,surrogate,files,pth)


function[] = compute_criticality(sleep,surrogate,files,pth)
switch sleep
    case true
        DIRRESULT = './2021_analysis/TD_output/CRIT/sleep/';
        ntrl = 2;
        if ~exist(DIRRESULT), mkdir(DIRRESULT), end
    case false
        DIRRESULT = './2021_analysis/TD_output/CRIT/wake/';
        ntrl = 1;
        if ~exist(DIRRESULT), mkdir(DIRRESULT), end
end

tic

startndx = 1;

if startndx ~= 1, warning('Starting index set to %i',startndx), end

for ifile = startndx:length(files)
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
    
    %%% DON'T DOWNSAMPLE, THIS WILL ALREADY BE DONE INSIDE CHAOS_TEST %%%
    
    % re-average-ref datamat (after interpolation and ICA)
    datamat = datamat - nanmean(datamat,1);
    
    % Take out flanking NaNs for filter
    frstchn = ~isnan(datamat(1,:));
    ndxA = find(frstchn,1,'first'); % find the first non-NaN element
    ndxB = find(frstchn,1,'last'); % find the last non-NaN element
    datamat = datamat(:,ndxA:ndxB);
    
    % Do some santiy checks
    try
        assert(size(datamat,1) == 19,'Wrong number of channels')
    catch
        if size(datamat,1) > 19
            ch = {'C3','C4','O1','O2','Cz','F3','F4','F7','F8','Fz','Fp1','Fp2',...
                'P3','P4','Pz','T3','T4','T5','T6'};
            chidx = contains(data.label,ch);
            datamat = datamat(chidx,:); % discard aux channels
        else
            error('Channels missing')
        end
    end
    % allow for small number of NaNs at beginning and end of file
    assert(sum(isnan(datamat(1:end)))/numel(datamat) < 0.01, 'Many NaNs in data!')
    %nmin = 30; % first n minutes of data
       
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
    
    %% Compute the criticality
    
    win_len = 10;
    % windowing parameters
    n_win = win_len*data.fsample;
    
    cfg.critalpha = 0.05; % Good value found from preliminary analysis
    cfg.method = 'minmax';
    
    window_shift=1; % 1.0 --> no overlap; 0.5 --> 50 percent
    n_shift = round(n_win*window_shift);
    L = length(1:n_shift:size(datamat,2)-n_win+1);
    Cr = nan(size(datamat,1),floor(durgood/win_len)); % preallocation
    valid = nan(size(datamat,1),floor(durgood/win_len)); % preallocation
    K = nan(size(datamat,1),floor(durgood/win_len)); % preallocation
    FLP = nan(size(datamat,1),floor(durgood/win_len)); % preallocation
    cnt = 0;
    tic
    for isection = 1:n_shift:size(datamat,2)-n_win+1
        for ich = 1:size(datamat,1)
            section = double(datamat(ich,isection:isection+n_win-1));
            if all(~isnan(section)) && ~all(section==0) % if there are no nans and the data are not just zeros
                if ich == 1, cnt = cnt + 1; end % increment counter if it's first channel
                if surrogate
                    surdat = surrogate_data(section,data.fsample);
                    % overwrite data section with the surrogate data
                    section = surdat;
                end
                                
                %%% Do filtering with FOOOF %%%
                n_butt_EEG = 5; % 5th order butterworth filter
                % Get cutoff frequency using FOOOF
                try
                    f_lp = select_low_pass_freq(section,data.fsample);
                    assert(~isnan(f_lp))
                catch
                    %fprintf('Failed to get frequency cutoff from FOOOF, skipping this one ...\n')
                    continue
                end
                fnyq = data.fsample/2;
                [Blp, Alp] = butter(n_butt_EEG,f_lp/fnyq,'low');
                section = filtfilt(Blp, Alp, double(section)')';
                section = single(section);
                
                n_valid = sum(~isnan(section));
                k = chaos_test(section,cfg.method);
                Cr(ich,cnt) = criticality(k,cfg.critalpha); 
                valid(ich,cnt) = n_valid;
                K(ich,cnt) = k;
                FLP(ich,cnt) = f_lp; 
                if ich == 1
                    fprintf('%2.1f%% complete (%i min %i sec elapsed)\n',...
                        cnt/(durgood/win_len)*100,floor(toc/60),floor(mod(toc,60)))
                end
            else
                %fprintf('Skipping NaNs ...\n')
                break % don't bother with other channels for this section
            end
        end
    end
    Crit.Cr = Cr;
    Crit.K = K;
    Crit.valid = valid;
    Crit.FLP = FLP;
    Crit.cfg = cfg;
    Crit.surrogate = surrogate; % Boolean, is this a surrogate?
    Crit.prc_used = percent;    % percent of data used
    Crit.dur_used = durgood;    % seconds of data used

    try
        outstr = data.fstr; % cutoff '.mat' suffix
    catch
        outstr = files(ifile).name(1:end-4);
    end

    if sleep && surrogate
        outname = sprintf('%s%s%s%s',DIRRESULT,outstr,'_sleep_CRIT_SURROGATE','.mat');
    elseif sleep && ~surrogate
        outname = sprintf('%s%s%s%s',DIRRESULT,outstr,'_sleep_CRIT','.mat');
    elseif ~sleep && surrogate
        outname = sprintf('%s%s%s%s',DIRRESULT,outstr,'_wake_CRIT_SURROGATE','.mat');
    elseif ~sleep && ~surrogate
        outname = sprintf('%s%s%s%s',DIRRESULT,outstr,'_wake_CRIT','.mat');
    end

    save(outname,'Crit','-v7.3');
    clear data datamat Crit
end

end
