% run function once for sleep data, once for awake data
clearvars
rng(45783) % seed added 02.26.19
OS = computer; % detect opperating system

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

gpu = true;

sleep = true; % asleep or awake?
surrogate = false; % scramble FFT phases?
compute_MSE(sleep,surrogate,files,pth,gpu)

sleep = false; % asleep or awake?
surrogate = false; % scramble FFT phases?
compute_MSE(sleep,surrogate,files,pth,gpu)

%%

function[] = compute_MSE(sleep,surrogate,files,pth,gpu)

switch sleep
    case true
        DIRRESULT = './2021_analysis/Dup15q_output/Xie/r=0.15/dynr/sleep/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end
        ntrl = 2;
    case false
        DIRRESULT = './2021_analysis/Dup15q_output/Xie/r=0.15/dynr/wake/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end
        ntrl = 1;
end

switch gpu
    case true
        fprintf('GPU option is set to TRUE\n')
    case false
        fprintf('Note: GPU option is set to false. Switch to true to run faster\n.')
end

start_ndx = 1;

if start_ndx ~= 1
    fprintf('Warning, start_ndx is set at %i',start_ndx)
end

for ifile = start_ndx:length(files)
    load(sprintf('%s%s',pth,files(ifile).name),'data')
    fprintf('Now loading %s \n',files(ifile).name)
    
    % check if this file had useable sleep
    % The second trial is always sleep, so this might not be necessary
    if ~isfield(data.cfg.dattype,'SLEEP_MONTI') || isempty(data.cfg.dattype.SLEEP_MONTI)
        fprintf('Skipping this one, no useable sleep ... \n')
        continue
    end
    
    datamat = data.trial{ntrl}; % sleep is stored in second trial
    
    assert(sum(isnan(mean(datamat)))~=length(mean(datamat)),'All data are NaNs!')
    
    % Do some santiy checks
    assert(size(datamat,1) == 19,'channel(s) missing!')
    
    % allow for small number of NaNs at beginning and end of file
    assert(sum(isnan(datamat(1:end)))/numel(datamat) < 0.01, 'Many NaNs in data!')
    
    % Set bad data to NaNs for purpose of taking freq transform
    badart = [data.cfg.dattype.bad; data.cfg.dattype.flash];
    
    % add muscle to artifacts if it is a field
    if isfield(data.cfg.dattype,'muscle')
        badart = [badart; data.cfg.dattype.muscle];
    end
    
    for iart = 1:size(badart,1)
        idpt = badart(iart,1) - data.sampleinfo(ntrl,1) + 1;
        if idpt > 0 && idpt+diff(badart(iart,:)) <= size(datamat,2)
            datamat(:,idpt:idpt+diff(badart(iart,:))) = nan;
        end
    end
    
    badart(badart==0) = 1; % do this so no 0 indicies that give bugs
    
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
    
    percent = c/length(datamat);
    durgood = c/data.fsample;
    
    %% partion into 30-sec windows with x percent overlap and compute MSE
    
    % longer windows are better (see Grandy et al 2016)
    
    % Do MSE
    srate = 200; % downsample to this rate to control for sampling rate
    cfg = [];
    cfg.verbose=true;
    cfg.tolerance = 0.15; % fraction of SD to use as state space radius
    cfg.m = 2;
    cfg.native_srate = data.fsample;
    %cfg.new_srate = cfg.native_srate;
    cfg.scales=1:20; % 20 time scales to include delta frequencies
    cfg.verbose = false;
    cfg.type='Xie'; % less sensitive to tolerance/r
    cfg.gpu = gpu;
    cfg.dynr = true;
    
    mse_3d = [];
    scale = [];
    n_valid = [];
    r = [];
    
    win_len = 30;
    datamat = datamat(:,1:data.fsample/srate:end); % downsample to 200 Hz, don't use filters
    cfg.new_srate = srate; % NOTE: this was missing before, older files won't show the new sampling rate in cfg
    % Take out flanking NaNs to speed things up
    frstchn = ~isnan(datamat(1,:));
    ndxA = find(frstchn,1,'first'); % find the first non-NaN element
    ndxB = find(frstchn,1,'last'); % find the last non-NaN element
    datamat = datamat(:,ndxA:ndxB);
    % windowing parameters
    n_win = win_len*srate;
    % 1.0 --> no overlap; 0.5 --> 50 percent
    window_shift=1; % do no overlap
    
    n_shift = round(n_win*window_shift);
    for isection = 1:n_shift:size(datamat,2)-n_win+1
        section = double(datamat(:,isection:isection+n_win-1));
        if surrogate
            if sum(isnan(mean(section))) == length(section) % if the whole thing is NaNs
                surdat = nan(size(section));
            else
                surdat = nan(size(section,1),sum(~isnan(mean(section)))); % allocation
                for ich = 1:size(section,1)
                    surdat(ich,:) = surrogate_data(section(ich,:));
                end
            end
            % overwrite data section with the surrogate data
            section = surdat;
        end
        [mse,scale_out,n_valid_out,r_out] = ro_mse(section,cfg);
        mse_3d = cat(3,mse_3d,mse);
        scale = [scale; scale_out];
        n_valid = [n_valid n_valid_out];
        r = cat(3,r,r_out);
        fprintf('%2.2f percent complete \n',isection/(size(datamat,2)-n_win+1)*100)
    end
    if gpu
        MSEout.mse = gather(mse_3d);
    else
        MSEout.mse = mse_3d;
    end
    MSEout.scale = scale;
    MSEout.n_valid = n_valid;
    MSEout.r = r;
    MSEout.win_len = win_len;
    MSEout.window_shift = window_shift;
    MSEout.cfg = cfg;
    MSEout.surrogate = surrogate; % Boolean, is this a surrogate?
    MSEout.prc_used = percent;    % percent of data used
    MSEout.dur_used = durgood;    % seconds of data used
    
    outstr = data.fstr; % cutoff '.mat' suffix
    
    if sleep && surrogate
        outname = sprintf('%s%s%s%s',DIRRESULT,outstr,'_sleep_MSE_SURROGATE','.mat');
    elseif sleep && ~surrogate
        outname = sprintf('%s%s%s%s',DIRRESULT,outstr,'_sleep_MSE','.mat');
    elseif ~sleep && surrogate
        outname = sprintf('%s%s%s%s',DIRRESULT,outstr,'_wake_MSE_SURROGATE','.mat');
    elseif ~sleep && ~surrogate
        outname = sprintf('%s%s%s%s',DIRRESULT,outstr,'_wake_MSE','.mat');
    end
    
    save(outname,'MSEout','-v7.3');
    clear data MSEout
end

end
