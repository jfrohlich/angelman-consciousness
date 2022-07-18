% run function once for sleep data, once for awake data
clearvars
dbstop if error
OS = computer; % detect opperating system

switch OS
    case 'GLNXA64'
        pth = './2021_analysis/postprocessed/AS/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
    case 'MACI64'
        pth = './2021_analysis/postprocessed/AS/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
    case 'PCWIN64'
        pth = './2021_analysis/postprocessed/AS/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
end

% these are the files we already know we want to use
% load('select_files','select_files')

sleep = true;
compute_FT(sleep,files,pth)

sleep = false;
compute_FT(sleep,files,pth)


%%

function[] = compute_FT(sleep,files,pth)

switch sleep
    case true
        DIRRESULT = './2021_analysis/freq_out/AS/nan=0.2/sleep/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end
    case false
        DIRRESULT = './2021_analysis/freq_out/AS/nan=0.2/wake/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end
end


start_ndx = 1;

if start_ndx ~= 1, warning('start_ndx is set to %i',start_ndx), end

for ifile = start_ndx:length(files)
    
    %    if logical(sum(ismember(select_files,files(ifile).name(1:end-4))))
    load(sprintf('%s%s',pth,files(ifile).name),'data')
    fprintf('Now loading %s \n',files(ifile).name)
    
    % Skip if no useable sleep
    if ~isfield(data.cfg.dattype,'SLEEP_MONTI') || isempty(data.cfg.dattype.SLEEP_MONTI)
        fprintf('\nSkipping this one, no useable sleep ... \n')
        continue
    end
    
    %Get datamat var for freq transform
    datamat = data.trial{1};
    
    
    % Do some santiy checks
    assert(size(datamat,1) == 19,'channel(s) missing!')
    % allow for small number of NaNs at beginning and end of file
    assert(sum(isnan(datamat(1:end)))/numel(datamat) < 0.01, 'Many NaNs in data!')
    
    if isfield(data,'clinON'), error('This is an overnight EEG'), end
    
    
    switch sleep
        case true
            % convert awake data to NaNs
            % Data to retain for analysis (NOT converted to NaNs!)
            
            keep = [data.cfg.dattype.SLEEP_ALL]; % only segments that we're sure are sleep
            
            datamat2 = ones(size(datamat)); % ones indicate data that will be switched to NaNs
            
            for islp = 1:size(keep,1)
                idpt = keep(islp,1) - data.sampleinfo(1,1) + 1;
                datamat2(:,idpt:idpt+diff(keep(islp,:))) = 0;
            end
            
            wkedx = datamat2 == 1;
            
            datamat(wkedx) = NaN; % replace awake data with NaNs
            badart = [];
        case false
            badart = data.cfg.dattype.SLEEP_ALL; % exclude anything that could possiblly be sleep
    end
    
    % Set bad data to NaNs for purpose of taking freq transform
    badart = [badart; data.cfg.dattype.bad; data.cfg.dattype.flash; data.cfg.dattype.drowsy];
    
    % add muscle to artifacts if it is a field
    if isfield(data.cfg.dattype,'muscle')
        badart = [badart; data.cfg.dattype.muscle];
    end
    
    badart(badart==0) = 1; % do this so no 0 indicies that give bugs
    
    for iart = 1:size(badart,1)
        idpt = badart(iart,1) - data.sampleinfo(1,1) + 1;
        datamat(:,idpt:idpt+diff(badart(iart,:))) = nan;
    end
    
    % determine amount of good data
    c = sum(~isnan(datamat(1,:)));
    durgood = c/data.fsample; % seconds of data used
    if durgood < 30
        fprintf('Skipping dataset, less than 30 seconds of good data ...\n')
        continue
    else
        display(strcat('Data has this many seconds of good data: ',num2str(round(durgood))))
    end
    
    percent = c/length(datamat); % proportion of data used
    
    
    %% do frequency transform
    
    cfg.foi_start = 0.5;
    cfg.foi_end   = 32;
    cfg.fsample = data.fsample;

    % Take out flanking NaNs to speed up frquency transform--this changes
    % the output for pow_full but not pow
    frstchn = ~isnan(datamat(1,:));
    ndxA = find(frstchn,1,'first'); % find the first non-NaN element
    ndxB = find(frstchn,1,'last'); % find the last non-NaN element
    datamat = datamat(:,ndxA:ndxB);

    % Allow NaNs for consistancy with other analysis
    cfg.allow_fraction_nan=0.2; % allow 20% NaNs in window
    [pow,pow_full,pow_median,pow_var,n,unit,foi_target,foi_delta,foi] = ro_freq_wavelet_JF(datamat,cfg);
    
    outstr = data.fstr;
    
    if sleep
        outname = sprintf('%s%s%s%s',DIRRESULT,outstr,'_sleep_freq','.mat');
    else
        outname = sprintf('%s%s%s%s',DIRRESULT,outstr,'_wake_freq','.mat');
    end
    
    save(outname,'pow','pow_full','pow_median','pow_var','n','unit',...
        'foi_target','foi_delta','foi','durgood','percent','-v7.3');
    clear data pow pow_full pow_median pow_var n unit foi_target ...
        foi_delta foi durgood percent
end

end

