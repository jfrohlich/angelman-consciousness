% run function once for sleep data, once for awake data
clearvars
dbstop if error
OS = computer; % detect opperating system

switch OS
    case 'GLNXA64'
        pth = './2021_analysis/postproc_2019/Dup15q/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
    case 'MACI64'
        pth = './2021_analysis/postprocessed/Dup15q/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
    case 'PCWIN64'
        pth = './2021_analysis/postprocessed/Dup15q/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
end

sleep = true;
compute_FT(sleep,files,pth)
sleep = false;
compute_FT(sleep,files,pth)


%%

function[] = compute_FT(sleep,files,pth)

switch sleep
    case true
        DIRRESULT = './2021_analysis/freq_out/Dup15q/nan=0.2/sleep/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end
        ntrl = 2;
    case false
        DIRRESULT = './2021_analysis/freq_out/Dup15q/nan=0.2/wake/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end
        ntrl = 1;
end

for ifile = 1:length(files)
    load(sprintf('%s%s',pth,files(ifile).name),'data')
    fprintf('Now loading %s \n',files(ifile).name)
    
    % Skip if no useable sleep
    if ~isfield(data.cfg.dattype,'SLEEP_MONTI') || isempty(data.cfg.dattype.SLEEP_MONTI)
        fprintf('\nSkipping this one, no useable sleep ... \n')
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
    cfg.allow_fraction_nan=0.2; % allow n% NaNs in window
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

