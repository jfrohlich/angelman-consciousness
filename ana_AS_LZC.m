clearvars
rng(45783) % seed added 02.26.19
OS = computer; % detect opperating system
cd 'C:\Users\joelf\OneDrive\Documents\Monti\Universal_scripts\EntRate-master'
mex COPTIMFLAGS="-O3" C:\Users\joelf\OneDrive\Documents\Monti\Universal_scripts\EntRate-master\LZ76.c % Compile LZ76 code
cd 'C:\Users\joelf\OneDrive\Documents\Monti\angelman_nih\Monti'

%%
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


%match = false; % power matched sections?

% sleep = true; % asleep or awake?
% surrogate = false; % scramble FFT phases?
% compute_LZC(sleep,surrogate,files,pth)

sleep = false; % asleep or awake?
surrogate = false; % scramble FFT phases?
compute_LZC(sleep,surrogate,files,pth)

% sleep = true; % asleep or awake?
% surrogate = true; % scramble FFT phases?
% compute_LZC(sleep,surrogate,files,pth)
% 
% sleep = false; % asleep or awake?
% surrogate = true; % scramble FFT phases?
% compute_LZC(sleep,surrogate,files,pth)

%%

function[] = compute_LZC(sleep,surrogate,files,pth)

switch sleep
    case true
        DIRRESULT = './2021_analysis/AS_output/LZCv/sleep/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end
    case false
        DIRRESULT = './2021_analysis/AS_output/LZCv/wake/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end
end



%% Parameters for signal smoothing (if needed, i.e., unused for LZCv)
% frequency start, stop, intervals
cfg.foi_start=0.5;
cfg.foi_end=30;
cfg.oct_delta=0.25;
srate = 200; % downsample to this rate to control for sampling rate
% desired center frequencies
tmp = 2.^[log2(cfg.foi_start):cfg.oct_delta:log2(cfg.foi_end)];
% number of sample to smooth
smooth = srate./tmp;
% frequencies that are possible, given the sampling rate
foi = srate./myround(smooth,'odd');

% upper and lower frequency boundaries for gLZC
lo = foi.*0.8; % for thresholding
hi = foi.*1.2; % for attenuating high frequencies
smooth_lo = myround(srate./lo,'odd');
smooth_hi = myround(srate./hi,'odd');
lo = srate./smooth_lo;
hi = srate./smooth_hi;
scale = 1:length(foi);
ftab = table(scale',foi',hi'-lo',lo',hi',smooth_lo',smooth_hi');
writetable(ftab,'gLZV_scales.csv');

%%

start_ndx = 1;

if start_ndx ~= 1, warning('start_ndx set to %i',start_ndx), end

for ifile = start_ndx:length(files)
    %if logical(sum(ismember(select_files,files(ifile).name(1:end-4))))
    load(sprintf('%s%s',pth,files(ifile).name),'data')
    fprintf('Now loading %s \n',files(ifile).name)
    scnt = 0; % section counter
    
    % check that it isn't a clinical overnight EEG
    if isfield(data,'clinON'), error('Clinical overnight file'), end
    
    % check if this file had useable sleep
    
    if ~isfield(data.cfg.dattype,'SLEEP_ALL') || isempty(data.cfg.dattype.SLEEP_ALL)
        fprintf('\nSkipping this one, no useable sleep ... \n')
        continue
    end
    
    %Get datamat 
    datamat = data.trial{1};
    
    % check that data are averaged referenced
    assert(all(nanmean(datamat,1) < 2e-3),'Data are not averaged referenced');
    
    % Do some santiy checks
    assert(size(datamat,1) == 19,'channel(s) missing!')
    
    % allow for small number of NaNs at beginning and end of file
    assert(sum(isnan(datamat(1:end)))/numel(datamat) < 0.01, 'Many NaNs in data!')

    switch sleep
        case true
            %% convert awake data to NaNs
           % Data to retain for analysis (NOT converted to NaNs!)
            keep = data.cfg.dattype.SLEEP_ALL; % only segments that we're sure are sleep
            datamat2 = ones(size(datamat)); % ones indicate data that will be switched to NaNs
            
            for islp = 1:size(keep,1)
                idpt = keep(islp,1) - data.sampleinfo(1,1) + 1;
                datamat2(:,idpt:idpt+diff(keep(islp,:))) = 0;
            end
            
            wkedx = datamat2 == 1;
            
            datamat(wkedx) = NaN; % replace awake data with NaNs
            badart = [];
        case false
            % SLEEP_ALL + Drowsy = sleep (old field), we exclude both
            badart = data.cfg.dattype.SLEEP_ALL; 
    end
    
    % Set bad data to NaNs for purpose of taking freq transform
    badart = [badart; data.cfg.dattype.bad; data.cfg.dattype.flash; data.cfg.dattype.drowsy ];
    
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
    
    %% partion into n point windows and compute LZC
    
    % make sure window length is expressed in terms of the *NEW* sampling
    % rate, otherwise we have wrong size for sections
    nswin = 60; % recommended to set same window size as mMSE
    win_len = nswin*srate;
    
    cfg = [];
    cfg.verbose=1;
    cfg.native_srate = data.fsample;
    cfg.new_srate = srate;
    cfg.verbose = false;
    cfg.foi = foi;
    cfg.smooth_lo = smooth_lo;
    cfg.smooth_hi = smooth_hi;
    cfg.lo = lo;
    cfg.hi = hi;
    

    % Take out flanking NaNs to speed things up
    frstchn = ~isnan(datamat(1,:));
    ndxA = find(frstchn,1,'first'); % find the first non-NaN element
    ndxB = find(frstchn,1,'last'); % find the last non-NaN element
    datamat = datamat(:,ndxA:ndxB);
    datamat = datamat(:,1:data.fsample/srate:end); % downsample to 200 Hz, don't use filters
    
    vLZC = [];
    gLZC = [];
    CTW = [];
    CSER = [];
    n_valid_vLZC = [];
    n_valid_gLZC = [];
    n_valid_CTW = [];
    n_valid_CSER = [];
    
    sections = cell(1,1); % hold sections for later reanalysis with the LZC decomp from Pedro
    
    % windowing parameters
    n_win = win_len;
    % 1.0 --> no overlap; 0.5 --> 50 percent
    window_shift=1; % otherwise do no overlap
    
    n_shift = round(n_win*window_shift);
    
    tic;
    for isection = 1:n_shift:size(datamat,2)-n_win+1
        section = double(datamat(:,isection:isection+n_win-1));        
        
        if surrogate
            surdat = nan(size(section,1),sum(~isnan(mean(section)))); % allocation
            for ich = 1:size(section,1)
                surdat(ich,:) = surrogate_data(section(ich,:));
            end
            % overwrite data section with the surrogate data
            section = surdat;
        end
        [cfg_out] = ro_LZCv(section,cfg);
        
        if ~surrogate % only do this for the real data
            scnt = scnt + 1;
            sections{scnt} = section;
        end
        
        vLZC = cat(2,cfg_out.vLZC,vLZC);
        gLZC = cat(3,cfg_out.gLZC,gLZC);
        CTW = cat(2,cfg_out.CTW,CTW);
        CSER = cat(2,cfg_out.CSER,CSER);

        n_valid_vLZC = cat(2,cfg_out.n_valid_vLZC,n_valid_vLZC);
        n_valid_gLZC = cat(3,cfg_out.n_valid_gLZC,n_valid_gLZC);
        n_valid_CTW = cat(2,cfg_out.n_valid_CTW,n_valid_CTW);
        n_valid_CSER = cat(3,cfg_out.n_valid_CSER,n_valid_CSER);
        fprintf('     %2.3f%% complete (%i min %i sec elapsed) \n',...
               (isection+n_win-1)/size(datamat,2)*100,floor(toc/60),floor(mod(toc,60)))
    end

    LZCout.cfg = cfg_out; % includes info on center frequencies and smoothing kernels
    LZCout.sections = sections; % data sections used for LZC, so we can use these later for the LZC decomp from Pedro
    LZCout.vLZC = vLZC;
    LZCout.gLZC = gLZC;
    LZCout.CTW = CTW;
    LZCout.CSER = CSER;
    LZCout.n_valid_vLZC = n_valid_vLZC;
    LZCout.n_valid_gLZC = n_valid_gLZC;
    LZCout.n_valid_CTW = n_valid_CTW;
    LZCout.n_valid_CSER = n_valid_CSER;
    LZCout.window_shift = window_shift; % window overlap
    
    LZCout.prc_used = percent; % percent of data used
    LZCout.dur_used = durgood; % seconds of data used
    LZCout.surrogate = surrogate; % Boolean, is this surrogate data?
    
    outstr = data.fstr; % cutoff '.mat' suffix
    
    if sleep && surrogate
        outname = sprintf('%s%s_sleep_LZCv_SURROGATE_wlen=%is',DIRRESULT,outstr,nswin);
    elseif sleep && ~surrogate
        outname = sprintf('%s%s_sleep_LZCv_wlen=%is',DIRRESULT,outstr,nswin);
    elseif ~sleep && surrogate
        outname = sprintf('%s%s_wake_LZCv_SURROGATE_wlen=%is',DIRRESULT,outstr,nswin);
    elseif ~sleep && ~surrogate
        outname = sprintf('%s%s_wake_LZCv_wlen=%is',DIRRESULT,outstr,nswin);
    end
    
    save(outname,'LZCout','-v7.3');
    clear data LZCout
    % end
end

end

function[out] = myround(num,option)

% This is a wrapper version of the native MATLAB round function that allows
% user to specify whether we should round towards nearest OOD integer
% (option = 'odd') or nearest EVEN integer (option = 'even'). Setting option to
% 'closest' just uses the native MATLAB round() function. Takes row vectors
% or scalars as input.

if nargin < 2
    option = 'closest';
end

assert(size(num,1) == 1);

out = nan(size(num)); % allocation

switch option
    case 'even'
        for i = 1:length(num)
            down = floor(num(i));
            if mod(down,2) == 0 % if rounding down gives up an even number
                out(i) = down;
            else
                out(i) = down+1; % round up instead!
            end
        end
    case 'odd'
        for i = 1:length(num)
            down = floor(num(i));
            if mod(down,2) == 1 % if rounding down gives up an odd number
                out(i) = down;
            else
                out(i) = down+1; % round up instead!
            end
        end
    case 'either'
        for i = 1:length(num)
            out(i) = round(num(i));
        end
    otherwise
        error('Option not detected')
end

end
