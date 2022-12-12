% Joel Frohlich
% University of California, Los Angeles (UCLA)
% Monti Lab, Psychology Department
%
% University of Tuebingen, Germany
% Institue for Neuromodulation and Neurotechnology 
%
% Last update: 12 Dec, 2022 (cleaned up code and comments)

%     The code compute_MSE() takes in a sleep parameter that specifies
% whether to use data from when the subjects are asleep or awake, and a
% surrogate parameter that specifies whether to scramble the FFT phases or
% not. The function creates a directory for the output based on the values
% of the sleep and surrogate parameters. The function loops through the
% .mat files in the specified directory, loading each file and extracting
% the data. The function performs some sanity checks on the data, converts
% the data to NaNs if necessary, and then applies  the mMSE algorithm to
% the data. The results are saved in the output directory.


% run function once for sleep data, once for awake data
clearvars
rng(45783) % seed added 02.26.19
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

gpu = true;

sleep = true; % asleep or awake?
surrogate = false; % scramble FFT phases?
compute_MSE(sleep,surrogate,files,pth,gpu)

sleep = false; % asleep or awake?
surrogate = false; % scramble FFT phases?
compute_MSE(sleep,surrogate,files,pth,gpu)

sleep = true; % asleep or awake?
surrogate = true; % scramble FFT phases?
compute_MSE(sleep,surrogate,files,pth,gpu)

sleep = false; % asleep or awake?
surrogate = true; % scramble FFT phases?
compute_MSE(sleep,surrogate,files,pth,gpu)


%%

function[] = compute_MSE(sleep,surrogate,files,pth,gpu)

switch sleep
    case true
         DIRRESULT = './2021_analysis/AS_output/Xie/r=0.15/dynr/sleep/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end
    case false
        DIRRESULT = './2021_analysis/AS_output/Xie/r=0.15/dynrwake/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end
end

switch gpu
    case true
        fprintf('GPU option is set to TRUE\n')
    case false
        fprintf('Note: GPU option is set to false. Switch to true to run faster\n.')
end

switch sleep
    case true
        start_ndx = 1;
    case false
        start_ndx = 1;
end

if start_ndx ~= 1
    fprintf('Warning, start_ndx is set at %i',start_ndx)
end

for ifile = start_ndx:length(files)
%    if logical(sum(ismember(select_files,sprintf('%s%s',files(ifile).name(1:end-28),'_wake_MSE.mat')))) ...
%            || logical(sum(ismember(select_files,sprintf('%s%s',files(ifile).name(1:end-4),'_wake_MSE.mat'))))
        load(sprintf('%s%s',pth,files(ifile).name),'data')
        fprintf('Now loading %s \n',files(ifile).name)
        
        % check if this file had useable sleep
        
        if ~isfield(data.cfg.dattype,'SLEEP_MONTI') || isempty(data.cfg.dattype.SLEEP_MONTI)
            fprintf('\nSkipping this one, no useable sleep ... \n')
            continue
        end
        
        %Get datamat var for freq transform
        datamat = data.trial{1};
        
        assert(sum(isnan(mean(datamat)))~=length(mean(datamat)),'All data are NaNs!')
        
        % Do some santiy checks
        assert(size(datamat,1) == 19,'channel(s) missing!')
        % allow for small number of NaNs at beginning and end of file
        assert(sum(isnan(datamat(1:end)))/numel(datamat) < 0.01, 'Many NaNs in data!')
       
        switch sleep
            case true
                %% convert awake data to NaNs
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
                badart = [data.cfg.dattype.SLEEP_ALL]; % exclude anything that could possiblly be sleep
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
                        surdat(ich,:) = surrogate_data(section(ich,:),srate);
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
%    end
end

end
