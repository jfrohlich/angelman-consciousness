function [] = ana_preproc_as_Dup15q_wsleep_butterHP_firLP()


% Load and pre-process the Dup15q data controls with sleep
% To do: change script so that everything saves in one single trial
% (individual .EEG files should be continuous)

dbstop if error

% Channel labels
%label_alt = {'C3','C4','O1','O2','Cz','F3','F4','F7','F8','Fz','FP1','FP2','P3','P4','Pz','T7','T8','P7','P8'};
label_load  = {'C3','C4','O1','O2','Cz','F3','F4','F7','F8','Fz','FP1','FP2','P3','P4','Pz','T3','T4','T5','T6'};
label_load_alt = {'EEG C3','EEG C4','EEG O1','EEG O2','EEG Cz','EEG F3','EEG F4',...
    'EEG F7','EEG F8','EEG Fz','EEG FP1','EEG FP2','EEG P3','EEG P4','EEG Pz','EEG T3','EEG T4','EEG T5','EEG T6'};
label_physio = {'EKG','EOG'};
label_physio_alt = {'ECG','ECG','LOC','ROC'};
load angelman_lay.mat lay elec

% Filter settings
% FIR filter parameter (for data)
f_hp = 0.4; % Hz
f_lp = 45; % Hz
% Order for Butterworth filter
n_butt = 4;
n_butt_EEG = 5;
% EOG filter parameter
f_hp_eog = 1; f_lp_eog = 15;
% ECG filter parameter
f_hp_ecg = 2; f_lp_ecg = 50;

% Filter defintion, use as
% for idef=1:length(filter_def), eval(filter_def{idef}), end
filter_def{1} = 'fnyq = fsample/2;'; % nyquest frequency
filter_def{2} = '[B_eog, A_eog] = butter(n_butt,[f_hp_eog,f_lp_eog]/fnyq);';            % EOG bandpass filter coefficients
filter_def{3} = '[B_ecg, A_ecg] = butter(n_butt,[f_hp_ecg,f_lp_ecg]/fnyq);';            % ECG bandpass filter coeffcients
filter_def{4} = 'for icnt=1:floor(fnyq/60); [B_line{icnt}, A_line{icnt}] = butter(n_butt,(icnt*60+[-1,1])/fnyq,''stop''); end'; % US line noise notch-filter coefficents
filter_def{5} = 'for icnt=1:floor(fnyq/30); [B_line2{icnt}, A_line2{icnt}] = butter(n_butt,(icnt*30+[-1,1])/fnyq,''stop''); end'; % US SUBHARMONIC line noise notch-filter coefficents
filter_def{6} = '[Bhp, Ahp] = butter(n_butt_EEG,[f_hp]/fnyq,''high'');';
filter_def{7} = '[Blp, Alp] = butter(n_butt_EEG,[f_lp]/fnyq,''low'');';

wakestart = '17:00:00'; % clock time to begin extracting wake EEG
wakeend = '17:30:00'; % clock time to end wake EEG

% wakestart = '15:30:00'; % clock time to begin extracting wake EEG
% wakeend = '16:00:00'; % clock time to end wake EEG

%% Look for files and import them

% set directory to save results
DIRRESULT  = './Dup15q_controls/Imported/Dup15q_butterHP_firLP/';  if ~exist(DIRRESULT), mkdir(DIRRESULT), end

% look in folder
subs_strc = dir('.\Sleep_EEG_Dup15q\EDFs for Joel\');

T_sleep = readtable('.\Dup15q_controls\Dup15q_sleep_times.csv');


% get just the subdirectory names;
subs = cell(1,1);
for i = 1:length(subs_strc)
    subs{i} = subs_strc(i).name;
end
subs = setdiff(subs,{'.','..'}); % prune junk

% create empty cell arrays (for table at the end)
lostchans = cell(length(subs),1); % missing channels
ntrials = cell(length(subs),1); % number of trials
whichtrials = cell(length(subs),1); % which trials were used
datalen = cell(length(subs),1); % total data length

start_ndx = 11;
if start_ndx ~= 1, fprintf('\nWarning, starting index = %i\n',start_ndx), end

for isub = start_ndx:length(subs)
    ID = subs{isub}; % subject ID
    destination = sprintf('./Sleep_EEG_Dup15q/EDFs for Joel/%s',ID);
    fprintf('\nNow importing data from subject %s\n',ID)
    
    tableidx = find(T_sleep.ID == str2num(ID)); % find matching entries in table
    eegdate = unique(T_sleep.EEGDate(tableidx));
    assert(length(eegdate)<=1,'More than one date for EEG visit detected')

    for ientry = 1:length(tableidx)
        condition = 'neither'; % initalize condition variable 
        
        % file start/stop times
        a = T_sleep.StartTime(tableidx(ientry));
        b = T_sleep.EndTime(tableidx(ientry));
        
        if a < duration('11:55:00') % if segment starts in the morning 
            % use 11:55 am as the earliest time of day that an assessment could
            % possibly start
            % add a day (because it's the next morning)
            A = datetime(year(eegdate),month(eegdate),day(eegdate)+1,...
                floor(hours(a)),floor(minutes(a)-floor(hours(a))*60),...
                floor(seconds(a)-floor(minutes(a))*60));
        else
            A = datetime(year(eegdate),month(eegdate),day(eegdate),...
                floor(hours(a)),floor(minutes(a)-floor(hours(a))*60),...
                floor(seconds(a)-floor(minutes(a))*60));
        end
        
        if b < duration('12:00:00') % if segment ends in the morning before 10 am
            % add a day (because it's the next morning)
            B = datetime(year(eegdate),month(eegdate),day(eegdate)+1,...
                floor(hours(b)),floor(minutes(b)-floor(hours(b))*60),...
                floor(seconds(b)-floor(minutes(b))*60));
        else
            B = datetime(year(eegdate),month(eegdate),day(eegdate),...
                floor(hours(b)),floor(minutes(b)-floor(hours(b))*60),...
                floor(seconds(b)-floor(minutes(b))*60));
        end
        
        N2a = T_sleep.GoodN2SleepStart(tableidx(ientry));
        N2b = T_sleep.GoodN2SleepStop(tableidx(ientry));
        
        if N2a < duration('12:00:00') % if N2 sleep starts in the morning
            % datetime for N2 sleep start: add one day because it's after
            % midnight (next morning)
            N2A = datetime(year(eegdate),month(eegdate),day(eegdate)+1,...
                floor(hours(N2a)),minutes(N2a)-floor(hours(N2a))*60,0);
        else
            % datetime for N2 sleep start
            N2A = datetime(year(eegdate),month(eegdate),day(eegdate),...
                floor(hours(N2a)),minutes(N2a)-floor(hours(N2a))*60,0);
        end
        
        if N2b < duration('12:00:00') % if N2 sleep ends in the morning
            % datetime for N2 sleep end: add one day because it's after
            % midnight (next morning)
            N2B = datetime(year(eegdate),month(eegdate),day(eegdate)+1,...
                floor(hours(N2b)),minutes(N2b)-floor(hours(N2b))*60,0);
        else
            % datetime for N2 sleep end
            N2B = datetime(year(eegdate),month(eegdate),day(eegdate),...
                floor(hours(N2b)),minutes(N2b)-floor(hours(N2b))*60,0);
        end
        
        % wake start/end times (assumes child awake 5:00 - 5:30 pm)
        wakeA = datetime(year(eegdate),month(eegdate),day(eegdate),...
            str2double(wakestart(1:2)),str2double(wakestart(4:5)),0);
        wakeB = datetime(year(eegdate),month(eegdate),day(eegdate),...
            str2double(wakeend(1:2)),str2double(wakeend(4:5)),0);
        
        assert(wakeB < N2A,'N2 sleep starts before awake peroid')
        
        % if this file contains BOTH the period of wake and N2 sleep
        if A < wakeA && B > N2B
            condition = 'both';
        % if this file starts before 5:00 pm and ends after 5:30 pm
        elseif A < wakeA &&  B > wakeB
            condition = 'wake';
        % if this file contains the period likely to be N2 sleep as judged
        % by spindles
        elseif A < N2A && B > N2B
            condition = 'sleep';
        end
        
        
        % if this file is either good for sleep or wake
        if ~strcmp(condition,'neither')
            importme = T_sleep.SegIDs{tableidx(ientry)};
            try
                [eegdata,header] = lab_read_edf(sprintf('%s/%s_1-1+.edf',destination,importme));
                %OUTEEG = pop_biosig(sprintf('%s/NKT/EEG2100/%s.eeg',destination,importme)); % import using biosig
            catch
                fprintf('Failed to load %s_1-1+.edf for %s condition... skipping this one\n',importme,condition)
                keyboard
                continue
            end
            
            %%% Check that times match %%%
            edf_time = duration(header.hour,header.minute,header.second);
            assert(edf_time == T_sleep.StartTime(tableidx(ientry)),'Wrong time');
            
            tmpdata.trial{1} = eegdata;
            tmpdata.time{1} = linspace(0,size(eegdata,2)/header.samplingrate,size(eegdata,2));
            for ich = 1:size(header.channels,1)
                tmpdata.label{ich} = char(strtrim(upper(header.channels(ich,:))));
            end
            tmpdata.fsample = header.samplingrate;
            clear eegdata header
          
            %%% extract section of interest
            switch condition
                case 'wake'
                    % number of seconds from beginning of recording until
                    % beginning of the (assumed) wake epoch beginning at 5 pm
                    nsec = seconds(duration(wakestart) - T_sleep.StartTime(tableidx(ientry)));
                    nsamp = nsec*tmpdata.fsample; % convert seconds to samples
                    
                    % number of seconds from beginning of recording until
                    % the end of the (assumed) wake epoch (5:30 pm)
                    msec = seconds(duration(wakeend) - T_sleep.StartTime(tableidx(ientry)));
                    % convert seconds to samples (plus one because indexing doesn't start from zero)
                    msamp = msec*tmpdata.fsample+1;
                    
                    tmpdata.trial{1} = tmpdata.trial{1}(:,nsamp:msamp);
                case 'sleep'
                    % number of seconds from beginning of recording until
                    % beginning of the N2 sleep section
                    nsec = seconds(T_sleep.GoodN2SleepStart(tableidx(ientry)) ....
                        - T_sleep.StartTime(tableidx(ientry)));
                    if nsec < 0 % if we switched from PM to AM during recording causing this bug
                        nsec = seconds(T_sleep.GoodN2SleepStart(tableidx(ientry))) ...
                            + (60*60*24 - seconds(T_sleep.StartTime(tableidx(ientry))));
                    end
                    nsamp = nsec*tmpdata.fsample; % convert seconds to samples
                    
                    % number of seconds from beginning of recording until
                    % the end of the N2 sleep section
                    msec = seconds(T_sleep.GoodN2SleepStop(tableidx(ientry)) ....
                        - T_sleep.StartTime(tableidx(ientry)));
                    if msec < 0 % if we switched from PM to AM during recording causing this bug
                        msec = seconds(T_sleep.GoodN2SleepStop(tableidx(ientry))) ...
                            + (60*60*24 - seconds(T_sleep.StartTime(tableidx(ientry))));
                    end
                    % convert seconds to samples (plus one because indexing doesn't start from zero)
                    msamp = msec*tmpdata.fsample+1;
                    
                    tmpdata.trial{1} = tmpdata.trial{1}(:,nsamp:msamp);
                case 'both'
                    
                    %%% AWAKE %%%
                    
                    % number of seconds from beginning of recording until
                    % beginning of the (assumed) wake epoch beginning at 5 pm
                    nsec = seconds(duration(wakestart) - T_sleep.StartTime(tableidx(ientry)));
                    nsamp = nsec*tmpdata.fsample; % convert seconds to samples
                    
                    % number of seconds from beginning of recording until
                    % the end of the (assumed) wake epoch (5:30 pm)
                    msec = seconds(duration(wakeend) - T_sleep.StartTime(tableidx(ientry)));
                    if msec < 0 % if we switched from PM to AM during recording causing this bug
                        msec = seconds(T_sleep.GoodN2SleepStop(tableidx(ientry))) ...
                            + (60*60*24 - seconds(T_sleep.StartTime(tableidx(ientry))));
                    end
                    % convert seconds to samples (plus one because indexing doesn't start from zero)
                    msamp = msec*tmpdata.fsample+1;
                    
                    tmpdata.trial{1} = tmpdata.trial{1}(:,nsamp:msamp);
                    
                    clear nsec nsamp msec msamp
                    
                    %%% ASLEEP %%%
                    
                    % number of seconds from beginning of recording until
                    % beginning of the N2 sleep section
                    nsec = seconds(T_sleep.GoodN2SleepStart(tableidx(ientry)) ....
                        - T_sleep.StartTime(tableidx(ientry)));
                    if nsec < 0 % if we switched from PM to AM during recording causing this bug
                        nsec = seconds(T_sleep.GoodN2SleepStart(tableidx(ientry))) ...
                            + (60*60*24 - seconds(T_sleep.StartTime(tableidx(ientry))));
                    end
                    nsamp = nsec*tmpdata.fsample; % convert seconds to samples
                    
                    % number of seconds from beginning of recording until
                    % the end of the N2 sleep section
                    msec = seconds(T_sleep.GoodN2SleepStop(tableidx(ientry)) ....
                        - T_sleep.StartTime(tableidx(ientry)));
                    % convert seconds to samples (plus one because indexing doesn't start from zero)
                    msamp = msec*tmpdata.fsample+1;
                    
                    tmpdata.trial{2} = tmpdata.trial{1}(:,nsamp:msamp);
                otherwise
                    error('Condition not recognized')
            end
            
            for itrl = 1:length(tmpdata.trial)
            
                label_in = tmpdata.label;
                dat = tmpdata.trial{itrl};
                fsample = tmpdata.fsample;
                if ~exist('data','var')
                    data = tmpdata;
                    data.cfg.raw.ctype.EEG.num = []; % initiate
                    data.cfg.raw.ctype.EEG.label = []; % initiate
                end

                % detect channels, reorder the EEG channels if neccessary
                [~,where] = ismember(upper(label_load),upper(label_in));

                % are we missing channels?
                nchan = sum(where>0);
                if nchan <= 16
                    % try alternative labels
                    [~,where] = ismember(upper(label_load_alt),upper(label_in));
                    nchan = sum(where>0);
                end
                
                if nchan < 19 && nchan > 16 % if number of missing channels is between 1 and 3

                    % reorder channels and replace missing channels with
                    % junk (can't use NaNs before filtering)
                    dat2 = ones(length(label_load),size(dat,2)).*1e6; % junk
                    for irow = 1:size(dat2,1)
                        if where(irow) ~= 0
                            dat2(irow,:) = dat(where(irow),:); % reorder channel
                        end
                    end
                    dat = dat2; % overwrite data

                    % which channels are missing?
                    [missing,missing_ndx] = setdiff(upper(label_load),upper(label_in)); % grab the index from label_load
                    % check if there were already missing channels from other trials
                    if ~isfield(data.cfg.raw.ctype,'bad') || isempty(data.cfg.raw.ctype.bad.num)
                        data.cfg.raw.ctype.bad.num = missing_ndx;
                        data.cfg.raw.ctype.bad.label = missing;
                    else
                        assert(all(data.cfg.raw.ctype.bad.num == missing_ndx)...
                            ,'missing channels don''t match with previous trials')
                        assert(all(strcmp(data.cfg.raw.ctype.bad.label, missing))...
                            ,'missing channels don''t match with previous trials')
                    end

                    % which channels are present?
                    [present,present_ndx] = intersect(upper(label_load),upper(label_in));
                    % check earlier trials
                    if isempty(data.cfg.raw.ctype.EEG.num)
                        data.cfg.raw.ctype.EEG.num = present_ndx;
                        data.cfg.raw.ctype.EEG.label = present;
                    else
                        assert(all(data.cfg.raw.ctype.EEG.num == present_ndx)...
                            ,'present channels don''t match with previous trials')
                        assert(all(strcmp(data.cfg.raw.ctype.EEG.label, present))...
                            ,'present channels don''t match with previous trials')
                    end

                elseif nchan <= 16
                    error('Too many missing channels')
                else
                    % reorder channels
                    dat = dat(where,:);

                    % which channels are present?
                    [present,present_ndx] = intersect(upper(label_load),upper(label_in));
                    % check earlier trials
                    if isempty(data.cfg.raw.ctype.EEG.num)
                        data.cfg.raw.ctype.EEG.num = present_ndx;
                        data.cfg.raw.ctype.EEG.label = present;
                    else
                        assert(all(data.cfg.raw.ctype.EEG.num == present_ndx)...
                            ,'present channels don''t match with previous trials')
                        assert(all(strcmp(data.cfg.raw.ctype.EEG.label, present))...
                            ,'present channels don''t match with previous trials')
                    end
                end

                time = (0:size(dat,2)-1)/tmpdata.fsample;
                label = label_load;

                % some files contain single, isolated NaNs -- these can be
                % removed using interp1 (spline)
                % Do each row seperately
                tsamp = 1:size(dat,2);
                for irow = 1:size(dat,1)
                    idx = ~isnan(dat(irow,:));
                    dat(irow,~idx) = interp1(tsamp(idx),dat(irow,idx),tsamp(~idx),'spline');
                end

                %%% Remove **mean** from signal**

                % Comment out line below if you don't want to do the detrending
                % Detrend the data to correct for DC offsets
                % Must use mean (median doesn't remove them very well)
                dat = single(dat) - repmat(mean(double(dat)')',1,length(dat)); % use nanmean if you get a bug

                %%% Remove **linear trend** from signal
                % not sure it has to be done in a loop, but detrend may prefer vectors over
                % matricies

                for ich = 1:size(dat,1)
                    % alternatively use mydetrend to handle nans if you get a bug
                    dat(ich,:) = detrend(dat(ich,:)); % doesn't work with NaNs
                end

                [where] = contains(upper(label_in),upper(label_physio));
                if sum(where>0)~=4
                    [where] = contains(upper(label_in),upper(label_physio_alt));
                end

                idx_chan(where) = true;
                NDX = idx_chan(where);
                ECG = [];
                EOG = [];

                for i = 1:length(NDX)
                    if contains(label_in(NDX(i)),'EKG') || contains(label_in(NDX(i)),'ECG')
                        if ~exist('ecg1','var')
                            ecg1 = data(NDX(i));
                            ECG = [ECG; ecg1];
                        else
                            ecg2 = data(NDX(i));
                            ECG = [ECG; ecg2];
                        end
                    elseif contains(label_in(NDX(i)),'EOG') || contains(label_in(NDX(i)),'LOC') ...
                            || contains(label_in(NDX(i)),'ROC')
                        if ~exist('eog1','var')
                            eog1 = data(NDX(i));
                            EOG = [EOG; eog1];
                        else
                            eog2 = data(NDX(i));
                            EOG = [EOG; eog2];
                        end
                    end
                end

                physio = [ECG; EOG];


                %% Do filtering

                % Define Filter Coefficients (see preproc_parameter.m)
                for idef=1:length(filter_def), eval(filter_def{idef}), end

                
                % Highpass Butterworth filter (better for removing drift)
                % filtfilt should correct for IIR phase distortion (double check) 
                % Note that the slow rolloff shouldn't be an issue using 0.4 Hz highpass,
                % the attenuation at 0.5 Hz is only 0.44 dB
                % (see Compare_filter_responses.m)

                % Lowpass filter using FIR 
                forder = 2; % FIR order 2x sampling rate
                lpFilt = designfilt('lowpassfir','cutoffFrequency',f_lp, ...
                     'filterOrder',fsample*forder,'sampleRate',fsample);

                tmp_hp = single(filtfilt(Bhp, Ahp, double(dat)')'); % 5th order Butterworth highpass
                dat = single(filtfilt(lpFilt, double(tmp_hp)')'); % FIR lowpass
                
                dat_hp = tmp_hp - dat;

                for icnt=1:1, dat_hp = single(filtfilt(B_line{icnt}, A_line{icnt}, double(dat_hp)')'); end % only notch 60 Hz to avoid numerical problems

                if exist('eog2','var')
                    eog2 = filtfilt(B_eog, A_eog, eog2')';
                end
                if exist('eog1','var')
                    eog1 = filtfilt(B_eog, A_eog, eog1')';
                end
                if exist('ecg1','var')
                    ecg1 = filtfilt(B_ecg, A_ecg, ecg1')';
                end
                if exist('ecg2','var')
                    ecg2 = filtfilt(B_ecg, A_ecg, ecg2')';
                end

                % average reference
                dat    = dat    - repmat(mean(dat,1),[size(dat,1),1]); % avg reference
                dat_hp = dat_hp - repmat(mean(dat_hp,1),[size(dat_hp,1),1]); % avg reference

                % get the gamma-band envelope (helps manual artifact rejection)
                gamma = nan(size(dat_hp));
                if mod(fnyq,2) == 0 % if even
                    n = fnyq + 1; % make odd number
                else
                    warning('Sampling rate is odd number!')
                    fprintf('\nDetected sampling rate = %i Hz\n',fsample)
                    n = fnyq;
                end
                N = (n-1)/2;
                for k = 1:size(dat_hp,1)
                    [up] = abs(hilbert(dat_hp(k,:))); % OLD CODE: [up,~] = envelope(dat_hp(k,:));
                    smoothEnv = conv(up,ones(1,n)./n); % smooth the envelope squared (inst. power)
                    gamma(k,:) = smoothEnv(N+1:end-N); % make same size as signal
                end

                %% Put into Fieldtrip structure

                switch condition
                    case 'wake'
                        data.time{1}      = single(time);
                        data.trial{1}     = single(dat);
                        data.orgfile{1}   = importme;
                    case 'sleep'
                        data.time{2}      = single(time);
                        data.trial{2}     = single(dat);
                        data.orgfile{2}   = importme;
                    case 'both'
                        data.time{itrl}      = single(time);
                        data.trial{itrl}     = single(dat);
                        data.orgfile{itrl}   = importme;
                end

                if ~isfield(data,'genotype')
                    data.genotype   = T_sleep.DupType{tableidx(ientry)};
                    data.age        = T_sleep.AgeAtSleepEEG(tableidx(ientry));
                    data.eegdate    = T_sleep.EEGDate(tableidx(ientry));
                    data.epilepsy   = T_sleep.EpilpesyStatus{tableidx(ientry)};
                    data.N2_start   = T_sleep.GoodN2SleepStart(tableidx(ientry));
                    data.N2_stop    = T_sleep.GoodN2SleepStop(tableidx(ientry));
                    data.wake_start = duration(wakestart);
                    data.wake_stop  = duration(wakeend);
                end

                if exist('eog1','var')
                    label{end+1} = 'eog1';
                    data.trial{ifile}(end+1,:) = single(eog1);
                end
                if exist('eog2','var')
                    label{end+1} = 'eog2';
                    data.trial{ifile}(end+1,:) = singel(eog2);
                end
                if exist('ECG1','var')
                    label{end+1} = 'ECG1';
                    data.trial{ifile}(end+1,:) = single(ECG1);
                end
                if exist('ECG2','var')
                    label{end+1} = 'ECG2';
                    data.trial{ifile}(end+1,:) = single(ECG2);
                end

                label{end+1} = 'GAMMA';

                % we've already taken the wave envelope so we can average without
                % cancelation
                switch condition
                    case 'wake'
                        data.trial{1}(end+1,:) = single(max(gamma));
                        %data.dat_hp{1}         = single(dat_hp);
                    case 'sleep'
                        data.trial{2}(end+1,:) = single(max(gamma));
                        %data.dat_hp{2}         = single(dat_hp);
                    case 'both'
                        data.trial{itrl}(end+1,:) = single(max(gamma));
                        %data.dat_hp{itrl}         = single(dat_hp);
                end

                if ~isfield(data,'site')
                    data.label        = label;
                    data.ID           = ID;
                    data.orgForm      = 'EDF';

                    %data.cfg.event    = event;
                    data.elec         = elec;
                    data.dimord       = 'chan_time';
                    %data.startdate    = datnum;

                    % Uncomment below if you have trouble with channel names not matching
                    % lay file
                    %data.label        = elec.label';

                    data.site = 'UCLA';

                    % add units and labels for auxillary channels
                    data.cfg.raw.ctypelabels      = {'EEG','physio','bad'};
                    data.cfg.raw.ctypemainlabel   = 'EEG';
                    data.cfg.raw.ctypephysiolabel = 'physio';
                    for icnt=1:length(label)
                        data.chaninfo.unit{icnt} = 'uV';
                        data.chaninfo.type{icnt} = 'EEG';
                    end
                    % carful, check the number of physio channels for each preproc script!
                    data.chaninfo.unit{icnt+1} = 'uV';
                    data.chaninfo.unit{icnt+2} = 'uV';
                    data.chaninfo.unit{icnt+3} = 'uV';
                    data.chaninfo.unit{icnt+4} = 'uV';
                    data.chaninfo.unit{icnt+5} = 'uV';
                    data.chaninfo.type{icnt+1} = 'physio';
                    data.chaninfo.type{icnt+2} = 'physio';
                    data.chaninfo.type{icnt+3} = 'physio';
                    data.chaninfo.type{icnt+4} = 'physio';
                    data.chaninfo.type{icnt+5} = 'physio';
                else
                    fprintf('Data fields have already been filled from first file\n')
                end
            end
        end
    end
    
    if exist('data','var') && length(data.trial)==2
        sampleinfo = nan(length(data.trial),2);
        trltrl = nan(length(data.trial),3);

        for itrl = 1:length(data.trial)
            if itrl == 1
                sampleinfo = [1,size(data.trial,2)];
                trltrl(itrl,:) = [1,size(data.trial{itrl},2),sampleinfo(itrl,1)];
            else
                sampleinfo(itrl,:) = [(1 + size(data.trial{itrl-1},2)),(size(data.trial{itrl-1},2) ...
                    + size(data.trial{itrl},2))];

                trltrl(itrl,:) = [(1 + size(data.trial{itrl-1},2)),(size(data.trial{itrl-1},2) ...
                    + size(data.trial{itrl},2)), sampleinfo(itrl,1)];
            end
        end

        data.cfg.trl.trl = trltrl;
        data.sampleinfo = sampleinfo;

        % extra channels get added as physio channels
        if length(data.label) > length(label_load)
            data.cfg.raw.ctype.physio.num = length(label_load)+1:length(data.label);
            data.cfg.raw.ctype.physio.label = data.label{length(label_load)+1:end};
        end

        %     % fill if no missing channels
        %     if isempty(data.cfg.raw.ctype.bad.num)
        %         data.cfg.raw.ctype.bad.num = [];
        %         data.cfg.raw.ctype.bad.label = cell(1,1);
        %     end

        % save data
        outname = sprintf('%sDup15q_%s_%ih%im.mat',DIRRESULT,ID,str2double(wakestart(1:2)),str2double(wakestart(4:5)));
        fprintf('\nSaving file %s\n',outname)
        save(outname,'data','-v7.3')
        
        return

        % do table entries
        if isfield( data.cfg.raw.ctype,'bad')
            lostchans{isub} = data.cfg.raw.ctype.bad.label; % missing channels
        else
            lostchans{isub} = cell(1,1);
        end
        ntrials{isub} = length(data.trial);
        whichtrials{isub} = data.orgfile; % which trials were used
        len = 0;
        for itrl = 1:ntrials{isub}
            len = len + size(data.trial{itrl},2)/fsample;
        end
        datalen{isub} = len; % total data length in sec
    else
        fprintf('Skipping this subject, missing data ... \n')
    end

    % clear variables from loop
    clearvars -except DIRRESULT f_hp_ecg filter_def lay FList f_hp_eog ...
        label_alt n_butt f_lp label_load n_butt_EEG elec f_lp_ecg label_physio ...
        n_times_fsample_fir f_hp f_lp_eog label_physio_alt start_ndx subs ...
        lostchans ntrials whichtrials datalen T_sleep wakestart wakeend
end

try
    % write table
    T = table(subs',lostchans,ntrials,whichtrials,datalen);
    writetable(T,'Dup15q_data_import.csv')
catch
    keyboard
end


end







