function [] = ana_preproc_as_boston_EDF_C_butterHP_firLP()

% Load and pre-process the AS data from Boston in EDF_C format

load angelman_lay.mat elec

%% Defines parameter for preprocessing centrally

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
filter_def{5} = '[Bhp, Ahp] = butter(n_butt_EEG,[f_hp]/fnyq,''high'');';
filter_def{6} = '[Blp, Alp] = butter(n_butt_EEG,[f_lp]/fnyq,''low'');';


DIRRESULT = './Monti/AS_butterHP_firLP/Boston_EDF_C/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end
DIRDATA = './data_raw/boston/October_EDF_C/';

%%
subject = dir([DIRDATA,'*.EDF']);

label_load = {'C3','C4','O1','O2','Cz','F3','F4','F7','F8','Fz','FP1','FP2','P3','P4','Pz','T7','T8','P7','P8'};
label_alt =  {'EC3REF','EC4REF','EO1REF','EO2REF','ECzREF','EF3REF','EF4REF','EF7REF','EF8REF','EFzREF',...
    'EFP1REF','EFP2REF','EP3REF','EP4REF','EPzREF','ET7REF','ET8REF','EP7REF','EP8REF'};
label_physio = {'ECGL','ECGR','LUE','RAE'};
label_physio_alt = {'PECGL','PECGR','PLUE','PRAE'};

for isub=1:length(subject)
    ID = subject(isub).name(1:6);
    filename = [DIRDATA,subject(isub).name];
    try
        [~,hdr]=lab_read_edf_jfh(filename); % allows reading the event information and other header information propperl
        assert(~isempty(hdr))
    catch
        fprintf('Skipping this one, can''t open file ... \n')
        continue
    end
    fsample=hdr.samplingrate;
    EVENTS=hdr.events;
    % create empty cell arrays
    clear event
    event.duration = cell(1,length(EVENTS.DUR));
    event.value = cell(1,length(EVENTS.TYP));
    event.sample = cell(1,length(EVENTS.POS));
    for ievt = 1:length(EVENTS.TYP)
        event(ievt).duration = EVENTS.DUR(ievt);
        event(ievt).value    = EVENTS.TYP{ievt};
        event(ievt).sample   = EVENTS.POS(ievt);
    end
    % sanity check
    if ~strcmp(ID,hdr.subject.name(1:6))
        isub
        filename
        %keyboard
        warning('check!')
    end
    [hdr2,data] = edfread(filename); % allows reading the data with correct scaling
    label_in = hdr2.label;
    datnum = datenum([hdr2.startdate,'.',hdr2.starttime],'dd.mm.yy.HH.MM.SS');
    
    % reorder the EEG channels if neccessary
    [~,where] = ismember(upper(label_in),upper(label_load));
    if sum(where>0)~=19
        [~,where] = ismember(upper(label_in),upper(label_alt));
    end
    
    
    assert(sum(where>0)==19,'missing channels')
    clear idx_chan 
    % reorder channels
    for icnt=1:length(label_load)
        idx_chan(icnt) = find(where==icnt);
    end
    dat = data(idx_chan,:);
    time = (0:size(dat,2)-1)/fsample;
    label = label_load;
    
    % Load EOG
    clear physio ecg1 ecg2 eogV eogH
    %physio = nan(4,length(data.trial{1}); % memory allocation
    [~,where] = ismember(upper(label_in),upper(label_physio));
    if sum(where>0)~=4
        [~,where] = ismember(upper(label_in),upper(label_physio_alt));
    end
    assert(sum(where>0)==4,'physio channels not found??')
    clear idx_chan % @Joel: make sure the channel order is correct!!! -> apply to other preproc scripts
    for icnt=1:length(label_physio)
        idx_chan(icnt) = find(where==icnt);
    end
    physio = data(idx_chan,:);
    ecg1=physio(1,:);
    ecg2=physio(2,:);
    eogV=physio(3,:);
    eogH=physio(4,:);
    
    % Comment out line below if you don't want to do the detrending
    % Detrend the data to correct for DC offsets
    % Must use mean (median doesn't remove them very well)
    dat = single(dat) - repmat(mean(double(dat)')',1,length(dat));

    % not sure it has to be done in a loop, but detrend may prefer vectors over
    % matricies

    for ich = 1:size(dat,1)
        dat(ich,:) = detrend(dat(ich,:));
    end


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

    dat_hp = single(filtfilt(Bhp, Ahp, double(dat)')'); % 5th order Butterworth highpass
    dat = single(filtfilt(lpFilt, double(dat_hp)')'); % FIR lowpass

    dat_gmm = dat_hp-dat; % subtract the bandpassed data from the higpassed data
    dat_gmm = single(filtfilt(B_line{icnt}, A_line{icnt}, double(dat_gmm)')'); % only notch 60 Hz to avoid numerical problems

    eogH = filtfilt(B_eog, A_eog, eogH')';
    eogV = filtfilt(B_eog, A_eog, eogV')';
    ecg1 = filtfilt(B_ecg, A_ecg, ecg1')';
    ecg2 = filtfilt(B_ecg, A_ecg, ecg2')';
    
    % average reference
    dat    = dat    - repmat(mean(dat,1),[size(dat,1),1]); % avg reference
    dat_gmm = dat_gmm - repmat(mean(dat_gmm,1),[size(dat_gmm,1),1]); % avg reference
    
    % get the gamma-band envelope (helps manual artifact rejection)
    gamma = nan(size(dat_gmm));
    if mod(fnyq,2) == 0 % if even
        n = fnyq + 1; % make odd number
    else
        warning('Sampling rate is odd number!')
        n = fnyq;
    end
    N = (n-1)/2;
    for k = 1:size(dat_gmm,1)
        [up] = abs(hilbert(dat_gmm(k,:))); % OLD CODE: [up,~] = envelope(dat_gmm(k,:));
        smoothEnv = conv(up,ones(1,n)./n); % smooth the envelope squared (inst. power)
        gamma(k,:) = smoothEnv(N+1:end-N); % make same size as signal
    end
    
    % Put into Fieldtrip structure
    clear data
    data.ID           = ID;
    data.cfg.event    = event;
    data.elec         = elec;
    %data.dat_gmm       = single(dat_gmm);
    data.dimord       = 'chan_time';
    data.startdate    = datnum;
    data.fsample      = fsample;
    % use these labels so they match the lay file!
    data.label        = elec.label'; 
    data.time{1}      = time;
    data.label{end+1} = 'EOGV';
    data.label{end+1} = 'EOGH';
    data.label{end+1} = 'ECG1';
    data.label{end+1} = 'ECG2';
    data.label{end+1} = 'GAMMA';
    data.orgForm      = 'EDF_C';
    data.time{1}      = single(time);
    data.trial{1}     = single(dat);
    data.trial{1}(end+1,:) = single(eogV);
    data.trial{1}(end+1,:) = single(eogH);
    data.trial{1}(end+1,:) = single(ecg1);
    data.trial{1}(end+1,:) = single(ecg2);
    % we've already taken the wave envelope so we can average without
    % cancelation
    data.trial{1}(end+1,:) = single(max(gamma));
    data.sampleinfo(1,:) = [1,size(dat,2)];
    data.cfg.trl.trl(1,:) = [1,size(dat,2),0];
    data.site = 'BOS';
    
    % add units and labels for auxillary channels 
    data.cfg.raw.ctypelabels      = {'EEG','physio'};
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
    
    outname = sprintf('%sAS_%s_%s.mat',DIRRESULT,ID,datestr(data.startdate,'yyyymmdd'));
    save(outname,'data','-v7.3')
end
