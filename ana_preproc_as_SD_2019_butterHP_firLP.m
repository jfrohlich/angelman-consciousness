function [] = ana_preproc_as_SD_2019_butterHP_firLP()

% Load and pre-process the AS data from Boston

%addpath('./AS_overnight/');
%addfieldtrip

load angelman_lay.mat lay elec

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

DIRDATA = './data_raw/san_diego_2019/';
DIRRESULT = './Monti/AS_butterHP_firLP/san_diego_2019/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end

%%

tmp = dir(sprintf('%s*.edf',DIRDATA));
cnt=0; clear subjects
for icnt=1:length(tmp)
    idx = strfind(tmp(icnt).name,'_');
    if length(idx)<2
        cnt=cnt+1;
        subjects{cnt} = tmp(icnt).name(1:end-4);
    elseif ~isempty(strfind(tmp(icnt).name,'_0001'))
        cnt=cnt+1;
        subjects{cnt} = tmp(icnt).name(1:idx(2)-1);
    end
end
for isubj=1:length(subjects)
    fprintf('%s\n',subjects{isubj})
end

label_load = {'C3','C4','O1','O2','Cz','F3','F4','F7','F8','Fz','FP1','FP2','P3','P4','Pz','T7','T8','P7','P8'};
label_alt =  {'C3','C4','O1','O2','Cz','F3','F4','F7','F8','Fz','FP1','FP2','P3','P4','Pz','T3','T4','T5','T6'};
label_physio = {'ECGL','ECGR','LUE','RAE'};
label_physio_alt = {'ECGL','ECGR','LOC','ROC'};

start_ndx = 1;

if start_ndx ~= 1, warning('Starting index is set to %i\n',start_ndx), end

%for isub = start_ndx:length(subjects)
for isub = length(subjects):-1:start_ndx
    fname = dir(sprintf('%s%s*.edf',DIRDATA,subjects{isub}));
    fprintf('Working on %s%s\n',DIRDATA,subjects{isub})
    if ~isempty(fname)
        clear cut datnum event
        for ifile=1:length(fname)
            ID = fname(ifile).name(1:6);
            filename = fname(ifile).name;
            [~,hdr]=lab_read_edf_jfh([DIRDATA,filename]); % reading the event information and other header information propperly
            [hdr2,data_tmp] = edfread([DIRDATA,filename]);
            if ifile==1
                fsample=hdr.samplingrate;
            else
                if fsample~=hdr.samplingrate; error('samplingrate'), end
            end
            
            % Load EEG data
            if ifile==1 % if we're taking the first chunk
                [hdr2,data] = edfread([DIRDATA,filename]);
                label_in = hdr2.label;
                datnum(ifile) = datenum([hdr2.startdate,'.',hdr2.starttime],'dd.mm.yy.HH.MM.SS');
                cut(ifile)=size(data,2);
            else % if we're not on the first data chunk for that subject
                data=cat(2,data,data_tmp);
                datnum(ifile) = datenum([hdr2.startdate,'.',hdr2.starttime],'dd.mm.yy.HH.MM.SS');
                cut(ifile)=size(data,2);
            end
            
            % Load Annotations
            if isfield(hdr,'events')
                if ~exist('event','var') % if we're taking the first chunk
                    EVENTS=hdr.events;
                    % create empty cell arrays
                    clear event
                    event.duration = cell(1,length(EVENTS.DUR));
                    event.value = cell(1,length(EVENTS.TYP));
                    event.sample = cell(1,length(EVENTS.POS));
                    for ievt = 1:length(EVENTS.TYP)
                        event(ievt).duration = EVENTS.DUR(ievt);
                        event(ievt).value    = EVENTS.TYP{ievt};
                        event(ievt).type     = EVENTS.TYP{ievt};
                        event(ievt).sample   = EVENTS.POS(ievt);
                    end
                else % if we're not on the first data chunk for that subject
                    EVENTS=hdr.events;
                    % add to events to subfields
                    len = length(event);
                    for ievt = 1:length(EVENTS.TYP)
                        event(len+ievt).duration = EVENTS.DUR(ievt);
                        event(len+ievt).value    = EVENTS.TYP{ievt};
                        event(len+ievt).type     = EVENTS.TYP{ievt};
                        %event(len+ievt).sample   =
                        %size(data,2)+EVENTS.POS(ievt); BUG!
                        event(len+ievt).sample   = cut(ifile-1)+EVENTS.POS(ievt);
                    end
                end
            end
        end
        
        % Re-order the EEG channels if neccessary
        [~,where] = ismember(upper(label_in),upper(label_load));
        if sum(where>0)~=19
            [~,where] = ismember(upper(label_in),upper(label_alt));
        end
        assert(sum(where>0)==19,'missing channels')
        clear idx_chan
        for icnt=1:length(label_load)
            idx_chan(icnt) = find(where==icnt);
        end
        dat = data(idx_chan,:);
        time = (0:size(dat,2)-1)/fsample;
        label = label_load;
        
        % Load ECG signals (if available)
        for ichan_in=1:length(label_in)
            isphysio(ichan_in) = ~isempty([strfind(upper(label_in{ichan_in}),'ECG'),strfind(upper(label_in{ichan_in}),'EKG')]);
        end
        idx_physio = find(isphysio);
        ecg = data(idx_physio,:);
        label_ecg = label_in(idx_physio);
        
        % Load EOG signals (if available)
        for ichan_in=1:length(label_in)
            isphysio(ichan_in) = ~isempty([strfind(upper(label_in{ichan_in}),'EOG')]);
        end
        idx_physio = find(isphysio);
        eog = data(idx_physio,:);
        label_eog = label_in(idx_physio);
        
        % Detrend
        dat = single(dat) - repmat(mean(double(dat)')',1,length(dat));
                
        % Define Filter Coefficients (see preproc_parameter.m)
        for idef=1:length(filter_def), eval(filter_def{idef}), end
        
        % Highpass Butterworth filter (better for removing drift)
        % filtfilt should correct for IIR phase distortion (double check)
        % Note that the slow rolloff shouldn't be an issue using 0.4 Hz highpass,
        % the attenuation at 0.5 Hz is only 0.44 dB
        % (see Compare_filter_responses.m)
        
        % Lowpass filter using FIR
        forder = 2; % FIR order 2x sampling rate
        
        try
            lpFilt = designfilt('lowpassfir','cutoffFrequency',f_lp, ...
                'filterOrder',fsample*forder,'sampleRate',fsample);
        catch
            keyboard
        end
        
        dat_hp = single(filtfilt(Bhp, Ahp, double(dat)')'); % 5th order Butterworth highpass
        dat = single(filtfilt(lpFilt, double(dat_hp)')'); % FIR lowpass
        
        dat_gmm = dat_hp-dat; % subtract the bandpassed data from the higpassed data
        dat_gmm = single(filtfilt(B_line{icnt}, A_line{icnt}, double(dat_gmm)')'); % only notch 60 Hz to avoid numerical problems
                
        % Average reference
        dat = dat-repmat(mean(dat,1),[size(dat,1),1]); % avg reference
        dat_gmm = dat_gmm-repmat(mean(dat_gmm,1),[size(dat_gmm,1),1]); % avg reference
        
        
        % get the gamma-band envelope (helps manual artifact rejection)
        gamma = nan(size(dat_hp));
        if mod(fnyq,2) == 0 % if even
            n = fnyq + 1; % make odd number
        else
            warning('Sampling rate is odd number!')
            n = fnyq;
        end
        N = (n-1)/2;
        for k = 1:size(dat_hp,1)
            [up] = abs(hilbert(dat_hp(k,:))); % OLD CODE: [up,~] = envelope(dat_hp(k,:));
            smoothEnv = conv(up,ones(1,n)./n); % smooth the envelope squared (inst. power)
            gamma(k,:) = smoothEnv(N+1:end-N); % make same size as signal
        end
        
        % Set FIR filter length arround the data fusion points to nan
        cut = [0,cut];
        for icnt=2:length(cut)
            idx_start = max(cut(icnt)-fsample*forder,1);
            idx_end = min(cut(icnt)+fsample*forder,size(dat,2));
            idx = idx_start:idx_end;
            dat(:,idx)=nan;
            dat_hp(:,idx)=nan;
            gamma(:,idx)=nan;
            eog(:,idx)=nan;
            ecg(:,idx)=nan;
        end
        
        % Put into Fieldtrip structure
        clear data
        data.cfg.event    = event;
        data.ID           = str2num(ID);
        data.elec         = elec;
        %data.dat_hp       = single(dat_hp);
        data.dimord       = 'chan_time';
        data.startdate    = datnum;
        data.fsample      = fsample;
        data.label        = label_alt;
        data.time{1}      = single(time);
        data.trial{1}     = single(dat);
        % physio channels
        data.label{end+1} = 'GAMMA';
        data.trial{1}(end+1,:) = single(max(gamma));
        for icnt=1:length(label_eog)
            data.label{end+1} = label_eog{icnt};
            data.trial{1}(end+1,:) = single(eog(icnt,:));
        end
        for icnt=1:length(label_ecg)
            data.label{end+1} = label_ecg{icnt};
            data.trial{1}(end+1,:) = single(ecg(icnt,:));
        end
        data.orgForm      = 'EDF';
        data.sampleinfo(1,:)   = [1,size(dat,2)];
        data.cfg.trl.trl(1,:)  = [1,size(dat,2),0];
        data.site = 'SD';
        % add units and labels for auxillary channels
        data.cfg.raw.ctypelabels      = {'EEG','physio'};
        data.cfg.raw.ctypemainlabel   = 'EEG';
        data.cfg.raw.ctypephysiolabel = 'physio';
        for icnt=1:length(label)
            data.chaninfo.unit{icnt} = 'uV';
            data.chaninfo.type{icnt} = 'EEG';
        end
        icnt=icnt+1; % "EMG"
        data.chaninfo.unit{icnt} = 'uV';
        data.chaninfo.type{icnt} = 'physio';
        for icnt2=1:length(label_eog)
            icnt=icnt+1;
            data.chaninfo.unit{icnt} = 'uV';
            data.chaninfo.type{icnt} = 'physio';
        end
        for icnt2=1:length(label_ecg)
            icnt=icnt+1;
            data.chaninfo.unit{icnt} = 'uV';
            data.chaninfo.type{icnt} = 'physio';
        end
        
        data.cfg.dattype.sleep = []; % clear existing sleep
        
        begin_sp = NaN;
        stop_sp = NaN;
        begin_fl = NaN;
        stop_fl = NaN;
        SLP = [];
        inSleep = 0;
        inFlash = 0;
        
        SLEEP_list = {'asleep','sleep','pt_asleep','asleep_','drowsy'};
        AWAKE_list = {'awake','wake','eo','Eyes_open','waking_the_pt_up',...
            'mom_waking_child_','waking_pt_up','awake_eo','waking_up_','waking_up',...
            'I_think_this_patient_is_awake!','awake_','WAKES_UP_'};
        
        if ~isfield(data.cfg,'dattype')
            data.cfg.dattype = [];
        end
        
        % empty bad data field
        data.cfg.dattype.bad = [];
        
        % Mark initial 60 seconds bad
        if isfield(data.cfg.dattype,'bad')
            data.cfg.dattype.bad = [data.cfg.dattype.bad; 1 data.fsample*60];
            data.cfg.dattype.bad(data.cfg.dattype.bad == 0) = 1; % fix this so ICA doesn't crash
        else
            data.cfg.dattype.bad = [1 data.fsample*60];
        end
        
        % % Mark NaN sections bad
        
        where = find(abs(diff(isnan(data.trial{1}(1,:))))>0);
        if mod(length(where),2) == 1
            where(end+1) = length(data.trial{1});
        end
        
        data.cfg.dattype.bad = [data.cfg.dattype.bad; reshape(where,2,length(where)/2)'];
        
        
        for i = 1:length(data.cfg.event)
            if ~inSleep % if we haven't started a sleep block yet
                if ismember(upper(data.cfg.event(i).value),upper(SLEEP_list))
                    begin_sp = data.cfg.event(i).sample;
                    inSleep = 1;
                end
            elseif inSleep % if we're already in a sleep block
                if ismember(upper(data.cfg.event(i).value),upper(AWAKE_list)) % if sleep ends
                    stop_sp = data.cfg.event(i).sample;
                    SLP = [SLP; begin_sp stop_sp];
                    inSleep = 0;
                    % reset begin_sp and stop_sp
                    begin_sp = NaN; stop_sp = NaN;
                end
            end
            
            if strcmp('2.0_Hz',data.cfg.event(i).value)
                begin_fl = data.cfg.event(i).sample;
                inFlash = 1; % we have entered a light flash block
            elseif strcmp('20.0_Hz',data.cfg.event(i).value)
                inFlash = 0; % we are exiting the light flash block
            elseif strcmp('Off',data.cfg.event(i).value) && ~inFlash
                %^ the first Off marked after the 30 Hz stimulus
                stop_fl = data.cfg.event(i).sample;
            end
        end
        
        if ~isnan(begin_sp) && isnan(stop_sp) % if we started sleep block but didn't finish
            stop_sp = length(data.trial{1});
            SLP = [SLP; begin_sp stop_sp];
            display('Subject never woke up. Sleep block to end of recording.')
        elseif isempty(SLP)
            display('No sleep found!')
        end
        
        data.cfg.dattype.sleep = double(SLP);
        
        if ~isnan(begin_fl) && isnan(stop_fl) % if we started flash block but didn't finish
            fprintf('\nNo ''Off'' tag after the 20 Hz light stimulus! Skipping this one ...\n')
            continue
        elseif isnan(begin_fl)
            display('No light flashes found!')
        else
            % Add two seconds of wiggle room
            begin_fl = begin_fl - data.fsample*2;
            stop_fl  = stop_fl + data.fsample*2;
            % Append light flash block to bad data
            data.cfg.dattype.flash= double([begin_fl stop_fl]);
        end
        
        % Now that we've marked NaN sections as bad, replace nans with 0s,
        % otherwise scroll GUI will crash
        replace = isnan(data.trial{1});
        data.trial{1}(replace) = 0;
        
        outname = sprintf('%sAS_%s_%s.mat',DIRRESULT,ID,datestr(data.startdate(1),'yyyymmdd'));
        %if exist(outname,'file'), fprintf('Warning: File already exists!\n'), keyboard, end
        save(outname,'data','-v7.3')
    end
end
