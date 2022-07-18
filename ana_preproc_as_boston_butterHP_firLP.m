function [] = ana_preproc_as_boston_butter()

% Load and pre-process the AS data from Boston

addfieldtrip

load angelman_lay.mat lay elec

in                 = load('./data_raw/dataset_AS.mat');
idx                = find(strcmp(in.site,'Bos'));
in.ID              = in.ID(idx);
in.age             = in.age(idx);
in.age_sz          = in.age_sz(idx);
in.birthday        = in.birthday(idx);
in.del_start_coord = in.del_start_coord(idx);
in.del_stop_coord  = in.del_stop_coord(idx);
in.genotype        = in.genotype(idx);
in.sex             = in.sex(idx);
in.site            = in.site(idx);
in.date_EEG        = in.date_EEG(idx);
in.sz              = in.sz(idx);
in.test_file       = in.test_file(idx); % for Hannah 

for icnt = 1:length(in.ID)
    subject(icnt).ID=in.ID(icnt);
    subject(icnt).sex=in.sex(icnt);
    subject(icnt).age=in.age(icnt);
    subject(icnt).genotype = in.genotype(icnt);
    subject(icnt).test_file = in.test_file(icnt);
    subject(icnt).DC_offset = in.DC_offset(icnt);
    subject(icnt).site='BOS';
    subject(icnt).fname=sprintf('./data_raw/boston/%.3i_%s_*.edf',subject(icnt).ID,datestr(datenum(in.date_EEG{icnt},'DD.MM.YYYY'),'YYYYMMDD'));
end

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


DIRRESULT = './Monti/AS_butterHP_firLP/Boston_EDF/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end

%%

label_load = {'C3','C4','O1','O2','Cz','F3','F4','F7','F8','Fz','FP1','FP2','P3','P4','Pz','T7','T8','P7','P8'};
label_alt =  {'C3','C4','O1','O2','Cz','F3','F4','F7','F8','Fz','FP1','FP2','P3','P4','Pz','T3','T4','T5','T6'};
label_physio = {'ECGL','ECGR','LUE','RAE'};
label_physio_alt = {'ECGL','ECGR','LOC','ROC'};

for isub=1:length(subject)
    filename_all = subject(isub).fname;
    fname = dir(filename_all);
    DIRDATA = subject(isub).fname(1:max(strfind(subject(isub).fname,'/')));
    if isempty(fname)
        % try replacing the underscore with a space
        filename_all(25) = ' ';
        fname = dir(filename_all);
    end
    if ~isempty(fname) 
        clear cut datnum
        for ifile=1:length(fname)
            ID = fname(ifile).name(1:6);
            filename = fname(ifile).name;
            [~,hdr]=lab_read_edf_jfh([DIRDATA,filename]); % allows reading the event information and other header information propperly
            if ifile==1 % if we're taking the first chunk
                fsample=hdr.samplingrate;
                EVENTS=hdr.events;
                % create empty cell arrays
                clear event 
                event.duration = cell(1,length(EVENTS.DUR));
                event.value    = cell(1,length(EVENTS.TYP));
                event.sample   = cell(1,length(EVENTS.POS));
                for ievt = 1:length(EVENTS.TYP)
                    event(ievt).duration = EVENTS.DUR(ievt);
                    event(ievt).value    = EVENTS.TYP{ievt};
                    event(ievt).sample   = EVENTS.POS(ievt);
                end
                % sanity check
                if ~strcmp(ID,hdr.subject.name(1:6))
                    isub
                    filename
                    keyboard
                    warning('check!')
                end
                [hdr2,data] = edfread([DIRDATA,filename]); % allows reading the data with correct scaling
                if subject(isub).DC_offset % check for DC offset
                    % correct DC offset
                    dat = dat - repmat(median(dat')',1,length(dat));
                end
                label_in = hdr2.label;
                datnum(ifile) = datenum([hdr2.startdate,'.',hdr2.starttime],'dd.mm.yy.HH.MM.SS');
                cut(ifile)=0;
            else % if we're not on the first data chunk for that subject
                if fsample~=hdr.samplingrate; error('samplingrate'), end
                cut(ifile)=size(data,2);
                EVENTS=hdr.events;
                % add to events to subfields
                len = length(event);
                for ievt = 1:length(EVENTS.TYP)
                    event(len+ievt).duration = EVENTS.DUR(ievt);
                    event(len+ievt).value    = EVENTS.TYP{ievt};
                    event(len+ievt).sample   = size(data,2)+EVENTS.POS(ievt);
                    %%#JOERG: Don't we have to add the sample number from the previous files!?
                    %event(len+ievt).sample   = EVENTS.POS(ievt);
                end
                sex=hdr.subject.sex;
                if ~strcmp(sex,hdr.subject.sex); error('sex'), end
                [hdr2,data_tmp] = edfread([DIRDATA,filename]);      % allows reading the data with correct scaling
                data=cat(2,data,data_tmp);
                datnum(ifile) = datenum([hdr2.startdate,'.',hdr2.starttime],'dd.mm.yy.HH.MM.SS');
            end
        end
        
        % reorder the EEG channels if neccessary
        [~,where] = ismember(upper(label_in),upper(label_load));
        if sum(where>0)~=19
            [~,where] = ismember(upper(label_in),upper(label_alt));
        end
        assert(sum(where>0)==19,'missing channels')
        clear idx_chan % @Joel: make sure the channel order is correct!!! -> apply to other preproc scripts
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
        clear idx_chan 
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
        
        % Set FIR filter length arround the data fusion points to nan
        for icnt=2:length(cut)
            idx_start = max(cut(icnt)-fsample*forder,1);
            idx_end = min(cut(icnt)+fsample*forder,size(dat,2));
            idx = idx_start:idx_end;
            dat(:,idx)=nan;
            dat_gmm(:,idx)=nan;
            gamma(:,idx)=nan;
            eogH(idx)=nan;
            eogV(idx)=nan;
            ecg1(idx)=nan;
            ecg2(idx)=nan;
        end
        
        % Put into Fieldtrip structure
        clear data
        data.cfg.event    = event;
        data.cfg.age      = subject(isub).age;
        data.cfg.sex      = subject(isub).sex;
        data.ID           = sprintf('%i',subject(isub).ID);
        data.elec         = elec;
        %data.dat_gmm       = single(dat_gmm);
        data.dimord       = 'chan_time';
        data.startdate    = datnum;
        data.fsample      = fsample;
        data.label        = label_alt;
        data.time{1}      = time;
        data.label{end+1} = 'EOGV';
        data.label{end+1} = 'EOGH';
        data.label{end+1} = 'ECG1';
        data.label{end+1} = 'ECG2';
        data.label{end+1} = 'GAMMA';
        data.orgForm      = 'EDF';
        data.time{1}      = single(time);
        data.trial{1}     = single(dat);
        data.trial{1}(end+1,:) = single(eogV);
        data.trial{1}(end+1,:) = single(eogH);
        data.trial{1}(end+1,:) = single(ecg1);
        data.trial{1}(end+1,:) = single(ecg2);
        % we've already taken the wave envelope so we can average without
        % cancelation
        data.trial{1}(end+1,:) = single(max(gamma));
        data.sampleinfo(1,:)   = [1,size(dat,2)];
        data.cfg.trl.trl(1,:)  = [1,size(dat,2),0];
        data.site = 'BOS';
        
         % add units and labels for auxillary channels 
        data.cfg.raw.ctypelabels      = {'EEG','physio'};
        data.cfg.raw.ctypemainlabel   = 'EEG';
        data.cfg.raw.ctypephysiolabel = 'physio';
        for icnt=1:length(label)
            data.chaninfo.unit{icnt} = 'uV';
            data.chaninfo.type{icnt} = 'EEG';
        end

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
  
        %outname = sprintf('%sAS_%.3i_%dm.mat',DIRRESULT,subject(isub).ID,subject(isub).age);
        outname = sprintf('%sAS_%.3i_%s.mat',DIRRESULT,subject(isub).ID,datestr(data.startdate(1),'yyyymmdd'));
        save(outname,'data','-v7.3')
    end
end
