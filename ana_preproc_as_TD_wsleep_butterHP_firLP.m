function [] = ana_preproc_as_TD_wsleep_2021()

% Load and pre-process the awake TD data from Catherine Chu
% Use this one

dbstop if error

% preprocessing parameter
preproc_parameter;

T = readtable('./TD_controls/MGH_EEG_code_abbrev.xlsx');

DIRDATA_W  = './TD_controls/Raw/Awake_EEG/';
DIRDATA_N2 = './TD_controls/Raw/Asleep_EEG/';
DIRRESULT  = './TD_controls/Imported/Detrend_butterHP_firLP/';  if ~exist(DIRRESULT), mkdir(DIRRESULT), end

startndx = 1;

if startndx ~= 1
    fprintf('\nWarning, startndx set to %i\n',startndx)
end

for ifile = startndx:size(T,1)
    % check if there's a corresponding sleep file
    if strcmp(T.Linked_to_patient__{ifile},'')
        fprintf('\nNo corresponding sleep file for %s, skipping this one ...\n',T.EDFCODE{ifile})
        continue
    else
        wakefile  = T.EDFCODE{ifile};
        sleepfile = T.Linked_to_patient__{ifile};
        SID = sscanf(T.Linked_to_patient__{ifile},'%f_N2'); % use the ID from the sleep file
        alt_ID = sscanf(T.EDFCODE{ifile},'%f_W');
    end
    
    ageout =  sscanf(T.Age{ifile},'%i y %i m');
    try
        % if there is data for year and month
        ageout =  sscanf(T.Age{ifile},'%i y %i m');
        assert(length(ageout)==2)
        age = ageout(1)*12 + ageout(2); % age in months
    catch
        % if there is only data for year (no month)
        assert(length(ageout)==1)
        age = ageout(1)*12;
    end
    
    % awake processed data
    data  = create_FT_struct(wakefile,DIRDATA_W,SID);
    data.age = age;
    data.comments = 'First trial is awake, second trial is asleep';
    
    % asleep processed data
    data2  = create_FT_struct(sleepfile,DIRDATA_N2,SID);
    data.trial{2} = data2.trial{1};
    
    % add time from sleep trial
    data.time{2} = data2.time{1};
    
    % add sleep to sample and trial fields
    data.sampleinfo(2,:) = [(1 + size(data.trial{1},2)),(size(data.trial{1},2) ...
        + size(data.trial{2},2))];
    data.cfg.trl.trl(2,:) = [(1 + size(data.trial{1},2)),(size(data.trial{1},2) ...
        + size(data.trial{2},2)), data.sampleinfo(2,1)];
    
    % mark all of second trial as sleep
    data.cfg.dattype.sleep = [data.sampleinfo(2,1) data.sampleinfo(2,2)];
    data.cfg.dattype.SLEEP_MONTI = [data.sampleinfo(2,1) data.sampleinfo(2,2)];
    
    % append IDs
    data.SID = SID;
    data.alt_ID = alt_ID;
    
    % Get data markings and ICA from previously cleaned files
    
    try
        olddat = load(sprintf('./TD_controls/Imported/Detrend_butter/Clean/%i_TD_%im.mat',SID,age));
        data.cfg.dattype = olddat.data.cfg.dattype; % get data annotations
        data.cfg.ica = olddat.data.cfg.ica; % get ICA
        data.cfg.raw.ctype.bad = olddat.data.cfg.raw.ctype.bad; % get bad chans
        data.cfg.raw.ctype.EEG = olddat.data.cfg.raw.ctype.EEG; % get EEG chans
    catch
        fprintf('\nSkipping this one, no clean data (we skipped it before, must be bad) ...\n')
        continue
    end
    
    % Account for filter edge artifacts (mark bad)
    edgeart = data.fsample*(1/f_hp)*3; % duration of edge artifact
    data.cfg.dattype.bad = [data.cfg.dattype.bad; 1 edgeart]; % beginning of trial 1
    data.cfg.dattype.bad = [data.cfg.dattype.bad; length(data.trial{1}) - edgeart length(data.trial{1})]; % end of trial 1
    data.cfg.dattype.bad = [data.cfg.dattype.bad; length(data.trial{1})+1 length(data.trial{1})+ edgeart]; % beginning of trial 2
    data.cfg.dattype.bad = [data.cfg.dattype.bad; length(data.trial{1}) + length(data.trial{2}) - edgeart length(data.trial{1}) + length(data.trial{2})]; % end of trial 2
    
%     % Also mark bad sections with disconnect amp
%     
%     % First trial (wake)
%     a=find(diff(mode(abs(diff(data.trial{1}))<0.01))>0)';
%     b=find(diff(mode(abs(diff(data.trial{1}))<0.01))<0)';
%     if length(a) ~= length(b)
%         c2 = sort([a; b]); % sort all jump points
%         c3 = c2(find(diff(round(c2,-2))~=0));
%         if mod(length(c3),2)==1
%             c3 = c3(1:end-1); % just discard the last point if length is odd
%         end
%         
%         c = reshape(c3,length(c3)/2,2);
%         a = c(:,1); b = c(:,2);
%     end
%     for irow = 1:size(a,1)
%         if b(irow) - a(irow) >= data.fsample
%             data.cfg.dattype.bad = [data.cfg.dattype.bad; ...
%                 a(irow)-data.fsample b(irow)+data.fsample];
%         end
%     end
%     
%     clear a b c % just to be safe
%     
%     % Second trial (sleep)
%     a=find(diff(mode(abs(diff(data.trial{2}))<0.01))>0)' + length(data.trial{1});
%     b=find(diff(mode(abs(diff(data.trial{2}))<0.01))<0)' + length(data.trial{1});
%     if length(a) ~= length(b)
%         c2 = sort([a; b]); % sort all jump points
%         c3 = c2(find(diff(round(c2,-2))~=0));
%         if mod(length(c3),2)==1
%             c3 = c3(1:end-1); % just discard the last point if length is odd
%         end
%         
%         c = reshape(c3,length(c3)/2,2);
%         a = c(:,1) + length(data.trial{1}); 
%         b = c(:,2) + length(data.trial{1});
%     end
%     for irow = 1:size(a,1)
%         if b(irow) - a(irow) >= data.fsample
%             data.cfg.dattype.bad = [data.cfg.dattype.bad; ...
%                 a(irow)-data.fsample b(irow)+data.fsample];
%         end
%     end
%     
%     clear a b c % just to be safe
    
    % save data
    outname = sprintf('%s%i_TD_%im.mat',DIRRESULT,SID,age);
    fprintf('\nSaving file %s\n',outname)
    save(outname,'data','-v7.3')
    clear data
end

end

function[data] = create_FT_struct(filename,path,SID)

%% Channels labels

label_load = {'C3','C4','O1','O2','Cz','F3','F4','F7','F8','Fz','FP1','FP2','P3','P4','Pz','T7','T8','P7','P8'};
label_alt  = {'C3','C4','O1','O2','Cz','F3','F4','F7','F8','Fz','FP1','FP2','P3','P4','Pz','T3','T4','T5','T6'};
label_physio = {'EKG','EOG'};
label_physio_alt = {'ECG','ECG','LOC','ROC'};
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

sex = filename(end);
filename = sprintf('%s.edf',filename);

% Do the awake file
fprintf('\nReading %s\n',sprintf('%s%s',path,filename))

[~,hdr]=lab_read_edf_jfh(sprintf('%s%s',path,filename)); % allows reading the event information and other header information properly
fsample=hdr.samplingrate; % use this one just to get the correct sampling rate

[hdr2,data] = edfread(filename); % allows reading the data with correct scaling
label_in = hdr2.label;
%datnum = datenum([hdr2.startdate,'.',hdr2.starttime],'dd.mm.yy.HH.MM.SS');

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

% Comment out line below if you don't want to do the detrending
% Detrend the data to correct for DC offsets
% Must use mean (median doesn't remove them very well)
dat = single(dat) - repmat(mean(double(dat)')',1,length(dat));

% not sure it has to be done in a loop, but detrend may prefer vectors over
% matricies

for ich = 1:size(dat,1)
    dat(ich,:) = detrend(dat(ich,:));
end

% Load EOG
clear physio ecg1 ecg2 eog1 eog2
%physio = nan(4,length(data.trial{1}); % memory allocation
[where] = contains(upper(label_in),upper(label_physio));
if sum(where>0)~=4
    [where] = contains(upper(label_in),upper(label_physio_alt));
end

%assert(sum(where>0)==4,'physio channels not found??')

clear idx_chan % @Joel: make sure the channel order is correct!!! -> apply to other preproc scripts
% for icnt=1:length(label_physio)
%     idx_chan(icnt) = find(where==icnt);
% end

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


%%

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

dat_gmm = dat_hp - dat;
dat_gmm = single(filtfilt(B_line{icnt}, A_line{icnt}, double(dat_gmm)')'); % only notch 60 Hz to avoid numerical problems

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
dat_gmm = dat_gmm - repmat(mean(dat_gmm,1),[size(dat_gmm,1),1]); % avg reference

% get the gamma-band envelope (helps manual artifact rejection)
gamma = nan(size(dat_gmm));
if mod(fnyq,2) == 0 % if even
    n = fnyq + 1; % make odd number
else
    warning('Sampling rate is odd number!')
    fprintf('\nDetected sampling rate = %i Hz\n',fsample)
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
data.ID           = SID;
% NO EVENTS IN THESE DATA??
% data.cfg.event    = event;
data.elec         = elec;
%data.dat_gmm       = single(dat_gmm);
data.dimord       = 'chan_time';
%data.startdate    = datnum;
data.fsample      = fsample;
% use these labels so they match the lay file!
data.label        = elec.label';
data.time{1}      = time;

if exist('eog1','var')
    data.label{end+1} = 'eog1';
    data.trial{1}(end+1,:) = single(eog1);
end
if exist('eog2','var')
    data.label{end+1} = 'eog2';
    data.trial{1}(end+1,:) = singel(eog2);
end
if exist('ECG1','var')
    data.label{end+1} = 'ECG1';
    data.trial{1}(end+1,:) = single(ECG1);
end
if exist('ECG2','var')
    data.label{end+1} = 'ECG2';
    data.trial{1}(end+1,:) = single(ECG2);
end

data.label{end+1} = 'GAMMA';
data.orgForm      = 'EDF';
data.orgName      = filename;
data.sex          = sex;
data.f_hp         = f_hp;
data.l_lp         = f_lp;
data.time{1}      = single(time);
data.trial{1}     = single(dat);
% we've already taken the wave envelope so we can average without
% cancelation
data.trial{1}(end+1,:) = single(max(gamma));
data.sampleinfo(1,:) = [1,size(dat,2)];
data.cfg.trl.trl(1,:) = [1,size(dat,2),data.sampleinfo(1,1)];
data.site = 'BOS';

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

end




