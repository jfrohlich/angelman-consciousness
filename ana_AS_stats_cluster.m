%function[wk_avg,sp_avg,pvals] = ana_AS_MSE_stats(scales)

% input: scales are the time scales we want to average
clearvars
rng(345800);  % For reproducibility (jack-kniefed correlations)

fsample = 200; % all files should be downsampled to this freq before computing complexity measures
scales = 1:20;
nyqfoi = (fsample/2)./scales; % MSE
strict = true; % do we filter out clinical EEG?
questionable = {'AS_100007_20180706','AS_104972_20170515'};

close all
dbstop if error
viewts = false;
match = true; % look at the data with 4 Hz low pass filter?
Nperm = 10^4; % number of permutations
surrtest = false; % surrogate data testing? 
load 1020_labels
if ~match
    load select_files
end

% MSE preferences
MSEmethod = 'Xie';
tolerance = 0.15;

switch strict
    case true
        switch match
            case true
                DIRFIGURE = './figures/match/';
            case false
                DIRFIGURE = './figures/';
        end
        
    case false
        switch match
            case true
                DIRFIGURE = './figures/match/withclinical/';
            case false
                DIRFIGURE = './figures/withclinical/';
        end
end

if ~exist(DIRFIGURE,'dir'), mkdir(DIRFIGURE), end

switch match
    case true
        pth = './MSE_March/match/Xie/r=0.15/dynr/';
        wake = dir(sprintf('%swake/*.mat',pth));
        sleep = dir(sprintf('%ssleep/*.mat',pth));
        frqpth = './freq_out/March/nan=0.2/match/';
        frqwake  = dir(sprintf('%swake/*.mat',frqpth));
        frqsleep = dir(sprintf('%ssleep/*.mat',frqpth)); 
        pthLZC = './gLZC_March/match/';
    case false
        pth = './MSE_March/Xie/r=0.15/dynr/';
        wake = dir(sprintf('%swake/*.mat',pth));
        sleep = dir(sprintf('%ssleep/*.mat',pth));
        frqpth = './freq_out/March/nan=0.2/';
        frqwake  = dir(sprintf('%swake/*.mat',frqpth));
        frqsleep = dir(sprintf('%ssleep/*.mat',frqpth)); 
        pthLZC = './gLZC_March/';
end
pthTD = './MSE_TD/dynr/';
wakeTD= dir(sprintf('%s/*.mat',pthTD));

mse_wk = cell(1,1);
mse_sp = cell(1,1);
win_wk = cell(1,1);
win_sp = cell(1,1);
dur_wk_all = cell(1,1);
dur_sp_all = cell(1,1);
SID_wk = cell(1,1);
SID_sp = cell(1,1);

lzc_wk = cell(1,1);
lzc_sp = cell(1,1);
win_wk_lzc = cell(1,1);
win_sp_lzc = cell(1,1);
dur_wk_all_lzc = cell(1,1);
dur_sp_all_lzc = cell(1,1);
SID_wk_lzc = cell(1,1);
SID_sp_lzc = cell(1,1);
trg_wk_all = cell(1,1);
trg_sp_all = cell(1,1);
trg_wk_segs = cell(1,1);
trg_sp_segs = cell(1,1);

IDs = [];
allage = cell(1,1);
allgender = [];
allgenotype = [];
allcog = cell(1,1);
alllng = cell(1,1);
allVCog = cell(1,1);
allNVCog = cell(1,1);
all_pow_sp = cell(1,1);
all_pow_wk = cell(1,1);
wk_pow = cell(1,1);
sp_pow = cell(1,1);
wk_rel = cell(1,1);
sp_rel = cell(1,1);
all_sec_wk = cell(1,1);
all_sec_sp = cell(1,1);
wk_sec = [];
sp_sec = [];
dur_wk = [];
dur_sp = [];
cog_score = [];
lng_score = [];



T = readtable('AS_table_2019_JF.csv','Delimiter',','); % exported from the master spreadsheet
B = readtable('AS_Bayley.xlsx');

% %% load TD data
% 
% mse_td = cell(1,1);
% win_td = [];
% dur_td_all = [];
% SID_td = [];
% age_td = [];
% td_scl = [];
% td_avg = [];
% td_all = [];
% 
% for ifile = 1:length(wakeTD)
%     
%     % get SID and age
%     underscore = strfind(wakeTD(ifile).name,'_');
%     stop = underscore(3)-1; % everything up to the third underscore 
%     fstr = wakeTD(ifile).name(1:stop);
%     SID = fstr(underscore(1)+1:underscore(2)-1);
%     age = fstr(underscore(2)+1:strfind(fstr,'m')-1);
%     
%     % load empirical values
%     We_TD = load(sprintf('%s/%s',pthTD,wakeTD(ifile).name));
%     fprintf('\nSuccessfully loaded MSE from TD %s\n',SID)
%  
%     % load surrogate values
%     try
%         Ws_TD = load(sprintf('%s/Surrogate/%s',pthTD,sprintf('%s%s',...
%             wakeTD(ifile).name(1:end-4),'_SURROGATE')));
%         fprintf('\nSuccessfully loaded surrogate MSE from TD %s\n',SID)
%     catch
%     %    Ws_TD = load(sprintf('%s/Surrogate/%s',pth,sprintf('%s%s',...
%     %        wakeTD(ifile).name(1:end-37),'_wake_MSE_SURROGATE')));    
%         fprintf('\n Surrogate not found, skipping this one ... \n')
%         continue
%     end
%     
%     
%      % find "buried" data structures (these got buried due to a glitch in
%     % the script use gather() to turn gpuArray back to regular array
%    
%     % recursively find the real structure
%     while isfield(We_TD.MSEout,'MSEout')
%         We_TD.MSEout = We_TD.MSEout.MSEout;
%     end
%     while isfield(Ws_TD.MSEout,'MSEout')
%         Ws_TD.MSEout = Ws_TD.MSEout.MSEout;
%     end
%     
%     if isfield(We_TD.MSEout,'surrogate')
%         assert(~We_TD.MSEout.surrogate,'This empirical file is actually a surrogate')
%     end
%     if isfield(Ws_TD.MSEout,'surrogate')
%         assert(Ws_TD.MSEout.surrogate,'This surrogate file is actually empirical data')
%     end
% 
%     assert(Ws_TD.MSEout.cfg.tolerance == We_TD.MSEout.cfg.tolerance,'Different tolerances set for empirical and surrogate data')
%     
%     % Build new structure by combining info from empirical and surrogate
%     % data
%     
%     W_TD.mse = We_TD.MSEout.mse-Ws_TD.MSEout.mse; % take difference
%     %W_TD.mse = We_TD.MSEout.mse; 
%     W_TD.cfg = We_TD.MSEout.cfg;
%     W_TD.scale = We_TD.MSEout.scale;
%     W_TD.n_valid = min([We_TD.MSEout.n_valid; Ws_TD.MSEout.n_valid]);
%     W_TD.prc_used = We_TD.MSEout.prc_used;
%     W_TD.dur_used = We_TD.MSEout.dur_used;
%     
%     min_valid_mse = 100; % Grandy et al., 2016
%     wk_gdwin = W_TD.n_valid >= min_valid_mse;
%     if sum(wk_gdwin) == 0 % at least one good window!
%         fprintf('Skipping this file, no good windows: %s\n', wakeTD(ifile).name)
%         continue
%     end
%     
%     wk_win_avg = nanmean(W_TD.mse(:,:,wk_gdwin),3);
%     mse_td{ifile} = wk_win_avg;
%     win_td(ifile) = sum(W_TD.n_valid(wk_gdwin));
%     dur_td_all(ifile) = W_TD.dur_used;
%     SID_td(ifile) = str2num(SID);
%     age_td(ifile) = str2num(age);
%     
%     td_avg(:,ifile) = mean(mse_td{ifile}(:,scales),2);
%     td_scl(:,ifile) = mean(mse_td{ifile}(:,scales),1);
%     td_all      = cat(3,td_all,mse_td{ifile});
%     clear SID fstr
% end

%% load Angelman data

start_idx = 1;

if start_idx ~= 1
    fprintf('\n Warning, start_idx is set to %i',start_idx)
end

for ifile = start_idx:length(wake)
    
    % Make sure we only load files that were used in the matched analysis
    
    underscore = strfind(wake(ifile).name,'_');
    stop = underscore(3)-1; % everything up to the third underscore 
    fstr = wake(ifile).name(1:stop);
    if ~match && ~ismember(fstr,select_files)
        continue
    else
        fprintf('\n    Extracting output for %s\n',fstr)
    end
    
    pass = true;
    
    %% Check for files we're not sure about (clinical EEGs, files with sedation, et.c)
    for i = 1:length(questionable)
        if contains(wake(ifile).name,questionable{i})
            switch strict
                case true
                    pass = false; % because we can't have a continue statement here
                case false
                    fprintf('\n Allowing questionable file \n')
            end
        end
    end
    
    if ~pass
        fprintf('\n Questionable file, skipping this one\n')
        continue
    end
    %%
    
    fprintf('\n Extracting power and MSE values from %s \n',wake(ifile).name(1:end-14))
    hasAge = ~isempty(find(ismember(wake(ifile).name(1:15),'m')));
    
    switch hasAge
        case true
            SID = str2num(wake(ifile).name(4:9));
            if strcmp(wake(ifile).name(14),'m')
                age =  str2num(wake(ifile).name(11:13));
            elseif strcmp(wake(ifile).name(13),'m')
                age =  str2num(wake(ifile).name(11:12));
            else
                error('Age not coded in file name')
            end
            date = []; % these files do not have date info in file name
        case false
            hyph = ~isempty(find(ismember(wake(ifile).name(1:15),'-')));
            switch hyph
                case true
                    SID = str2num(wake(ifile).name(1:6));
                    age = []; % these files do not have age info in file name
                    year = str2num(wake(ifile).name(8:11));
                    month = str2num(wake(ifile).name(13:14));
                    day = str2num(wake(ifile).name(16:17));
                    date = day + month * 10^2 + year * 10^4;
                case false
                    SID = str2num(wake(ifile).name(4:9));
                    age = []; % these files do not have age info in file name
                    year = str2num(wake(ifile).name(11:14));
                    month = str2num(wake(ifile).name(15:16));
                    day = str2num(wake(ifile).name(17:18));
                    date = day + month * 10^2 + year * 10^4;
            end
    end
    
    if ~isempty(date)
        isubj= find(T.TimeRecording == date & T.ID == SID);
        age = T.AgeRecording(isubj);
    elseif ~isempty(age)
        isubj = find(T.AgeRecording == age & T.ID == SID);
        date = T.TimeRecording(isubj);
    else
        error(fprintf('Both age and date are missing for %s\n',wake(ifile).name))
    end
    
    assert(length(isubj)==1,'Error looking up dataset in table ...')
    
    % Check in table if subject has useable data
    if ~T.Use(isubj)
        fprintf('Skipping bad subject %s\n',num2str(SID))
        continue
    else
        if isempty(age)
            age = T.AgeRecording(isubj);
        end
    end
    
    % skip data from adult subjects
    if age > 18*12
        fprintf('\n Adult subject, skipping data ... \n')
        continue
    end
    
    BCog = B.bayley_cognitive_raw(B.ID == SID & round(B.Age) == age);
    if isempty(BCog) % if empty, take the closest age
        [~,idx] = min(abs(B.Age(B.ID==SID) - age));
        BCogall = B.bayley_cognitive_raw(B.ID == SID);
        BCog = BCogall(idx);
    end
    EL = B.bayley_expressive_raw(B.ID == SID & round(B.Age) == age);
    if isempty(EL) % if empty, take the closest age
        [~,idx] = min(abs(B.Age(B.ID==SID) - age));
        ELall = B.bayley_expressive_raw(B.ID == SID);
        EL = ELall(idx);
    end
    RL = B.bayley_receptive_raw(B.ID == SID & round(B.Age) == age);
    if isempty(RL) % if empty, take the closest age
        [~,idx] = min(abs(B.Age(B.ID==SID) - age));
        RLall = B.bayley_receptive_raw(B.ID == SID);
        RL = RLall(idx);
    end
    BLang = mean([EL; RL]);
    VCog  = T.Verbal(T.ID == SID & T.AgeRecording == age);
    NVCog = T.Nonverbal(T.ID == SID & T.AgeRecording == age);
    
    if strcmpi(T.Gender(isubj),'M')
        gender = 1;
    else
        gender = 0;
    end
    
    if strcmpi(T.Deletion(isubj),'Yes')
        genotype = 1;
    else
        genotype = 0;
    end
    
    %% Load power (do this first so we can skip EEGs without valid power output)    
    
    % minimum number of good freq transform windows for lowest frequency 
    n_thresh = 30;
    
    if ~ismember(SID,IDs)
        IDs = [IDs; SID];
    end
    idx = find(IDs == SID);
    
    % SLEEP POWER
    for jfile = 1:length(frqsleep)
        datestr = num2str(date);
        altdatestr = sprintf('%s-%s-%s',datestr(1:4),datestr(5:6),datestr(7:8));
        if contains(frqsleep(jfile).name,num2str(SID)) && ... % find the power output from the same file
                ( contains(frqsleep(jfile).name,datestr) ||  ...
                contains(frqsleep(jfile).name,altdatestr) || ...
                contains(frqsleep(jfile).name,sprintf('%im',age)) )
            fprintf('     Loading power...\n')
            load(sprintf('%s%s',frqpth,'sleep/',frqsleep(jfile).name),'pow','pow_full','pow_var','foi','n','cfg')
            if exist('cfg','var')
                assert(cfg.allow_fraction_nan==0.2,'Wrong fraction NaNs')
            end
            break
        end
    end
    
    % load power vars into structure
    if exist('pow','var') && n(1) >= n_thresh % if there's power and at least n good window for lowest frequency 
        sp_freq.pow = pow;
        sp_freq.pow_full = pow_full;
        sp_freq.pow_var = pow_var;
        sp_freq.foi = foi;
        sp_freq.n = n;
        % relative power 
        pow_rel = nan(size(pow));
        for ich = 1:size(pow,1)
            TOT = trapz(log2(foi),pow(ich,:));
            pow_rel(ich,:) = pow(ich,:)./TOT; % devide by total power to normalize 
        end
        sp_freq.pow_rel = pow_rel;
        
        clear pow pow_full pow_rel freq cfg
    else
        fprintf('\n Bad sleep power, only %i good windows, skipping this one ...\n',n(1))
        continue
    end
    
    % AWAKE POWER
    for jfile = 1:length(frqwake)
        datestr = num2str(date);
        altdatestr = sprintf('%s-%s-%s',datestr(1:4),datestr(5:6),datestr(7:8));
        if contains(frqwake(jfile).name,num2str(SID)) && ... % find the power output from the same file
                ( contains(frqwake(jfile).name,datestr) ||  ...
                contains(frqwake(jfile).name,altdatestr) || ...
                contains(frqwake(jfile).name,sprintf('%im',age)) )
            fprintf('     Loading power...\n')
            load(sprintf('%s%s',frqpth,'wake/',frqwake(jfile).name),'pow','pow_full','pow_var','foi','n','cfg')
            if exist('cfg','var')
                assert(cfg.allow_fraction_nan==0.2,'Wrong fraction NaNs')
            end                
            break
        end
    end
    
    if exist('pow','var') && n(1) >= n_thresh % is there's power and at least n good window for lowest frequency 
        % load power vars into structure
        wk_freq.pow = pow;
        wk_freq.pow_full = pow_full;
        wk_freq.pow_var = pow_var;
        wk_freq.foi = foi;
        wk_freq.n = n;
        % relative power 
        pow_rel = nan(size(pow));
        for ich = 1:size(pow,1)
            TOT = trapz(log2(foi),pow(ich,:));
            pow_rel(ich,:) = pow(ich,:)./TOT; % devide by total power to normalize 
        end
        wk_freq.pow_rel = pow_rel;

        clear pow pow_full pow_rel n cfg
    else
        fprintf('\n Bad awake power, only %i good windows, skipping this one ...\n',n(1))
        continue
    end
    
    % store sleep power under the cell index for that subject
    if idx > 1 && idx <= length(all_pow_sp)
        all_pow_sp{idx}(length(all_pow_sp{idx})+1) = sp_freq;
        %all_sec_sp{idx}(length(all_sec_sp{idx})+1) = Nsecsp;
    else
        all_pow_sp{idx} = sp_freq;
        %all_sec_sp{idx} = Nsecsp;
    end
    
    % store wake power under the cell index for that subject
    if idx > 1 && idx <= length(all_pow_wk)
        all_pow_wk{idx}(length(all_pow_wk{idx})+1) = wk_freq;
        %all_sec_wk{idx}(length(all_sec_wk{idx})+1) = Nsecwk;
    else
        all_pow_wk{idx} = wk_freq;
        %all_sec_wk{idx} = Nsecwk;
    end

    fprintf('\nSuccessfully pulled sleep and wake from this subject: %s\n',wake(ifile).name(1:end-12))
    
    clear sp_freq wk_freq 
    
    %% Load MSE
    
    % load empirical values
    We = load(sprintf('%swake/%s%s',pth,fstr,'_wake_MSE'));
    
    % find "buried" data structures (these got buried due to a glitch in
    % the script use gather() to turn gpuArray back to regular array
   
    % recursively find the real structure
    while isfield(We.MSEout,'MSEout')
        We.MSEout = We.MSEout.MSEout;
    end
    
    % check fields
    assert(We.MSEout.cfg.tolerance == tolerance,'Wrong tolerance')
    assert(strcmp(We.MSEout.cfg.type,MSEmethod),'Wrong method')
    assert(We.MSEout.cfg.dynr,'This file does not have dynamic tolerance')
    assert(We.MSEout.cfg.new_srate == fsample,'File wasn''t downsampled?') 
    assert(~We.MSEout.surrogate,'This empirical file is actually a surrogate')
    
    if surrtest
        % load surrogate values
        try
            Ws = load(sprintf('%swake/Surrogate/%s',pth,sprintf('%s%s',...
                fstr,'_wake_MSE_SURROGATE')));
        catch
            fprintf('\n Surrogate not found, skipping this one ... \n')
            continue
        end

        while isfield(Ws.MSEout,'MSEout')
            Ws.MSEout = Ws.MSEout.MSEout;
        end
        
        % check fields
        assert(Ws.MSEout.cfg.tolerance == tolerance,'Wrong tolerance')
        assert(strcmp(Ws.MSEout.cfg.type,MSEmethod),'Wrong method')
        assert(Ws.MSEout.cfg.dynr,'This file does not have dynamic tolerance')
        assert(Ws.MSEout.cfg.new_srate == fsample,'File wasn''t downsampled?') 
        assert(Ws.MSEout.surrogate,'This surrogate file is actually empirical data')
        
    end
    
    
    % Build new structure by combining info from empirical and surrogate
    % data
    
    switch surrtest
        case true
            W.mse = We.MSEout.mse-Ws.MSEout.mse; % take difference
            W.n_valid = min([We.MSEout.n_valid; Ws.MSEout.n_valid]);
        case false
            W.mse = We.MSEout.mse; 
            W.n_valid = We.MSEout.n_valid; 
    end
    
    W.cfg = We.MSEout.cfg;
    W.scale = We.MSEout.scale;
    W.prc_used = We.MSEout.prc_used;
    W.dur_used = We.MSEout.dur_used;
    if match
        % targeted comparison
        W.wk_len = We.MSEout.cfg.wk_len;
        W.wk_segs_used = We.MSEout.cfg.wk_segs_used;
    end
    
    fprintf('\nLoading awake file %s \n',wake(ifile).name)
    try % attempt to load the corresponding sleep file
        sleepfile = strcat(fstr,'_sleep_MSE');
        
        % load empirical values
        Se = load(sprintf('%ssleep/%s.mat',pth,sleepfile));
        
        % find "buried" data structures (these got buried due to a glitch in
        % the script use gather() to turn gpuArray back to regular array

        % recursively find the real structure
        while isfield(Se.MSEout,'MSEout')
            Se.MSEout = Se.MSEout.MSEout;
        end

        catch
            fprintf('\n No sleep data for this subject ....\n')
            continue
    end
    
    % check fields
    assert(Se.MSEout.cfg.tolerance == tolerance,'Wrong tolerance')
    assert(strcmp(Se.MSEout.cfg.type,MSEmethod),'Wrong method')
    assert(Se.MSEout.cfg.dynr,'This file does not have dynamic tolerance')
    assert(Se.MSEout.cfg.new_srate == fsample,'File wasn''t downsampled?') 
    assert(~Se.MSEout.surrogate,'This empirical file is actually a surrogate')
    
    if surrtest 
        try
            % load surrogate values
            Ss = load(sprintf('%ssleep/Surrogate/%s.mat',pth,sprintf('%s%s',...
                sleepfile,'_SURROGATE')));
        catch
%             Ss = load(sprintf('%ssleep/Surrogate/%s.mat',pth,sprintf('%s%s',...
%                 sleepfile(1:end-34),'_sleep_MSE_SURROGATE')));
            fprintf('\n Surrogate not found, skipping this one ... \n')
            continue
        end
        while isfield(Ss.MSEout,'MSEout')
            Ss.MSEout = Ss.MSEout.MSEout;
        end
        % check fields
        assert(Ss.MSEout.cfg.tolerance == tolerance,'Wrong tolerance')
        assert(strcmp(Ss.MSEout.cfg.type,MSEmethod),'Wrong method')
        assert(Ss.MSEout.cfg.dynr,'This file does not have dynamic tolerance')
        assert(Ss.MSEout.cfg.new_srate == fsample,'File wasn''t downsampled?') 
        assert(Ss.MSEout.surrogate,'This empirical file is actually empirical data')
    end
    
    % Build new structure by combining info from empirical and surrogate
    % data

    switch surrtest
        case true
            S.mse = Se.MSEout.mse-Ss.MSEout.mse; % take difference
            S.n_valid = min([Se.MSEout.n_valid; Ss.MSEout.n_valid]);
        case false
            S.mse = Se.MSEout.mse;
            S.n_valid = Se.MSEout.n_valid; 
    end
    
    S.cfg = Se.MSEout.cfg;
    S.scale = Se.MSEout.scale;
    S.prc_used = Se.MSEout.prc_used;
    S.dur_used = Se.MSEout.dur_used;
    if match
        % targeted comparison
        S.sp_len = Se.MSEout.cfg.sp_len;
        S.sp_segs_used = Se.MSEout.cfg.sp_segs_used;
    end
    
    fprintf('\n')
    fprintf('Successfully pulled sleep and wake from this subject: %s\n',wake(ifile).name(1:end-12))
    fprintf('\n')
    W.mse(isinf(W.mse)) = NaN;
    S.mse(isinf(S.mse)) = NaN;
    
    assert(W.cfg.new_srate == S.cfg.new_srate, ...
       'Sleep and wake have different sampling rates ...')
    
    assert(W.cfg.m == S.cfg.m, 'Sleep and wake have different embedding dimension ...')
    
    % extract good windows, if applicable 
    
    min_valid_mse = 100; % Grandy et al., 2016
    
    wk_gdwin = W.n_valid >= min_valid_mse;
    sp_gdwin = S.n_valid >= min_valid_mse;
    if sum(wk_gdwin) == 0 || sum(sp_gdwin) == 0 % at least one good window!
        fprintf('Skipping this file, no good windows: %i\n', SID)
        continue
    end
    
    %% Load LZC

    % load empirical values
    %try
        Welzc = load(sprintf('%swake/%s%s',pthLZC,fstr,'_wake_LZC'));
    %catch
    %    fprintf('\n This subject is missing Lempel-Ziv\n')
    %    continue
    %end
    if isfield(Welzc.LZCout,'surrogate') % not all files will have this field
        assert(~Welzc.LZCout.surrogate,'This empirical file is actually a surrogate')
    end
    assert(Welzc.LZCout.cfg.new_srate == fsample,'File wasn''t downsampled?')
    
    
    if surrtest
        % load surrogate values
        Wslzc = load(sprintf('%swake/Surrogate/%s',pthLZC,sprintf('%s%s',...
            fstr,'_wake_LZC_SURROGATE')));
        if isfield(Wslzc.LZCout,'surrogate') % not all files will have this field
            assert(Wslzc.LZCout.surrogate,'This surrogate file is actually empirical data')
        end
        assert(Wslzc.LZCout.cfg.new_srate == fsample,'File wasn''t downsampled?')
    end
    
    % Build new structure by combining info from empirical and surrogate
    % data
    
    switch surrtest
        case true
            Wlzc.gLZC = Welzc.LZCout.gLZC-Wslzc.LZCout.gLZC; % take difference
            Wlzc.vLZC = Welzc.LZCout.vLZC-Wslzc.LZCout.vLZC; % take difference
            Wlzc.n_valid_gLZC = min([Welzc.LZCout.n_valid_gLZC; Wslzc.LZCout.n_valid_gLZC]);
            Wlzc.n_valid_vLZC = min([Welzc.LZCout.n_valid_vLZC; Wslzc.LZCout.n_valid_vLZC]);
        case false
            Wlzc.gLZC = Welzc.LZCout.gLZC; 
            Wlzc.vLZC = Welzc.LZCout.vLZC;
            Wlzc.n_valid_gLZC = Welzc.LZCout.n_valid_gLZC;
            Wlzc.n_valid_vLZC = Welzc.LZCout.n_valid_vLZC;
    end
    Wlzc.cfg = Welzc.LZCout.cfg;
    Wlzc.prc_used = Welzc.LZCout.prc_used;
    Wlzc.dur_used = Welzc.LZCout.dur_used;
    if match
        % targeted comparison sanity check
        assert(Welzc.LZCout.cfg.wk_len == W.wk_len,'Different sections matched')
        assert(Welzc.LZCout.cfg.wk_segs_used == W.wk_segs_used,'Different sections matched')
    end
    fprintf('\nLoading awake file %s \n',fstr)
    
    try % attempt to load the corresponding sleep file
        sleepfile = strcat(fstr,'_sleep_LZC');
        
        % Load empirical values
        Selzc = load(strcat(sprintf('%ssleep/',pthLZC),sleepfile));
        
        fprintf('\n Loading sleep file %s \n',sleepfile)
    catch
        fprintf('\n No sleep data for this subject ....\n')
        continue
    end
    
    if isfield(Selzc.LZCout,'surrogate') % not all files will have this field
        assert(~Selzc.LZCout.surrogate,'This empirical file is actually a surrogate')
    end
    assert(Selzc.LZCout.cfg.new_srate == fsample,'File wasn''t downsampled?')
    
    if surrtest
        % Load surrogate values
        Sslzc = load(sprintf('%ssleep/Surrogate/%s%s.mat',pthLZC,...
            sleepfile,'_SURROGATE'));     
        if isfield(Sslzc.LZCout,'surrogate') % not all files will have this field
            assert(Sslzc.LZCout.surrogate,'This surrogate file is actually empirical data')
        end
        assert(Sslzc.LZCout.cfg.new_srate == fsample,'File wasn''t downsampled?')
    end

    % Build new structure by combining info from empirical and surrogate
    % data
    
    switch surrtest
        case true
            Slzc.gLZC = Selzc.LZCout.gLZC-Sslzc.LZCout.gLZC; % take difference
            Slzc.vLZC = Selzc.LZCout.vLZC-Sslzc.LZCout.vLZC; % take difference
            Slzc.n_valid_gLZC = min([Selzc.LZCout.n_valid_gLZC; Sslzc.LZCout.n_valid_gLZC]);
            Slzc.n_valid_vLZC = min([Selzc.LZCout.n_valid_vLZC; Sslzc.LZCout.n_valid_vLZC]);
        case false
            Slzc.gLZC = Selzc.LZCout.gLZC; 
            Slzc.vLZC = Selzc.LZCout.vLZC; 
            Slzc.n_valid_gLZC = Selzc.LZCout.n_valid_gLZC;
            Slzc.n_valid_vLZC = Selzc.LZCout.n_valid_vLZC;
    end
    Slzc.cfg = Selzc.LZCout.cfg;
    Slzc.prc_used = Selzc.LZCout.prc_used;
    Slzc.dur_used = Selzc.LZCout.dur_used;
    if match
        % targeted comparison sanity check
        assert(Selzc.LZCout.cfg.sp_len == S.sp_len,'Different sections matched')
        assert(Selzc.LZCout.cfg.sp_segs_used == S.sp_segs_used,'Different sections matched')
    end
    
    if ~ismember(SID,IDs)
        IDs = [IDs; SID];
    end
    idx = find(IDs == SID);
    fprintf('\n')
    fprintf('Successfully pulled sleep and wake from this subject: %s\n',fstr)
    fprintf('\n')
    
    % remove infinities 
    Wlzc.vLZC(isinf(Wlzc.vLZC)) = NaN;
    Slzc.vLZC(isinf(Slzc.vLZC)) = NaN;
    Wlzc.gLZC(isinf(Wlzc.gLZC)) = NaN;
    Slzc.gLZC(isinf(Slzc.gLZC)) = NaN;
    
    % discard extraneous entries for n_valid, since already know it's the
    % same across channels
    Wlzc.n_valid_vLZC = Wlzc.n_valid_vLZC(1:size(Wlzc.vLZC,1):end,:);
    Wlzc.n_valid_gLZC = Wlzc.n_valid_vLZC(1:size(Wlzc.gLZC,1):end,:);
    Slzc.n_valid_vLZC = Slzc.n_valid_vLZC(1:size(Slzc.vLZC,1):end,:);
    Slzc.n_valid_gLZC = Slzc.n_valid_vLZC(1:size(Slzc.gLZC,1):end,:);
    % check frequency vector and smoothing for awake
    if ~exist('LZC_foi','var')
        LZC_foi = Wlzc.cfg.foi;
        % scales we want to look at (don't go far below 1 Hz)
        scales = find(LZC_foi>= 0.9); 
    else
        assert(sum(LZC_foi==Wlzc.cfg.foi) == length(LZC_foi),'Frequency vectors don''t match');
    end
    if ~exist('smooth_lo','var')
        smooth_lo = Wlzc.cfg.smooth_lo;
    else
        assert(sum(smooth_lo==Wlzc.cfg.smooth_lo) == length(smooth_lo),'Smooth_lo vectors don''t match');
    end
    if ~exist('smooth_hi','var')
        smooth_hi = Wlzc.cfg.smooth_hi;
    else
        assert(sum(smooth_hi==Wlzc.cfg.smooth_hi) == length(smooth_hi),'Smooth_hi vectors don''t match');
    end
    % check frequency vector and smoothing for asleep
    assert(sum(LZC_foi==Slzc.cfg.foi) == length(LZC_foi),'Frequency vectors don''t match');
    assert(sum(smooth_lo==Slzc.cfg.smooth_lo) == length(smooth_lo),'Smooth_lo vectors don''t match');
    assert(sum(smooth_hi==Slzc.cfg.smooth_hi) == length(smooth_hi),'Smooth_hi vectors don''t match');
    
    % sanity checks below
    assert(Wlzc.cfg.new_srate == Slzc.cfg.new_srate, ...
        'Sleep and wake have different sampling rates ...')
    
    method = 'gMLZ';
    min_valid_lzc = 2000; % (shortest number of data points we'll accept)
    
    % find good windows for multiscale (gMLZ or dLZC)
    switch method
        case 'gMLZ'
            wk_lzc_gdwin = Wlzc.n_valid_gLZC >= min_valid_lzc;
            sp_lzc_gdwin = Slzc.n_valid_gLZC >= min_valid_lzc;
            wk_lzc_win_avg = nanmean(Wlzc.gLZC(:,:,wk_lzc_gdwin),3);
            sp_lzc_win_avg = nanmean(Slzc.gLZC(:,:,sp_lzc_gdwin),3);
        case 'dLZC'
            wk_lzc_gdwin = Wlzc.n_valid_dLZC(:,1) >= min_valid_lzc;
            sp_lzc_gdwin = Slzc.n_valid_dLZC(:,1) >= min_valid_lzc;
            wk_lzc_win_avg = nanmean(Wlzc.dLZC(:,:,wk_lzc_gdwin),3);
            sp_lzc_win_avg = nanmean(Slzc.dLZC(:,:,sp_lzc_gdwin),3);
        otherwise
            error('Method unknown')
    end
    
     % find good windows for Vanilla LZC
     wk_lzc_ss_gdwin = Wlzc.n_valid_vLZC >= min_valid_lzc;
     sp_lzc_ss_gdwin = Slzc.n_valid_vLZC >= min_valid_lzc;
    
    if sum(wk_lzc_gdwin) == 0 || sum(sp_lzc_gdwin) == 0 % at least one good window!
        fprintf('Skipping this file, no good windows: %i\n', SlzcID)
        bcnt = bcnt + 1;
        continue
    end
    
    %% Store file IDs so we know which files to pull in unmatched analysis
    % find underscores in file name
    underscore = strfind(wake(ifile).name,'_');
    stop = underscore(3)-1; % everything up to the third underscore 
    fstr = wake(ifile).name(1:stop);
    if idx > 1 && idx <= length(fileID)
        fileID{idx}{length(fileID{idx})+1} = fstr;
    else
        fileID{idx} = cell(1);
        fileID{idx}{1} = fstr;
    end
    clear fstr
      
    wk_win_avg = nanmean(W.mse(:,:,wk_gdwin),3);
    sp_win_avg = nanmean(S.mse(:,:,sp_gdwin),3);
    if idx <= length(mse_wk) && ~isempty(mse_wk{idx})
        mse_wk{idx} = cat(3,mse_wk{idx},wk_win_avg);
        allage{idx} = cat(2,allage{idx},age);
        dur_wk_all{idx} = cat(2,dur_wk_all{idx},W.dur_used);
        dur_sp_all{idx} = cat(2,dur_sp_all{idx},S.dur_used);
        if match
            trg_wk_all{idx} = cat(2,trg_wk_all{idx},W.wk_len);
            trg_sp_all{idx} = cat(2,trg_sp_all{idx},S.sp_len);
            trg_wk_segs{idx} = cat(2,trg_wk_segs{idx},W.wk_segs_used);
            trg_sp_segs{idx} = cat(2,trg_sp_segs{idx},S.sp_segs_used);
        end
        win_wk{idx} = cat(2,win_wk{idx},sum(W.n_valid(wk_gdwin)));
        win_sp{idx} = cat(2,win_sp{idx},sum(S.n_valid(sp_gdwin)));
        fprintf('Adding another visit\n')
    else
        mse_wk{idx} = wk_win_avg;
        allage{idx} = age;
        allgender(idx) = gender;
        allgenotype(idx) = genotype;
        win_wk{idx} = sum(W.n_valid(wk_gdwin));
        win_sp{idx} = sum(S.n_valid(sp_gdwin));
        dur_wk_all{idx} = W.dur_used;
        dur_sp_all{idx} = S.dur_used;
        if match
            trg_wk_all{idx} = W.wk_len;
            trg_sp_all{idx} = S.sp_len;
            trg_wk_segs{idx} = W.wk_segs_used;
            trg_sp_segs{idx} = S.sp_segs_used;
        end
        SID_wk{idx} = SID;
        SID_sp{idx} = SID;
    end
    if idx > length(allcog) || isempty(allcog(idx))
        if ~isempty(BCog)
            allcog{idx} = BCog;
        else
            allcog{idx} = NaN;
        end
    else
        idx2 = length(allcog{idx})+1;
        if ~isempty(BCog)
            allcog{idx}(idx2) = BCog;
        else
            allcog{idx}(idx2) = NaN;
        end
    end
    if idx > length(alllng) || isempty(alllng(idx))
        alllng{idx} = BLang;
        if ~isempty(BLang)
            alllng{idx} = BLang;
        else
            alllng{idx} = NaN;
        end
    else
        idx2 = length(alllng{idx})+1;
        if ~isempty(BLang)
            alllng{idx}(idx2) = BLang;
        else
            alllng{idx}(idx2) = NaN;
        end
    end

    % do for LZC
    
    min_valid_lzc = 2000; % 2015 mLZC paper
    
    wk_lzc_win_avg = nanmean(Wlzc.gLZC(:,:,wk_lzc_gdwin),3);
    sp_lzc_win_avg = nanmean(Slzc.gLZC(:,:,sp_lzc_gdwin),3);
    if idx <= length(lzc_wk) && ~isempty(lzc_wk{idx})
        lzc_wk{idx} = cat(3,lzc_wk{idx},wk_lzc_win_avg);
        lzc_sp{idx} = cat(3,lzc_sp{idx},sp_lzc_win_avg);
        dur_wk_all_lzc{idx} = cat(2,dur_wk_all_lzc{idx},Wlzc.dur_used);
        dur_sp_all_lzc{idx} = cat(2,dur_sp_all_lzc{idx},Slzc.dur_used);
        win_wk_lzc{idx} = cat(2,win_wk_lzc{idx},sum(Wlzc.n_valid_gLZC(wk_lzc_gdwin)));
        win_sp_lzc{idx} = cat(2,win_sp_lzc{idx},sum(Slzc.n_valid_gLZC(sp_lzc_gdwin)));
        fprintf('Adding LZC another visit\n')
    else
        lzc_wk{idx} = wk_lzc_win_avg;
        lzc_sp{idx} = sp_lzc_win_avg;
        win_wk_lzc{idx} = sum(Wlzc.n_valid_gLZC(wk_lzc_gdwin));
        win_sp_lzc{idx} = sum(Slzc.n_valid_gLZC(sp_lzc_gdwin));
        dur_wk_all_lzc{idx} = Wlzc.dur_used;
        dur_sp_all_lzc{idx} = Slzc.dur_used;
        SID_wk_lzc{idx} = SID;
        SID_sp_lzc{idx} = SID;
    end
    
    % Now look at the other cognitive scale (from paper supplement)
    if idx > 1 && idx <= length(allVCog)
        allVCog{idx}  = [allVCog{idx} VCog];
    else
        allVCog{idx}  = VCog;
    end
    
    if idx > 1 && idx <= length(allNVCog)
        allNVCog{idx}  = [allVCog{idx} NVCog];
    else
        allNVCog{idx}  = NVCog;
    end
    
    if idx <= length(mse_sp) && ~isempty(mse_sp{idx})
        mse_sp{idx} = cat(3,mse_sp{idx},sp_win_avg);
        fprintf('Stacking another visit\n')
    else
        mse_sp{idx} = sp_win_avg;
    end
    clear isubj idate iage date age year month day
end

% clear genotype, age, etc. where we didn't get power

badidx = [];
% find emtpty indicies 
for i = 1:length(mse_wk)
    if isempty(mse_wk{i})
        % sanity check
        assert(isempty(all_pow_wk{i}),'Wake power is not empty for this subject!')
        badidx = [badidx i];
    end
end

%clear entries
select_files = cell(1,1); % list of the files we actually uses
allgenotype(badidx) = [];
allgender(badidx) = [];
allage(badidx) = [];
allNVCog(badidx) = [];
allVCog(badidx) = [];
allcog(badidx) = [];
alllng(badidx) = [];
all_pow_wk(badidx) = [];
all_pow_sp(badidx) = [];
%all_sec_wk(badidx) = [];
%all_sec_sp(badidx) = [];

% MSE
mse_wk(badidx) = [];
mse_sp(badidx) = [];
win_wk(badidx) = [];
win_sp(badidx) = [];
dur_wk_all(badidx) = [];
dur_sp_all(badidx) = [];
SID_sp(badidx) = [];
SID_wk(badidx) = [];
% LZC
lzc_wk(badidx) = [];
lzc_sp(badidx) = [];
win_wk_lzc(badidx) = [];
win_sp_lzc(badidx) = [];
dur_wk_all_lzc(badidx) = [];
dur_sp_all_lzc(badidx) = [];
SID_sp_lzc(badidx) = [];
SID_wk_lzc(badidx) = [];



% average across scales
% select visit with greatest number of good windows

wk_avg = nan(19,length(mse_wk)); % memory allocation for average across scales
sp_avg = nan(19,length(mse_sp)); % memory allocation for average across scales

wk_scl = nan(length(scales),length(mse_wk)); % preserve scales (average across channels instead)
sp_scl = nan(length(scales),length(mse_sp)); % preserve scales (average across channels instead)

wk_avg_lzc = nan(19,length(lzc_wk)); % memory allocation for average across scales
sp_avg_lzc = nan(19,length(lzc_sp)); % memory allocation for average across scales

wk_scl_lzc = nan(length(scales),length(lzc_wk)); % preserve scales (average across channels instead)
sp_scl_lzc = nan(length(scales),length(lzc_sp)); % preserve scales (average across channels instead)

wk_all = []; % should probably be allocated ... fix later
sp_all = [];
wk_lzc_all = [];
sp_lzc_all = [];

assert(length(mse_wk)==length(lzc_wk),'LZC and MSE are different sizes')
lcnt = 0; % number of subjects with longitudinal data

for i = 1:length(mse_wk)
    if ~isempty(mse_wk{i})
        if ndims(mse_wk{i}) == 3 % if there's longitudinal data
            lcnt = lcnt + 1;
            assert(ndims(all_pow_wk{i})==2,'Only one power file')
            % get number of good 1 Hz windows
            Nwn = [];
            for j = 1:size(all_pow_wk{i},2)
                Nwn = [Nwn all_pow_wk{i}(j).n(1)];
            end
            [g,gdx] = max(Nwn);
            if sum(Nwn==g) > 1 % if there's more than one max, take youngest
                gage = allage{i}(Nwn==g);
                minage = min(gage);
                gdx = find(allage{i}(:)==minage); % overwrite good index
            end
            tmp_wk = mse_wk{i}(:,:,gdx); % take visit with either best power or youngest age (if power tie)
            tmp_sp = mse_sp{i}(:,:,gdx);
            tmp_wk_lzc = lzc_wk{i}(:,:,gdx);
            tmp_sp_lzc = lzc_sp{i}(:,:,gdx);
            all_ages(i) = allage{i}(gdx);
            cog_score(i) = allcog{i}(gdx);
            lng_score(i) = alllng{i}(gdx);
            VC_score(i) = allVCog{i}(gdx);
            NVC_score(i) = allNVCog{i}(gdx);
            wk_pow{i} = all_pow_wk{i}(gdx);
            sp_pow{i} = all_pow_sp{i}(gdx);
            %wk_sec(i) = all_sec_wk{i}(gdx);
            %sp_sec(i) = all_sec_sp{i}(gdx);
            dur_wk(i) = dur_wk_all{i}(gdx);
            dur_sp(i) = dur_sp_all{i}(gdx);
            if match
                trg_wk(i) = trg_wk_all{i}(gdx);
                trg_sp(i) = trg_sp_all{i}(gdx);
                trg_wk_seg(i) = trg_wk_segs{i}(gdx);
                trg_sp_seg(i) = trg_sp_segs{i}(gdx);
            end
            select_files{i} = fileID{i}{gdx};
            % MSE
            wk_avg(:,i) = mean(tmp_wk(:,scales),2); % mean across scales
            wk_scl(:,i) = mean(tmp_wk(:,scales),1); % mean across channels
            wk_all      = cat(3,wk_all,tmp_wk);
            sp_avg(:,i) = mean(tmp_sp(:,scales),2);
            sp_scl(:,i) = mean(tmp_sp(:,scales),1);
            sp_all      = cat(3,sp_all,tmp_sp);
            % LZC
            wk_avg_lzc(:,i) = mean(tmp_wk_lzc(:,scales),2); % mean across scales
            wk_scl_lzc(:,i) = mean(tmp_wk_lzc(:,scales),1); % mean across channels
            wk_lzc_all      = cat(3,wk_lzc_all,tmp_wk_lzc);
            sp_avg_lzc(:,i) = mean(tmp_sp_lzc(:,scales),2);
            sp_scl_lzc(:,i) = mean(tmp_sp_lzc(:,scales),1);
            sp_lzc_all      = cat(3,sp_lzc_all,tmp_sp_lzc);
        else
            % MSE 
            wk_avg(:,i) = mean(mse_wk{i}(:,scales),2);
            wk_scl(:,i) = mean(mse_wk{i}(:,scales),1);
            wk_all      = cat(3,wk_all,mse_wk{i});
            sp_avg(:,i) = mean(mse_sp{i}(:,scales),2);
            sp_scl(:,i) = mean(mse_sp{i}(:,scales),1);
            sp_all      = cat(3,sp_all,mse_sp{i});
            % LZC
            wk_avg_lzc(:,i) = mean(lzc_wk{i}(:,scales),2);
            wk_scl_lzc(:,i) = mean(lzc_wk{i}(:,scales),1);
            wk_lzc_all      = cat(3,wk_lzc_all,lzc_wk{i});
            sp_avg_lzc(:,i) = mean(lzc_sp{i}(:,scales),2);
            sp_scl_lzc(:,i) = mean(lzc_sp{i}(:,scales),1);
            sp_lzc_all      = cat(3,sp_lzc_all,lzc_sp{i});
            all_ages(i) = allage{i};
            cog_score(i) = allcog{i};
            lng_score(i) = alllng{i};
            VC_score(i) = allVCog{i};
            NVC_score(i) = allNVCog{i};
            wk_pow{i} = all_pow_wk{i};
            sp_pow{i} = all_pow_sp{i};
            %wk_sec(i) = all_sec_wk{i};
            %sp_sec(i) = all_sec_sp{i};
            dur_wk(i) = dur_wk_all{i};
            dur_sp(i) = dur_sp_all{i};
            if match
                trg_wk(i) = trg_wk_all{i};
                trg_sp(i) = trg_sp_all{i};
                trg_wk_seg(i) = trg_wk_segs{i};
                trg_sp_seg(i) = trg_sp_segs{i};
            end
            select_files{i} = fileID{i}{1};
        end
    end
end

wk_prnan = isnan(sum(wk_scl));
wk_scl(:,wk_prnan) = []; % prune nans;
sp_prnan = isnan(sum(sp_scl));
sp_scl(:,sp_prnan) = []; % prune nans

%% Report sample size

if match % otherwise not all longitudinal data gets loaded
    fprintf('\n %i participants have longitudinal data \n',lcnt)
end

fprintf('\n N = %i participants\n',size(wk_scl,2))

if match
    % save the files we used 
    save select_files select_files
end


%% plot average MSE scales

%figure('Color','w')
figure('Color','w')
hold on
[~,~,ci] = ttest(wk_scl');
plot(log2(nyqfoi),mean(wk_scl,2),'r','linewidth',4)
patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
    1.0,'EdgeColor','none')
[~,~,ci] = ttest(sp_scl');
plot(log2(nyqfoi),mean(sp_scl,2),'b','linewidth',4)
patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [ci(1,:) fliplr(ci(2,:))], 'b','facealpha', ...
    1.0,'EdgeColor','none')
legend('wake','wake 95% CI','sleep','95% sleep','AutoUpdate','off',...
    'FontSize',14,'location','northeast')
legend boxoff
%title('MSE values over all time scales (no regression)','fontsize',18)
xlim([log2(5) log2(100)])
ticks = 2:1:7;
xticks(ticks)
ticklabels = cell(1,length(ticks));
for i = 1:length(ticklabels)
    ticklabels{i} = num2str(2^ticks(i));
end
xticklabels(ticklabels)
xlabel('Nyquist frequency (Hz)')
ylabel('SampEn')
makefigpretty
% ADD TOP X-AXIS FOR TAU
ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation',...
    'right', 'Color','none');
makefigpretty
line(log2(nyqfoi),mean(wk_scl,2),'Parent',ax2,'Color','none')
xlim([log2(5) log2(100)])
pltscl = fliplr([1 2 5 10 20]);
ticks = log2(nyqfoi(pltscl));
xticks(ticks)
ticklabels = cell(1,length(ticks));
for i = 1:length(ticklabels)
    ticklabels{i} = num2str(pltscl(i));
end
xticklabels(ticklabels)
xlabel('time scale')
set(gca,'ytick',[]);
yAX = get(gca,'YAxis');
yAX.Visible = 'off';
% SAVE FIGURE
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'MSE_scales_no_reg'))
print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,'MSE_scales_no_reg'))

%% plot MSE percent change

figure('Color','w')
hold on
wsdelta = ((mean(wk_scl,2) - mean(sp_scl,2))./mean(sp_scl,2)).*100;
[~,~,ci] = ttest(wk_scl,sp_scl,'dim',2); ci = ci';
plot(log2(nyqfoi),wsdelta ,'m','linewidth',4)
patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [(ci(1,:)./mean(sp_scl,2)').*100 ...
    fliplr((ci(2,:)./mean(sp_scl,2)').*100)], 'm','facealpha', 1.0,'EdgeColor','none')
legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',14,'location','southeast')
legend boxoff
plot(1:length(wsdelta),zeros(1,length(wsdelta)),'-.k')
xlim([log2(5) log2(100)])
ticks = 2:1:7;
xticks(ticks)
ticklabels = cell(1,length(ticks));
for i = 1:length(ticklabels)
    ticklabels{i} = num2str(2^ticks(i));
end
xticklabels(ticklabels)
%title('MSE difference over all time scales (no regression)','fontsize',18)
xlabel('Nyquist frequency (Hz)')
ylabel('SampEn change from sleep (%)')
makefigpretty
% ADD TOP X-AXIS FOR TAU
ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation',...
    'right', 'Color','none');
makefigpretty
line(log2(nyqfoi),wsdelta,'Parent',ax2,'Color','none')
xlim([log2(5) log2(100)])
pltscl = fliplr([1 2 5 10 20]);
ticks = log2(nyqfoi(pltscl));
xticks(ticks)
ticklabels = cell(1,length(ticks));
for i = 1:length(ticklabels)
    ticklabels{i} = num2str(pltscl(i));
end
xticklabels(ticklabels)
xlabel('time scale')
set(gca,'ytick',[]);
yAX = get(gca,'YAxis');
yAX.Visible = 'off';
% SAVE FIGURE
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'MSE_scales_delta_no_reg'))
print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,'MSE_scales_delta_no_reg'))


%% MSE permutation clustering 

clst = RMpermclusttest(wk_all,sp_all,Nperm);

sgcl_pos = [];
sgcl_neg = [];

if sum(clst.n_pos > 0) > 0 
    sgcl_pos = clst.P_val_pos < 0.05;
end

if sum(clst.n_neg > 0) > 0 
    sgcl_neg = clst.P_val_neg < 0.05;
end
cnt = 0;


for icls = 1:length(sgcl_pos)
    if sgcl_pos(icls)
        cnt = cnt +1;
        figure('Color','w')
        plot(log2(nyqfoi),sum(clst.lbl_pos==icls,2),'linewidth',2)
        ylim([0 20])
        xlim([log2(5) log2(100)])
        ticks = 2:1:7;
        xticks(ticks)
        ticklabels = cell(1,length(ticks));
        for i = 1:length(ticklabels)
            ticklabels{i} = num2str(2^ticks(i));
        end
        xticklabels(ticklabels)
        xlabel('Nyquist frequency (Hz)')
        ylabel('Number of channels')
        makefigpretty
        ax1 = gca; % current axes
        ax1_pos = ax1.Position; % position of first axes
        ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation',...
            'right', 'Color','none');
        makefigpretty
        line(log2(nyqfoi),wsdelta,'Parent',ax2,'Color','none')
        xlim([log2(5) log2(100)])
        pltscl = fliplr([1 2 5 10 20]);
        ticks = log2(nyqfoi(pltscl));
        xticks(ticks)
        ticklabels = cell(1,length(ticks));
        for i = 1:length(ticklabels)
            ticklabels{i} = num2str(pltscl(i));
        end
        xticklabels(ticklabels)
        xlabel('time scale')
        set(gca,'ytick',[]);
        yAX = get(gca,'YAxis');
        yAX.Visible = 'off';
        title(sprintf('MSE postivie cluster, p = %1.3f',clst.P_val_pos(icls)),'fontsize',18)
        figname = sprintf('MSE pos clust spect #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,figname))
        
        figure('Color','w')
        plot_topo_AS(sum(clst.lbl_pos==icls)')
        title(sprintf('MSE pos clust, p = %1.3f',clst.P_val_pos(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(wk_lzc_all,2)])
        mycolorbar
        figname = sprintf('MSE pos clust toop #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,figname))
    end
end

cnt = 0;
for icls = 1:length(sgcl_neg)
    if sgcl_neg(icls)
        cnt = cnt + 1;
        figure('Color','w')
        plot(log2(nyqfoi),sum(clst.lbl_neg==icls,2),'linewidth',2)
        ylim([0 20])
        xlim([log2(5) log2(100)])
        ticks = 2:1:7;
        xticks(ticks)
        ticklabels = cell(1,length(ticks));
        for i = 1:length(ticklabels)
            ticklabels{i} = num2str(2^ticks(i));
        end
        xticklabels(ticklabels)
        xlabel('Nyquist frequency (Hz)')
        ylabel('Number of channels')
        makefigpretty
        ax1 = gca; % current axes
        ax1_pos = ax1.Position; % position of first axes
        ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation',...
            'right', 'Color','none');
        line(log2(nyqfoi),wsdelta,'Parent',ax2,'Color','none')
        xlim([log2(5) log2(100)])
        pltscl = fliplr([1 2 5 10 20]);
        ticks = log2(nyqfoi(pltscl));
        xticks(ticks)
        ticklabels = cell(1,length(ticks));
        for i = 1:length(ticklabels)
            ticklabels{i} = num2str(pltscl(i));
        end
        xticklabels(ticklabels)
        xlabel('time scale')
        set(gca,'ytick',[]);
        yAX = get(gca,'YAxis');
        yAX.Visible = 'off';
        makefigpretty
        title(sprintf('MSE negative cluster, p = %1.3f',clst.P_val_neg(icls)),'fontsize',18)
        figname = sprintf('MSE neg clust spect #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,figname))
        
        figure('Color','w')
        plot_topo_AS(sum(clst.lbl_neg==icls)')
        title(sprintf('MSE neg clust, p = %1.3f',clst.P_val_neg(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(wk_lzc_all,2)])
        mycolorbar
        figname = sprintf('MSE negs clust topo #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,figname))
    end
end

%% plot average gMLZ scales

figure('Color','w')
h=subplot(1,1,1);
hold on
[~,~,ci] = ttest(wk_scl_lzc');
plot(log2(LZC_foi(scales)),mean(wk_scl_lzc,2),'r','linewidth',4)
patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
    1.0,'EdgeColor','none')
[~,~,ci] = ttest(sp_scl_lzc');
plot(log2(LZC_foi(scales)),mean(sp_scl_lzc,2),'b','linewidth',4)
patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], 'b','facealpha', ...
    1.0,'EdgeColor','none')
legend('wake','wake 95% CI','sleep','95% sleep','AutoUpdate','off','FontSize',14,'location','northwest')
legend boxoff
ticks = 0:1:5;
xticks(ticks)
ticklabels = cell(1,length(ticks));
for i = 1:length(ticklabels)
    ticklabels{i} = num2str(2^ticks(i));
end
xticklabels(ticklabels)
%title('MSE difference over all time scales (no regression)','fontsize',18)
xlabel('Center frequency (Hz)')
ylabel('LZ complexity')
makefigpretty
% ADD TOP X-AXIS FOR TAU
ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation',...
    'right', 'Color','none');
makefigpretty
line(log2(LZC_foi),wsdelta,'Parent',ax2,'Color','none')
xlim([log2(0.9) log2(LZC_foi(end))])
pltscl = [1 5 10 15 20];
ticks = log2(LZC_foi(pltscl));
xticks(ticks)
ticklabels = cell(1,length(ticks));
for i = 1:length(ticklabels)
    ticklabels{i} = num2str(pltscl(i));
end
xticklabels(ticklabels)
xlabel('time scale')
set(gca,'ytick',[]);
yAX = get(gca,'YAxis');
yAX.Visible = 'off';
xlim([log2(0.9) log2(LZC_foi(end))])
title('gMLZ values  covary nuisance over all time scales (no regression)','fontsize',18)
xlabel('Center frequency (Hz)')
ylabel('LZ complexity')
set(h,'XTick',-1:1:5,'XTickLabel',2.^(-1:1:5))
xlim([log2(0.9) 5])
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_scales_no_reg'))

%% gMLZ percent change from sleep 

figure('Color','w')
h=subplot(1,1,1);
hold on
wsdelta = (mean(wk_scl_lzc,2) - mean(sp_scl_lzc,2))./mean(sp_scl_lzc,2).*100;
[~,~,ci] = ttest(wk_scl_lzc,sp_scl_lzc,'dim',2); ci = ci';
plot(log2(LZC_foi(scales)),wsdelta ,'m','linewidth',4)
patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , ...
    [ci(1,:)./mean(sp_scl_lzc,2)'.*100 fliplr(ci(2,:)./mean(sp_scl_lzc,2)'.*100)], 'r','facealpha', ...
    1.0,'EdgeColor','none')
legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',14,'location','northwest')
legend boxoff
plot(log2(LZC_foi(scales)),zeros(1,length(wsdelta)),'-.k')
xlim([0 5])
%ylim([-0.04 0.04])
title('gMLZ difference over all time scales (no regression)','fontsize',18)
xlabel('Center frequency (Hz)')
ylabel('gLZ change from sleep (%)')
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_scales_delta_no_reg'))


%% gMLZ permutation clustering

clst_LZC = RMpermclusttest(wk_lzc_all,sp_lzc_all,Nperm);

sgcl_pos = [];
sgcl_neg = [];

if sum(clst_LZC.n_pos > 0) > 0 
    sgcl_pos = clst_LZC.P_val_pos < 0.05;
end

if sum(clst_LZC.n_neg > 0) > 0 
    sgcl_neg = clst_LZC.P_val_neg < 0.05;
end
cnt = 0;

for icls = 1:length(sgcl_pos)
    if sgcl_pos(icls)
        cnt = cnt + 1;
        figure('Color','w')
        h=subplot(1,1,1);
        plot(log2(LZC_foi(scales)),sum(clst_LZC.lbl_pos==icls,2),'linewidth',2)
        ylim([0 20])
        xlim([0 5])
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        xlim([0 5])
        xlabel('Center frequency (Hz)')
        ylabel('Number of channels')
        makefigpretty
        title(sprintf('gMLZ postivie cluster covary nuisance, p = %1.3f',clst_LZC.P_val_pos(icls)),'fontsize',18)
        figname = sprintf('LZC pos clust spect #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        
        figure('Color','w')
        plot_topo_AS(sum(clst_LZC.lbl_pos==icls)')
        title(sprintf('gMLZ pos clust, p = %1.3f',clst_LZC.P_val_pos(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(wk_lzc_all,2)])
        mycolorbar
        figname = sprintf('LZC pos clust topo #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
    end
end
cnt = 0;
for icls = 1:length(sgcl_neg)
    if sgcl_neg(icls)
        cnt = cnt+1;
        figure('Color','w')
        h=subplot(1,1,1);
        plot(log2(LZC_foi(scales)),sum(clst_LZC.lbl_neg==icls,2),'linewidth',2)
        ylim([0 20])
        xlim([0 5])
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        xlim([0 5])
        xlabel('Center frequency (Hz)')
        ylabel('Number of channels')
        makefigpretty
        title(sprintf('gMLZ negative cluster covary nuisance, p = %1.3f',clst_LZC.P_val_neg(icls)),'fontsize',18)
        figname = sprintf('gMLZ neg clust spect #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        
        figure('Color','w')
        plot_topo_AS(sum(clst_LZC.lbl_neg==icls)')
        title(sprintf('gMLZ neg clust, p = %1.3f',clst_LZC.P_val_neg(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(WP,2)])
        mycolorbar
        figname = sprintf('gMLZ neg clust topo #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
    end
end



%% Do scatter plots of age 

figure('Color','w')
scatter(log2(all_ages),mean(wk_avg),'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 1])
[r,p] = corrcoef(all_ages,mean(wk_avg));
xlabel('age log2(months)')
ylabel('awake MSE')
mylsline
%ylim([min(mean(wk_avg)) max(mean(wk_avg))])
title(sprintf('r = %1.3f, p = %1.3f, scales = %i-%i',r(2,1),p(2,1),scales(1),scales(end)),'fontsize',18)
makefigpretty
figname = sprintf('age_wk_mse_%i_%i',scales(1),scales(end));
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

figure('Color','w')
scatter(log2(all_ages),mean(sp_avg),'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 1])
hold on
[r,p] = corrcoef(all_ages,mean(sp_avg));
xlabel('age log2(months)')
ylabel('asleep MSE')
mylsline
%ylim([min(mean(sp_avg)) max(mean(sp_avg))])
title(sprintf('r = %1.3f, p = %1.3f, scales = %i-%i',r(2,1),p(2,1),scales(1),scales(end)),'fontsize',18)
makefigpretty
figname = sprintf('age_sp_mse_%i_%i',scales(1),scales(end));
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

% plot R^2 as a function of time scale
Rwk = nan(1,size(wk_avg,1));
Rsp = nan(1,size(sp_avg,1));
for iscale = 1:size(wk_avg,1)
    rwk = corrcoef(wk_avg(iscale,:),all_ages,'rows','complete');
    rsp = corrcoef(sp_avg(iscale,:),all_ages,'rows','complete');
    Rwk(iscale) = rwk(2,1)^2;
    Rsp(iscale) = rsp(2,1)^2;
end

figure('Color','w')
plot([Rwk; Rsp]','linewidth',2)
xlabel('Time scale')
ylabel('R^{2}')
legend('Wake','Sleep')
legend boxoff
title('Variance of MSE explained by age')
makefigpretty
figname = 'R2MSEage';
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

%% Do scatter plots of age and MSE 

figure('Color','w')
scatter(log2(all_ages),mean(wk_avg_lzc),'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 1])
[r,p] = corrcoef(all_ages,mean(wk_avg_lzc));
xlabel('age log2(months)')
ylabel('awake LZC')
mylsline
%ylim([min(mean(wk_avg_lzc)) max(mean(wk_avg_lzc))])
title(sprintf('r = %1.3f, p = %1.3f, scales = %i-%i',r(2,1),p(2,1),scales(1),scales(end)),'fontsize',18)
makefigpretty
figname = sprintf('age_wk_LZC_%i_%i',scales(1),scales(end));
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

figure('Color','w')
scatter(log2(all_ages),mean(sp_avg_lzc),'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 1])
hold on
[r,p] = corrcoef(all_ages,mean(sp_avg_lzc));
xlabel('age log2(months)')
ylabel('asleep LZC')
mylsline
%ylim([min(mean(sp_avg_lzc)) max(mean(sp_avg_lzc))])
title(sprintf('r = %1.3f, p = %1.3f, scales = %i-%i',r(2,1),p(2,1),scales(1),scales(end)),'fontsize',18)
makefigpretty
figname = sprintf('age_sp_LZC_%i_%i',scales(1),scales(end));
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

% plot R^2 as a function of time scale
Rwk = nan(1,size(wk_avg_lzc,1));
Rsp = nan(1,size(sp_avg_lzc,1));
for iscale = 1:size(wk_avg_lzc,1)
    rwk = corrcoef(wk_avg_lzc(iscale,:),all_ages,'rows','complete');
    rsp = corrcoef(sp_avg_lzc(iscale,:),all_ages,'rows','complete');
    Rwk(iscale) = rwk(2,1)^2;
    Rsp(iscale) = rsp(2,1)^2;
end

figure('Color','w')
plot([Rwk; Rsp]','linewidth',2)
xlabel('Time scale')
ylabel('R^{2}')
legend('Wake','Sleep')
legend boxoff
title('Variance of LZC explained by age')
makefigpretty
figname = 'R2LZCage';
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

%% Do scatter plots of age and gMLZ

figure('Color','w')
scatter(log2(all_ages),mean(wk_avg_lzc),'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 1])
[r,p] = corrcoef(all_ages,mean(wk_avg_lzc));
xlabel('age log2(months)')
ylabel('awake LZC')
mylsline
%ylim([min(mean(wk_avg_lzc)) max(mean(wk_avg_lzc))])
title(sprintf('r = %1.3f, p = %1.3f, scales = %i-%i',r(2,1),p(2,1),scales(1),scales(end)),'fontsize',18)
makefigpretty
figname = sprintf('age_wk_lzc_%i_%i',scales(1),scales(end));
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

figure('Color','w')
scatter(log2(all_ages),mean(sp_avg_lzc),'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 1])
hold on
[r,p] = corrcoef(all_ages,mean(sp_avg_lzc));
xlabel('age log2(months)')
ylabel('asleep LZC')
mylsline
%ylim([min(mean(sp_avg_lzc)) max(mean(sp_avg_lzc))])
title(sprintf('r = %1.3f, p = %1.3f, scales = %i-%i',r(2,1),p(2,1),scales(1),scales(end)),'fontsize',18)
makefigpretty
figname = sprintf('age_sp_lzc_%i_%i',scales(1),scales(end));
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

% plot R^2 as a function of time scale
Rwk = nan(1,size(wk_avg_lzc,1));
Rsp = nan(1,size(sp_avg_lzc,1));
for iscale = 1:size(wk_avg_lzc,1)
    rwk = corrcoef(wk_avg_lzc(iscale,:),all_ages,'rows','complete');
    rsp = corrcoef(sp_avg_lzc(iscale,:),all_ages,'rows','complete');
    Rwk(iscale) = rwk(2,1)^2;
    Rsp(iscale) = rsp(2,1)^2;
end

figure('Color','w')
plot([Rwk; Rsp]','linewidth',2)
xlabel('Time scale')
ylabel('R^{2}')
legend('Wake','Sleep')
legend boxoff
title('Variance of LZC explained by age')
makefigpretty
figname = 'R2LZCage';
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))


% %% DO SEPERATE PLOTS BY GENOTYPE 
% 
% %% plot average MSE scales - deletion genotype only
% 
% figure('Color','w')
% hold on
% [~,~,ci] = ttest(wk_scl(:,logical(allgenotype))');
% plot(mean(wk_scl(:,logical(allgenotype)),2),'r','linewidth',4)
% patch( [1:length(scales) fliplr(1:length(scales))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
%     1.0,'EdgeColor','none')
% [~,~,ci] = ttest(sp_scl(:,logical(allgenotype))');
% plot(mean(sp_scl(:,logical(allgenotype)),2),'b','linewidth',4)
% patch( [1:length(scales) fliplr(1:length(scales))] , [ci(1,:) fliplr(ci(2,:))], 'b','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('wake','wake 95% CI','sleep','95% sleep','AutoUpdate','off','FontSize',14,'location','southeast')
% xlim([1 length(scales)])
% title('Deletion genotype MSE values over all time scales','fontsize',18)
% xlabel('time scale')
% ylabel('sample entropy')
% %axis([1 20 0 2.5])
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'Deletion_MSE_scales'))
% 
% figure('Color','w')
% hold on
% wsdelta = mean(wk_scl(:,logical(allgenotype)),2) - mean(sp_scl(:,logical(allgenotype)),2);
% [~,~,ci] = ttest(wk_scl(:,logical(allgenotype)),sp_scl(:,logical(allgenotype)),'dim',2); ci = ci';
% plot(wsdelta ,'m','linewidth',4)
% patch( [1:length(scales) fliplr(1:length(scales))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',14,'location','northeast')
% plot(1:length(wsdelta),zeros(1,length(wsdelta)),'-.k')
% xlim([1 length(scales)])
% %ylim([-0.1 0.5])
% title('Deletion genotype MSE difference over all time scales','fontsize',18)
% xlabel('time scale')
% ylabel('sample entropy')
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'Deletion_MSE_scales_delta'))
% 
% [h,pval,ci,stat] = ttest(wk_scl(:,logical(allgenotype)),sp_scl(:,logical(allgenotype)),'dim',2);
% Q = mafdr(pval,'BHFDR',true); % Bejamini Hochberg FDR
% ngood = sum(Q < 0.05);
% fprintf('\n%i deletion time scales have significant increase in MSE after FDR correction\n',ngood)
% 
% 
% 
% %% plot average MSE scales - nondeletion genotype only
% 
% figure('Color','w')
% hold on
% [~,~,ci] = ttest(wk_scl(:,logical(~allgenotype))');
% plot(mean(wk_scl(:,logical(~allgenotype)),2),'r','linewidth',4)
% patch( [1:length(scales) fliplr(1:length(scales))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
%     1.0,'EdgeColor','none')
% [~,~,ci] = ttest(sp_scl(:,logical(~allgenotype))');
% plot(mean(sp_scl(:,logical(~allgenotype)),2),'b','linewidth',4)
% patch( [1:length(scales) fliplr(1:length(scales))] , [ci(1,:) fliplr(ci(2,:))], 'b','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('wake','wake 95% CI','sleep','95% sleep','AutoUpdate','off','FontSize',14,'location','southeast')
% xlim([1 length(scales)])
% title('Non-Deletion genotype MSE values over all time scales','fontsize',18)
% xlabel('time scale')
% ylabel('sample entropy')
% %axis([1 20 0 2.5])
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'Non-Deletion_MSE_scales'))
% 
% figure('Color','w')
% hold on
% wsdelta = mean(wk_scl(:,logical(~allgenotype)),2) - mean(sp_scl(:,logical(~allgenotype)),2);
% [~,~,ci] = ttest(wk_scl(:,logical(~allgenotype)),sp_scl(:,logical(~allgenotype)),'dim',2); ci = ci';
% plot(wsdelta ,'m','linewidth',4)
% patch( [1:length(scales) fliplr(1:length(scales))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',14,'location','northeast')
% plot(1:length(wsdelta),zeros(1,length(wsdelta)),'-.k')
% xlim([1 length(scales)])
% %ylim([-0.1 0.5])
% title('Non-Deletion genotype MSE difference over all time scales','fontsize',18)
% xlabel('time scale')
% ylabel('sample entropy')
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'Non-Deletion_MSE_scales_delta'))
% 
% [h,pval,ci,stat] = ttest(wk_scl(:,logical(~allgenotype)),sp_scl(:,logical(~allgenotype)),'dim',2);
% Q = mafdr(pval,'BHFDR',true); % Bejamini Hochberg FDR
% ngood = sum(Q < 0.05);
% fprintf('\n%i non-deletion time scales have significant increase in MSE after FDR correction\n',ngood)
% 
% %% plot Awake MSE separately for each genotype
% 
% figure('Color','w')
% hold on
% [~,~,ci] = ttest(wk_scl(:,logical(allgenotype))');
% plot(mean(wk_scl(:,logical(allgenotype)),2),'color',[0.5 1 0.5],'linewidth',4)
% patch( [1:length(scales) fliplr(1:length(scales))] , [ci(1,:) fliplr(ci(2,:))], [0.5 1 0.5],'facealpha', ...
%     1.0,'EdgeColor','none')
% [~,~,ci] = ttest(wk_scl(:,logical(~allgenotype))');
% plot(mean(wk_scl(:,logical(~allgenotype)),2),'color',[0 0.5 0.5],'linewidth',4)
% patch( [1:length(scales) fliplr(1:length(scales))] , [ci(1,:) fliplr(ci(2,:))],[0 0.5 0.5],'facealpha', ...
%     1.0,'EdgeColor','none')
% legend('DEL','DEL 95% CI','Non-DEL','Non-DEL 95% sleep','AutoUpdate','off','FontSize',14,'location','southeast')
% xlim([1 length(scales)])
% title('Awake MSE values over all time scales','fontsize',18)
% xlabel('time scale')
% ylabel('sample entropy')
% %axis([1 20 0 2.5])
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'Genotype_awake_MSE_scales'))
% 
% figure('Color','w')
% hold on
% wsdelta = mean(wk_scl(:,logical(allgenotype)),2) - mean(wk_scl(:,logical(~allgenotype)),2);
% [~,~,ci] = ttest2(wk_scl(:,logical(allgenotype)),wk_scl(:,logical(~allgenotype)),'dim',2); ci = ci';
% plot(wsdelta ,'k','linewidth',4)
% patch( [1:length(scales) fliplr(1:length(scales))] , [ci(1,:) fliplr(ci(2,:))], 'k','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('DEL-Non-DEL','DEL-Non-DEL 95% CI','AutoUpdate','off','FontSize',14,'location','southeast')
% plot(1:length(wsdelta),zeros(1,length(wsdelta)),'-.k')
% xlim([1 length(scales)])
% title('DEL - Non-DEL Awake MSE difference over all time scales','fontsize',18)
% xlabel('time scale')
% ylabel('sample entropy')
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'Genotype_awake_MSE_scales_delta'))
% 
% [h,pval,ci,stat] = ttest2(wk_scl(:,logical(allgenotype)),wk_scl(:,logical(~allgenotype)),'dim',2);
% Q = mafdr(pval,'BHFDR',true); % Bejamini Hochberg FDR
% ngood = sum(Q < 0.05);
% fprintf('Genotype effect awake, \n%i time scales have significant increase in MSE after FDR correction\n',ngood)
% 
% %% plot sleep MSE separately for each genotype
% 
% figure('Color','w')
% hold on
% [~,~,ci] = ttest(sp_scl(:,logical(allgenotype))');
% plot(mean(sp_scl(:,logical(allgenotype)),2),'color',[0.5 1 0.5],'linewidth',4)
% patch( [1:length(scales) fliplr(1:length(scales))] , [ci(1,:) fliplr(ci(2,:))], [0.5 1 0.5],'facealpha', ...
%     1.0,'EdgeColor','none')
% [~,~,ci] = ttest(sp_scl(:,logical(~allgenotype))');
% plot(mean(sp_scl(:,logical(~allgenotype)),2),'color',[0 0.5 0.5],'linewidth',4)
% patch( [1:length(scales) fliplr(1:length(scales))] , [ci(1,:) fliplr(ci(2,:))],[0 0.5 0.5],'facealpha', ...
%     1.0,'EdgeColor','none')
% legend('DEL','DEL 95% CI','Non-DEL','Non-DEL 95% sleep','AutoUpdate','off','FontSize',14,'location','southeast')
% xlim([1 length(scales)])
% title('sleep MSE values over all time scales','fontsize',18)
% xlabel('time scale')
% ylabel('sample entropy')
% %axis([1 20 0 2.5])
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'Genotype_sleep_MSE_scales'))
% 
% figure('Color','w')
% hold on
% wsdelta = mean(sp_scl(:,logical(allgenotype)),2) - mean(sp_scl(:,logical(~allgenotype)),2);
% [~,~,ci] = ttest2(sp_scl(:,logical(allgenotype)),sp_scl(:,logical(~allgenotype)),'dim',2); ci = ci';
% plot(wsdelta ,'k','linewidth',4)
% patch( [1:length(scales) fliplr(1:length(scales))] , [ci(1,:) fliplr(ci(2,:))], 'k','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('DEL-Non-DEL','DEL-Non-DEL 95% CI','AutoUpdate','off','FontSize',14,'location','southeast')
% plot(1:length(wsdelta),zeros(1,length(wsdelta)),'-.k')
% xlim([1 length(scales)])
% title('DEL - Non-DEL sleep MSE difference over all time scales','fontsize',18)
% xlabel('time scale')
% ylabel('sample entropy')
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'Genotype_sleep_MSE_scales_delta'))
% 
% [h,pval,ci,stat] = ttest2(sp_scl(:,logical(allgenotype)),sp_scl(:,logical(~allgenotype)),'dim',2);
% Q = mafdr(pval,'BHFDR',true); % Bejamini Hochberg FDR
% ngood = sum(Q < 0.05);
% fprintf('\nGenotype effect sleep: %i time scales have significant increase in MSE after FDR correction\n',ngood)
% 
% 
% %% plot average gMLZ scales DELETION
% 
% figure('Color','w')
% h=subplot(1,1,1);
% hold on
% [~,~,ci] = ttest(wk_scl_lzc(:,logical(allgenotype))');
% plot(log2(LZC_foi(scales)),mean(wk_scl_lzc(:,logical(allgenotype)),2),'r','linewidth',4)
% patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
%     1.0,'EdgeColor','none')
% [~,~,ci] = ttest(sp_scl_lzc(:,logical(allgenotype))');
% plot(log2(LZC_foi(scales)),mean(sp_scl_lzc(:,logical(allgenotype)),2),'b','linewidth',4)
% patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], 'b','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('wake','wake 95% CI','sleep','95% sleep','AutoUpdate','off','FontSize',14,'location','southeast')
% xlim([1 length(scales)])
% title('Deletion gMLZ values over all time scales','fontsize',18)
% xlabel('Center frequency (Hz)')
% ylabel('LZ complexity')
% set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
% xlim([0 5])
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'Deletion_gMLZ_scales'))
% 
% figure('Color','w')
% h=subplot(1,1,1);
% hold on
% wsdelta = mean(wk_scl_lzc(:,logical(allgenotype)),2) - mean(sp_scl_lzc(:,logical(allgenotype)),2);
% [~,~,ci] = ttest(wk_scl_lzc(:,logical(allgenotype)),sp_scl_lzc(:,logical(allgenotype)),'dim',2); ci = ci';
% plot(log2(LZC_foi(scales)),wsdelta ,'m','linewidth',4)
% patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',14,'location','southeast')
% plot(log2(LZC_foi(scales)),zeros(1,length(wsdelta)),'-.k')
% xlim([0 5])
% %ylim([-0.04 0.04])
% title('Deletion gMLZ difference over all time scales','fontsize',18)
% xlabel('Center frequency (Hz)')
% ylabel('LZ complexity')
% set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'Deletion_gMLZ_scales_delta'))
% 
% [h,pval,ci,stat] = ttest(wk_scl_lzc(:,logical(allgenotype)),sp_scl_lzc(:,logical(allgenotype)),'dim',2);
% Q = mafdr(pval,'BHFDR',true); % Bejamini Hochberg FDR
% ngood = sum(Q < 0.05);
% fprintf('\n%i time scales have significant increase in DELETION gMLZ after FDR correction\n',ngood)
% 
% %% plot average gMLZ scales NONDELETION
% 
% figure('Color','w')
% h=subplot(1,1,1);
% hold on
% [~,~,ci] = ttest(wk_scl_lzc(:,logical(~allgenotype))');
% plot(log2(LZC_foi(scales)),mean(wk_scl_lzc(:,logical(~allgenotype)),2),'r','linewidth',4)
% patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
%     1.0,'EdgeColor','none')
% [~,~,ci] = ttest(sp_scl_lzc(:,logical(~allgenotype))');
% plot(log2(LZC_foi(scales)),mean(sp_scl_lzc(:,logical(~allgenotype)),2),'b','linewidth',4)
% patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], 'b','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('wake','wake 95% CI','sleep','95% sleep','AutoUpdate','off','FontSize',14,'location','southeast')
% xlim([1 length(scales)])
% title('NONDELETION gMLZ values over all time scales','fontsize',18)
% xlabel('Center frequency (Hz)')
% ylabel('LZ complexity')
% set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
% xlim([0 5])
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'NONDELETION_gMLZ_scales'))
% 
% figure('Color','w')
% h=subplot(1,1,1);
% hold on
% wsdelta = mean(wk_scl_lzc(:,logical(~allgenotype)),2) - mean(sp_scl_lzc(:,logical(~allgenotype)),2);
% [~,~,ci] = ttest(wk_scl_lzc(:,logical(~allgenotype)),sp_scl_lzc(:,logical(~allgenotype)),'dim',2); ci = ci';
% plot(log2(LZC_foi(scales)),wsdelta ,'m','linewidth',4)
% patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',14,'location','southeast')
% plot(log2(LZC_foi(scales)),zeros(1,length(wsdelta)),'-.k')
% xlim([0 5])
% %ylim([-0.04 0.04])
% title('NONDELETION gMLZ difference over all time scales','fontsize',18)
% xlabel('Center frequency (Hz)')
% ylabel('LZ complexity')
% set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'NONDELETION_gMLZ_scales_delta'))
% 
% [h,pval,ci,stat] = ttest(wk_scl_lzc(:,logical(~allgenotype)),sp_scl_lzc(:,logical(~allgenotype)),'dim',2);
% Q = mafdr(pval,'BHFDR',true); % Bejamini Hochberg FDR
% ngood = sum(Q < 0.05);
% fprintf('\n%i time scales have significant increase in NONDELETION gMLZ after FDR correction\n',ngood)
% 
% %% plot Awake gMLZ separately for each genotype
% 
% figure('Color','w')
% hold on
% [~,~,ci] = ttest(wk_scl_lzc(:,logical(allgenotype))');
% plot(log2(LZC_foi(scales)),mean(wk_scl_lzc(:,logical(allgenotype)),2),'color',[0.5 1 0.5],'linewidth',4)
% patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], [0.5 1 0.5],'facealpha', ...
%     1.0,'EdgeColor','none')
% [~,~,ci] = ttest(wk_scl_lzc(:,logical(~allgenotype))');
% plot(log2(LZC_foi(scales)),mean(wk_scl_lzc(:,logical(~allgenotype)),2),'color',[0 0.5 0.5],'linewidth',4)
% patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))],[0 0.5 0.5],'facealpha', ...
%     1.0,'EdgeColor','none')
% legend('DEL','DEL 95% CI','Non-DEL','Non-DEL 95% sleep','AutoUpdate','off','FontSize',14,'location','southeast')
% xlim([0 5])
% title('Awake gMLZ values over all time scales','fontsize',18)
% xlabel('time scale')
% ylabel('LZ complexity')
% %axis([0 5 0 0.5])
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'Genotype_awake_gMLZ_scales'))
% 
% figure('Color','w')
% hold on
% wsdelta = mean(wk_scl_lzc(:,logical(allgenotype)),2) - mean(wk_scl_lzc(:,logical(~allgenotype)),2);
% [~,~,ci] = ttest2(wk_scl_lzc(:,logical(allgenotype)),wk_scl_lzc(:,logical(~allgenotype)),'dim',2); ci = ci';
% plot(log2(LZC_foi),wsdelta ,'k','linewidth',4)
% patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], 'k','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('DEL-Non-DEL','DEL-Non-DEL 95% CI','AutoUpdate','off','FontSize',14,'location','southeast')
% plot(1:length(wsdelta),zeros(1,length(wsdelta)),'-.k')
% xlim([1 5])
% title('DEL - Non-DEL Awake gMLZ difference over all time scales','fontsize',18)
% xlabel('time scale')
% ylabel('LZ complexity')
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'Genotype_awake_gMLZ_scales_delta'))
% 
% [h,pval,ci,stat] = ttest2(wk_scl_lzc(:,logical(allgenotype)),wk_scl_lzc(:,logical(~allgenotype)),'dim',2);
% Q = mafdr(pval,'BHFDR',true); % Bejamini Hochberg FDR
% ngood = sum(Q < 0.05);
% fprintf('\nGenotype effect awake: %i time scales have significant increase in gMLZ after FDR correction\n',ngood)
% 
% %% plot sleep gMLZ separately for each genotype
% 
% figure('Color','w')
% hold on
% [~,~,ci] = ttest(sp_scl_lzc(:,logical(allgenotype))');
% plot(log2(LZC_foi(scales)),mean(sp_scl_lzc(:,logical(allgenotype)),2),'color',[0.5 1 0.5],'linewidth',4)
% patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], [0.5 1 0.5],'facealpha', ...
%     1.0,'EdgeColor','none')
% [~,~,ci] = ttest(sp_scl_lzc(:,logical(~allgenotype))');
% plot(log2(LZC_foi(scales)),mean(sp_scl_lzc(:,logical(~allgenotype)),2),'color',[0 0.5 0.5],'linewidth',4)
% patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))],[0 0.5 0.5],'facealpha', ...
%     1.0,'EdgeColor','none')
% legend('DEL','DEL 95% CI','Non-DEL','Non-DEL 95% sleep','AutoUpdate','off','FontSize',14,'location','southeast')
% title('sleep gMLZ values over all time scales','fontsize',18)
% xlabel('time scale')
% ylabel('LZ complexity')
% %axis([0 5 0 0.5])
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'Genotype_sleep_gMLZ_scales'))
% 
% figure('Color','w')
% hold on
% wsdelta = mean(sp_scl_lzc(:,logical(allgenotype)),2) - mean(sp_scl_lzc(:,logical(~allgenotype)),2);
% [~,~,ci] = ttest2(sp_scl_lzc(:,logical(allgenotype)),sp_scl_lzc(:,logical(~allgenotype)),'dim',2); ci = ci';
% plot(log2(LZC_foi),wsdelta ,'k','linewidth',4)
% patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], 'k','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('DEL-Non-DEL','DEL-Non-DEL 95% CI','AutoUpdate','off','FontSize',14,'location','southeast')
% plot(1:length(wsdelta),zeros(1,length(wsdelta)),'-.k')
% xlim([0 5])
% title('DEL - Non-DEL sleep gMLZ difference over all time scales','fontsize',18)
% xlabel('time scale')
% ylabel('LZ complexity')
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'Genotype_sleep_gMLZ_scales_delta'))
% 
% [h,pval,ci,stat] = ttest2(sp_scl_lzc(:,logical(allgenotype)),sp_scl_lzc(:,logical(~allgenotype)),'dim',2);
% Q = mafdr(pval,'BHFDR',true); % Bejamini Hochberg FDR
% ngood = sum(Q < 0.05);
% fprintf('\n Genotype effect sleep: %i time scales have significant increase in gMLZ after FDR correction\n',ngood)
% 
% %% plot TD figure: average MSE scales
% 
% % Find the appropriate age match
% keep = age_td > mean(all_ages)-std(all_ages)*2 & age_td < ...
%     mean(all_ages)+std(all_ages)*2 & age_td > 12 & age_td < 12*18;
% 
% figure('Color','w')
% hold on
% [~,~,ci] = ttest(td_scl(:,keep)');
% plot(mean(td_scl(:,keep),2),'color',[0 0.5 0.5],'linewidth',4)
% patch( [1:length(scales) fliplr(1:length(scales))] , [ci(1,:) fliplr(ci(2,:))],[0 0.5 0.5],'facealpha', ...
%     1.0,'EdgeColor','none')
% legend('TD','TD 95% CI','AutoUpdate','off','FontSize',14,'location','northeast')
% xlim([1 length(scales)])
% title('MSE values over all time scales (no regression)','fontsize',18)
% xlabel('time scale')
% ylabel('sample entropy')
% %axis([1 20 -0.06 0])
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'MSE_scales_no_reg'))
% 
% fprintf('\n TD sample, n = %i\n',sum(keep))
% 
% [h,pval,ci,stat] = ttest(td_scl(:,keep));
% Q = mafdr(pval,'BHFDR',true); % Bejamini Hochberg FDR
% ngood = sum(Q < 0.05);
% fprintf('\n%i TD time scales have significant increase in MSE after FDR correction\n',ngood)


% 
% 
% %% Regress out effects of log2(age) and genotype on MSE
% 
% B_wk = nan(3,size(wk_all,1),size(wk_all,2));
% R_wk = nan(size(wk_all,3),size(wk_all,1),size(wk_all,2));
% 
% B_sp = nan(3,size(wk_all,1),size(wk_all,2));
% R_sp = nan(size(wk_all,3),size(wk_all,1),size(wk_all,2));
% 
% for ichan = 1:size(wk_all,1)
%     for jscale = 1:size(wk_all,2)
%         [b,~,r] = regress(squeeze(wk_all(ichan,jscale,:)),[log2(all_ages)' allgenotype' ones(1,size(wk_all,3))']);
%         B_wk(:,ichan,jscale) = b;
%         R_wk(:,ichan,jscale) = r;
%     end
% end
% 
% for ichan = 1:size(sp_all,1)
%     for jscale = 1:size(sp_all,2)
%         [b,~,r] = regress(squeeze(sp_all(ichan,jscale,:)),[log2(all_ages)' allgenotype' ones(1,size(sp_all,3))']);
%         B_sp(:,ichan,jscale) = b;
%         R_sp(:,ichan,jscale) = r;
%     end
% end
% 
% % overwrite variables with residual + y-intercept
% 
% wk_all = R_wk+B_wk(end,:,:); % regress out effect of age and genotype (residuals + y-intercept)
% sp_all = R_sp+B_sp(end,:,:); % regress out effect of age and genotype (residuals + y-intercept)

% % swap dimesions so it's channels x scales x subjects
% wk_all = permute(wk_all,[2 3 1]);
% sp_all = permute(sp_all,[2 3 1]);
% 
% wk_scl_lzc = squeeze(mean(wk_lzc_all));
% sp_scl_lzc = squeeze(mean(sp_lzc_all));
% 
% wk_avg_lzc = squeeze(mean(wk_all,2));
% sp_avg_lzc = squeeze(mean(sp_all,2));

% %% Regress out effects of log2(age) and genotype on gMLZ
% 
% B_wk_lzc = nan(3,size(wk_lzc_all,1),size(wk_lzc_all,2));
% R_wk_lzc = nan(size(wk_lzc_all,3),size(wk_lzc_all,1),size(wk_lzc_all,2));
% 
% B_sp_lzc = nan(3,size(wk_lzc_all,1),size(wk_lzc_all,2));
% R_sp_lzc = nan(size(wk_lzc_all,3),size(wk_lzc_all,1),size(wk_lzc_all,2));
% 
% for ichan = 1:size(wk_lzc_all,1)
%     for jscale = 1:size(wk_lzc_all,2)
%         [b,~,r] = regress(squeeze(wk_lzc_all(ichan,jscale,:)),[log2(all_ages)' allgenotype' ones(1,size(wk_lzc_all,3))']);
%         B_wk_lzc(:,ichan,jscale) = b;
%         R_wk_lzc(:,ichan,jscale) = r;
%     end
% end
% 
% for ichan = 1:size(sp_lzc_all,1)
%     for jscale = 1:size(sp_lzc_all,2)
%         [b,~,r] = regress(squeeze(sp_lzc_all(ichan,jscale,:)),[log2(all_ages)' allgenotype' ones(1,size(sp_lzc_all,3))']);
%         B_sp_lzc(:,ichan,jscale) = b;
%         R_sp_lzc(:,ichan,jscale) = r;
%     end
% end
% 
% % overwrite variables with residual + y-intercept
% 
% wk_lzc_all = R_wk_lzc+B_wk_lzc(end,:,:); % regress out effect of age and genotype (residuals + y-intercept)
% sp_lzc_all = R_sp_lzc+B_sp_lzc(end,:,:); % regress out effect of age and genotype (residuals + y-intercept

% % swap dimesions so it's channels x scales x subjects
% wk_lzc_all = permute(wk_lzc_all,[2 3 1]);
% sp_lzc_all = permute(sp_lzc_all,[2 3 1]);
% 
% wk_scl = squeeze(mean(wk_all));
% sp_scl = squeeze(mean(sp_all));
% 
% wk_avg = squeeze(mean(wk_all,2));
% sp_avg = squeeze(mean(sp_all,2));

% %% plot average MSE scales
% 
% figure('Color','w')
% hold on
% [~,~,ci] = ttest(wk_scl');
% plot(mean(wk_scl,2),'r','linewidth',4)
% patch( [1:length(scales) fliplr(1:length(scales))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
%     1.0,'EdgeColor','none')
% [~,~,ci] = ttest(sp_scl');
% plot(mean(sp_scl,2),'b','linewidth',4)
% patch( [1:length(scales) fliplr(1:length(scales))] , [ci(1,:) fliplr(ci(2,:))], 'b','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('wake','wake 95% CI','sleep','95% sleep','AutoUpdate','off','FontSize',14,'location','southeast')
% xlim([1 length(scales)])
% title('MSE values over all time scales','fontsize',18)
% xlabel('time scale')
% ylabel('sample entropy')
% %axis([1 20 0 2.5])
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'MSE_scales'))
% 
% figure('Color','w')
% hold on
% wsdelta = mean(wk_scl,2) - mean(sp_scl,2);
% [~,~,ci] = ttest(wk_scl,sp_scl,'dim',2); ci = ci';
% plot(wsdelta ,'m','linewidth',4)
% patch( [1:length(scales) fliplr(1:length(scales))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',14,'location','northeast')
% plot(1:length(wsdelta),zeros(1,length(wsdelta)),'-.k')
% xlim([1 length(scales)])
% %ylim([-0.1 0.5])
% title('MSE difference over all time scales','fontsize',18)
% xlabel('time scale')
% ylabel('sample entropy')
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'MSE_scales_delta'))
% 
% remove columns that have nans in either condition (sleep or wake)

both = wk_avg + sp_avg;
nanidx = sum(isnan(both))>0;

wk_avg(:,nanidx) = [];
sp_avg(:,nanidx) = [];

assert(size(wk_avg,2)==size(sp_avg,2),'Conditions have different sizes!')

pvals = nan(1,19);
tvals = nan(1,19);

for ichan = 1:19
    [~,p,~,stat] = ttest(wk_avg(ichan,:),sp_avg(ichan,:));
    pvals(ichan) = p;
    tvals(ichan) = stat.tstat;
end

Q = mafdr(pvals,'BHFDR',true); % Bejamini Hochberg FDR

nsig = sum(Q<0.05);

fprintf('%i MSE channels significant after FDR correction, alpha = 0.05\n',nsig)
% 
% 
% %% plot average gMLZ scales
% 
% figure('Color','w')
% h=subplot(1,1,1);
% hold on
% [~,~,ci] = ttest(wk_scl_lzc');
% plot(log2(LZC_foi(scales)),mean(wk_scl_lzc,2),'r','linewidth',4)
% patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
%     1.0,'EdgeColor','none')
% [~,~,ci] = ttest(sp_scl_lzc');
% plot(log2(LZC_foi(scales)),mean(sp_scl_lzc,2),'b','linewidth',4)
% patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], 'b','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('wake','wake 95% CI','sleep','95% sleep','AutoUpdate','off','FontSize',14,'location','southeast')
% xlim([1 length(scales)])
% title('gMLZ values over all time scales','fontsize',18)
% xlabel('Center frequency (Hz)')
% ylabel('LZ complexity')
% set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
% xlim([0 5])
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_scales'))
% 
% figure('Color','w')
% h=subplot(1,1,1);
% hold on
% wsdelta = mean(wk_scl_lzc,2) - mean(sp_scl_lzc,2);
% [~,~,ci] = ttest(wk_scl_lzc,sp_scl_lzc,'dim',2); ci = ci';
% plot(log2(LZC_foi(scales)),wsdelta ,'m','linewidth',4)
% patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',14,'location','southeast')
% plot(log2(LZC_foi(scales)),zeros(1,length(wsdelta)),'-.k')
% xlim([0 5])
% %ylim([-0.04 0.04])
% title('gMLZ difference over all time scales','fontsize',18)
% xlabel('Center frequency (Hz)')
% ylabel('LZ complexity')
% set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_scales_delta'))
% 
% [h,pval,ci,stat] = ttest(wk_scl_lzc,sp_scl_lzc,'dim',2);
% Q = mafdr(pval,'BHFDR',true); % Bejamini Hochberg FDR
% ngood = sum(Q < 0.05);
% fprintf('\n%i time scales have significant increase in gMLZ after FDR correction\n',ngood)
% 
% remove columns that have nans in either condition (sleep or wake)

both = wk_avg_lzc + sp_avg_lzc;
nanidx = sum(isnan(both))>0;

wk_avg_lzc(:,nanidx) = [];
sp_avg_lzc(:,nanidx) = [];

assert(size(wk_avg_lzc,2)==size(sp_avg_lzc,2),'Conditions have different sizes!')

pvals = nan(1,19);
tvals = nan(1,19);

for ichan = 1:19
    [~,p,~,stat] = ttest(wk_avg_lzc(ichan,:),sp_avg_lzc(ichan,:));
    pvals(ichan) = p;
    tvals(ichan) = stat.tstat;
end

Q = mafdr(pvals,'BHFDR',true); % Bejamini Hochberg FDR

nsig = sum(Q<0.05);

fprintf('%i LZC channels significant  after FDR correction, alpha = 0.05\n',nsig)
% 
% %keyboard

%% Extract delta power 

% absolute power
WP = nan(19,length(foi),length(wk_pow)); 
SP = nan(19,length(foi),length(sp_pow));
% relative power
WPr = nan(19,length(foi),length(wk_pow));
SPr = nan(19,length(foi),length(sp_pow));


for i = 1:length(wk_pow)
    if isstruct(wk_pow{i})
        WP(:,:,i) = wk_pow{i}.pow;
        WPv(:,:,i) = wk_pow{i}.pow_var;
        WPr(:,:,i) = wk_pow{i}.pow_rel;
    end
end


for i = 1:length(sp_pow)
    if isstruct(sp_pow{i})
        SP(:,:,i) = sp_pow{i}.pow;
        SPv(:,:,i) = sp_pow{i}.pow_var;
        SPr(:,:,i) = sp_pow{i}.pow_rel;
    end
end

f1=1; f2 = 4;
idf = foi> f1 & foi<f2;

delta_wk = nan(size(WP,1),size(WP,3));
delta_sp = nan(size(SP,1),size(SP,3));
delta_rel_wk = nan(size(WPr,1),size(WPr,3));
delta_rel_sp = nan(size(SPr,1),size(SPr,3));

for ich = 1:size(WP,1)
    for isb = 1:size(WP,3)
        delta_wk(ich,isb) = trapz(foi(idf),WP(ich,idf,isb));
        delta_sp(ich,isb) = trapz(foi(idf),SP(ich,idf,isb));
    end
end

for ich = 1:size(WPr,1)
    for isb = 1:size(WPr,3)
        delta_rel_wk(ich,isb) = trapz(foi(idf),WPr(ich,idf,isb));
        delta_rel_sp(ich,isb) = trapz(foi(idf),SPr(ich,idf,isb));
    end
end

for ich = 1:size(WPv,1)
    for isb = 1:size(WPv,3)
        delta_var_wk(ich,isb) = trapz(foi(idf),WPv(ich,idf,isb));
        delta_var_sp(ich,isb) = trapz(foi(idf),SPv(ich,idf,isb));
    end
end

%% PSD plots absolute power

figure('Color','w')
h=subplot(1,1,1);
[~,~,wkci] = ttest( log10(squeeze(nanmean(WP)))' );
[~,~,spci] = ttest( log10(squeeze(nanmean(SP)))' );
foi_hd = 2.^[0:0.01:5]; % interpolation resolution for plotting
wk = log10(squeeze(nanmean(WP)))';
sp = log10(squeeze(nanmean(SP)))';
wkclr = [1 0 0];
spclr = [0 0 1];
lbl1 = 'Wake';
lbl2 = 'Sleep';
plot(log2(foi_hd),interp1(foi,nanmean(wk)',foi_hd,'spline'),'Color',wkclr,'LineWidth',4)
hold on
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,wkci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,wkci(2,:),foi_hd,'spline'))],wkclr,'facealpha',1.0,'EdgeColor','none')
plot(log2(foi_hd),interp1(foi,nanmean(sp)',foi_hd,'spline'),'Color',spclr,'LineWidth',4)
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,spci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,spci(2,:),foi_hd,'spline'))],spclr,'facealpha',1.0,'EdgeColor','none')
legend({lbl1,sprintf('%s%s',lbl1,' 95% CI'),lbl2,sprintf('%s%s',lbl2,' 95% CI')},...
    'AutoUpdate','off','FontSize',14)
xlim([0 5])
xlabel('Frequency (Hz)','FontSize',12)
ylabel('Power log_1_0(\muV^{2}/Hz)','FontSize',14)
title('Absolute Power','FontSize',20)
set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
box off
%axis square
legend boxoff 
ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
    'LineStyle','none','FontSize',14);
set(gca,'linewidth',3)
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 22)
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 22)
set(gca, 'TickDir', 'out')
figname = sprintf('PSDs_abs_power_sp_wk');
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
makefigpretty

clear wkci spci wk sp

%% PSD difference plot

figure('Color','w')
h=subplot(1,1,1);
foi_hd = 2.^[0:0.01:5]; % interpolation resolution for plotting
wk = log10(squeeze(nanmean(WP)))';
sp = log10(squeeze(nanmean(SP)))';
wsd = mean(wk)-mean(sp);
[~,~,ci] = ttest(wk,sp);
lbl1 = 'Wake - Sleep';
plot(log2(foi_hd),interp1(foi,wsd',foi_hd,'spline'),'m','LineWidth',4)
hold on
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,ci(2,:),foi_hd,'spline'))],'m','facealpha',1.0,'EdgeColor','none')
legend({lbl1,sprintf('%s%s',lbl1,' 95% CI')},'AutoUpdate','off','FontSize',14)
plot(log2(foi_hd),zeros(1,length(foi_hd)),'k--')
xlim([0 5])
xlabel('Frequency (Hz)','FontSize',12)
%ylabel('Power log_1_0(\muV^{2}/Hz)','FontSize',14)
ylabel('Power change from sleep (dB)')
title('Absolute Power','FontSize',20)
set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
box off
%axis square
legend boxoff 
ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
    'LineStyle','none','FontSize',14);
set(gca,'linewidth',3)
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 22)
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 22)
set(gca, 'TickDir', 'out')
figname = sprintf('PSDs_abs_power_DIFFERNCE');
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
makefigpretty 

%% power permutation clustering 

clst_pow = RMpermclusttest(WP,SP,Nperm);

sgcl_pos = [];
sgcl_neg = [];

if sum(clst_pow.n_pos > 0) > 0 
    sgcl_pos = clst_pow.P_val_pos < 0.05;
end

if sum(clst_pow.n_neg > 0) > 0 
    sgcl_neg = clst_pow.P_val_neg < 0.05;
end
%%
cnt = 0;
for icls = 1:length(sgcl_pos)
    if sgcl_pos(icls)
        cnt = cnt + 1;
        figure('Color','w')
        h=subplot(1,1,1);
        plot(log2(foi),sum(clst_pow.lbl_pos==icls,2),'linewidth',2)
        ylim([0 20])
        xlim([0 5])
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        xlim([0 5])
        xlabel('Frequency (Hz)')
        ylabel('Number of channels')
        makefigpretty
        title(sprintf('power postivie cluster, p = %1.3f',clst_pow.P_val_pos(icls)),'fontsize',18)
        figname = sprintf('Power pos clust spect #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        
        figure('Color','w')
        plot_topo_AS(sum(clst_pow.lbl_pos==icls)')
        title(sprintf('power pos clust, p = %1.3f',clst_pow.P_val_pos(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(WP,2)])
        mycolorbar
        figname = sprintf('Power pos clust toop #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
    end
end

cnt = 0;
for icls = 1:length(sgcl_neg)
    if sgcl_neg(icls)
        cnt = cnt + 1;
        figure('Color','w')
        h=subplot(1,1,1);
        plot(log2(foi),sum(clst_pow.lbl_neg==icls,2),'linewidth',2)
        ylim([0 20])
        xlim([0 5])
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        xlim([0 5])
        xlabel('Frequency (Hz)')
        ylabel('Number of channels')
        makefigpretty
        title(sprintf('power negative cluster, p = %1.3f',clst_pow.P_val_neg(icls)),'fontsize',18)
        figname = sprintf('Power neg clust spectrum #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        
        figure('Color','w')
        plot_topo_AS(sum(clst_pow.lbl_neg==icls)')
        title(sprintf('power neg clust, p = %1.3f',clst_pow.P_val_neg(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(WP,2)])
        mycolorbar
        figname = sprintf('Power neg clust toop #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
    end
end


%% Effect sizes, plot both in same figure (not sure if this is right way to show it though) 

WPscl = squeeze(mean(WP));
SPscl = squeeze(mean(SP));
cohend_pow = nan(1,size(WPscl,1));
for ifrq = 1:length(cohend_pow)
    cohend_pow(ifrq) = effect_size(mean(WPscl(ifrq,:),2),mean(SPscl(ifrq,:),2),...
        std(WPscl(ifrq,:),[],2),std(SPscl(ifrq,:),[],2),...
        size(WPscl(ifrq,:),2),size(SPscl(ifrq,:),2));
end
cohend_lzc = nan(1,size(wk_scl_lzc,1));
for iscl = 1:length(cohend_lzc)
    cohend_lzc(iscl) = effect_size(mean(wk_scl_lzc(iscl,:),2),mean(sp_scl_lzc(iscl,:),2),...
        std(wk_scl_lzc(iscl,:),[],2),std(sp_scl_lzc(iscl,:),[],2),...
        size(wk_scl_lzc(iscl,:),2),size(sp_scl_lzc(iscl,:),2));
end
figure
h=subplot(1,1,1);
hold on
plot(log2(LZC_foi(scales)),cohend_lzc ,'linewidth',4)
plot(log2(foi),cohend_pow,'linewidth',4)
plot(log2(foi),zeros(1,length(foi)),'k--')
title('Effect sizes, wake-sleep','fontsize',18)
legend({'gMLZ wake-sleep','Spectral power wake-sleep'},'AutoUpdate','off',...
    'FontSize',14,'location','southeast')
legend boxoff
xlabel('(Center) frequency (Hz)')
ylabel('Effect size (Cohen''s D)')
makefigpretty 
xlim([0 5])
ylim([-1 1])
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'CohenD_ALL'))

%% do same but averaged across scale rather than channel

% take average across delta band
fmax = 4;
WPavg = squeeze(mean(WP(:,foi<=fmax,:),2));
SPavg = squeeze(mean(SP(:,foi<=fmax,:),2));
wk_LZ_DELTA = squeeze(mean(wk_lzc_all(:,LZC_foi<=fmax,:),2));
sp_LZ_DELTA = squeeze(mean(sp_lzc_all(:,LZC_foi<=fmax,:),2));

cd_pow_top = nan(1,size(WPavg,1));
for ich = 1:length(cd_pow_top)
    cd_pow_top(ich) = effect_size(mean(WPavg(ich,:),2),mean(SPavg(ich,:),2),...
        std(WPavg(ich,:),[],2),std(SPavg(ich,:),[],2),...
        size(WPavg(ich,:),2),size(SPavg(ich,:),2));
end
cd_lzc_top = nan(1,size(wk_LZ_DELTA,1));
for ich = 1:length(cd_lzc_top)
    cd_lzc_top(ich) = effect_size(mean(wk_LZ_DELTA(ich,:),2),mean(sp_LZ_DELTA(ich,:),2),...
        std(wk_LZ_DELTA(ich,:),[],2),std(sp_LZ_DELTA(ich,:),[],2),...
        size(wk_LZ_DELTA(ich,:),2),size(sp_LZ_DELTA(ich,:),2));
end

figure('Color','w')
plot_topo_AS(cd_pow_top')
title('Spectral power effect sizes \delta frequencies','fontsize',18)
mycolorbar
caxis([-1 1])
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'CohenD_pow_topo'))

figure('Color','w')
plot_topo_AS(cd_lzc_top')
title('gMLZ effect sizes \delta scales','fontsize',18)
mycolorbar
caxis([-1 1])
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'CohenD_lzc_topo'))
%%
% 
% 
% %% ASBOLUE DELTA POWER TOPOS
% 
% % Delta power -log10P topoplots
% 
% figure('Color','w')
% plot_topo_AS(-log10(Qa))
% title(sprintf('Delta absolute power -log10(p), FDR corrected'),'fontsize',18)
% colormap jet
% caxis([0 7])
% mycolorbar
% figname = sprintf('delta_power_pvals_log10p');
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% %% Delta power t-stats topoplots
% 
% figure('Color','w')
% plot_topo_AS(stata.tstat)
% title(sprintf('Delta absolute power t-stats'),'fontsize',18)
% colormap jet
% caxis([0 7])
% mycolorbar
% figname = sprintf('delta_power_tstats');
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% % Delta power awake topoplot
% 
% figure('Color','w')
% plot_topo_AS(nanmean(log10(delta_wk),2))
% title(sprintf('Mean delta absolute power, wake condition'),'fontsize',18)
% caxis([0 3])
% colormap jet
% mycolorbar
% figname = sprintf('delta_power_wake');
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% % Delta power sleep topoplot
% 
% figure('Color','w')
% plot_topo_AS(nanmean(log10(delta_sp),2))
% caxis([0 3])
% mycolorbar
% title(sprintf('Mean delta absolute power, sleep condition'),'fontsize',18)
% figname = sprintf('delta_power_sleep');
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% %% Delta power wake vs sleep mean diff topoplot
% 
% figure('Color','w')
% delta = nanmean(log10(delta_sp),2) - nanmean(log10(delta_wk),2);
% plot_topo_AS(delta)
% title(sprintf('Mean delta absolute power difference, sleep - wake'),'fontsize',18)
% colormap jet
% caxis([0 0.25])
% mycolorbar
% figname = sprintf('delta_power_sleep_minus_wake');
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% %% DELTA RELATIVE POWER TOPOS
% 
% % Delta power -log10P topoplots
% 
% figure('Color','w')
% plot_topo_AS(-log10(Qr))
% title(sprintf('Delta relative power -log10(p), FDR corrected'),'fontsize',18)
% colormap jet
% caxis([0 7])
% mycolorbar
% figname = sprintf('Rel_delta_power_pvals_log10p');
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% % Delta power t-stats topoplots
% 
% figure('Color','w')
% plot_topo_AS(statr.tstat)
% title(sprintf('Delta relative power t-stats'),'fontsize',18)
% colormap jet
% caxis([0 7])
% mycolorbar
% figname = sprintf('Rel_delta_power_tstats');
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% % Delta power awake topoplot
% 
% figure('Color','w')
% plot_topo_AS(nanmean(log10(delta_rel_wk),2))
% title(sprintf('Mean delta relative power, wake condition'),'fontsize',18)
% caxis([0 3])
% colormap jet
% mycolorbar
% figname = sprintf('Rel_delta_power_wake');
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% % Delta power sleep topoplot
% 
% figure('Color','w')
% plot_topo_AS(nanmean(log10(delta_rel_sp),2))
% caxis([0 3])
% mycolorbar
% title(sprintf('Mean delta relative power, sleep condition'),'fontsize',18)
% figname = sprintf('Rel_delta_power_sleep');
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% % Delta power wake vs sleep mean diff topoplot
% 
% figure('Color','w')
% delta = nanmean(log10(delta_rel_sp),2) - nanmean(log10(delta_rel_wk),2);
% plot_topo_AS(delta)
% title(sprintf('Mean delta relative power difference, sleep - wake'),'fontsize',18)
% colormap jet
% caxis([0 0.25])
% mycolorbar
% figname = sprintf('Rel_delta_power_sleep_minus_wake');
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

% %% Transpose wk_scl and sp_scl for plotting correlations
% 
% wk_scl = wk_scl';
% sp_scl = sp_scl';

%% MSE Correlations

r_pow_mse_wk = corr(log10(mean(delta_sp))',wk_scl');
r_pow_mse_sp = corr(log10(mean(delta_sp))',sp_scl');
rpb_wk = mes(repmat(allgenotype',1,size(wk_scl,1)),wk_scl','requiv');
rpb_sp = mes(repmat(allgenotype',1,size(sp_scl,1)),sp_scl','requiv');
r_age_mse_wk = corr(all_ages',wk_scl');
r_age_mse_sp = corr(all_ages',sp_scl');
r_pow_cog_wk = corr(cog_score',wk_scl','rows','complete');
r_pow_cog_sp = corr(cog_score',sp_scl','rows','complete');
r_pow_lng_wk = corr(lng_score',wk_scl','rows','complete');
r_pow_lng_sp = corr(lng_score',sp_scl','rows','complete');

figure('Color','w')
hold on
plot(log2(nyqfoi),r_pow_mse_wk,'r','linewidth',2)
plot(log2(nyqfoi),r_pow_mse_sp,'r:','linewidth',2)
plot(log2(nyqfoi),rpb_wk.requiv,'c','linewidth',2)
plot(log2(nyqfoi),rpb_sp.requiv,'c:','linewidth',2)
plot(log2(nyqfoi),r_age_mse_wk,'b','linewidth',2)
plot(log2(nyqfoi),r_age_mse_sp,'b:','linewidth',2)
plot(log2(nyqfoi),r_pow_cog_wk,'m','linewidth',2)
plot(log2(nyqfoi),r_pow_cog_sp,'m:','linewidth',2)
plot(log2(nyqfoi),r_pow_lng_wk,'k','linewidth',2)
plot(log2(nyqfoi),r_pow_lng_sp,'k:','linewidth',2)
ylim([-1 1])
legend({'\delta power awake','\delta power sleep','Genotype awake','Genotype sleep','Age awake','Age asleep','Cognition awake', ...
    'Cognition sleep','Language awake','Langauge sleep'},'location','best')
legend boxoff
ylabel('r')
xlabel('Nyquist frequency (Hz)')
makefigpretty
xlim([log2(5) log2(100)])
ticks = 2:1:7;
xticks(ticks)
ticklabels = cell(1,length(ticks));
for i = 1:length(ticklabels)
    ticklabels{i} = num2str(2^ticks(i));
end
xticklabels(ticklabels)
xlabel('Nyquist frequency (Hz)')
ylabel('Correlation (r)')
makefigpretty
% ADD TOP X-AXIS FOR TAU
ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation',...
    'right', 'Color','none');
makefigpretty
line(log2(nyqfoi),mean(wk_scl,2),'Parent',ax2,'Color','none')
xlim([log2(5) log2(100)])
pltscl = fliplr([1 2 5 10 20]);
ticks = log2(nyqfoi(pltscl));
xticks(ticks)
ticklabels = cell(1,length(ticks));
for i = 1:length(ticklabels)
    ticklabels{i} = num2str(pltscl(i));
end
xticklabels(ticklabels)
xlabel('time scale')
set(gca,'ytick',[]);
yAX = get(gca,'YAxis');
yAX.Visible = 'off';
% SAVE FIGURE
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'MSE_correlations'))
print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,'MSE_correlations'))

%% LZC Correlations

r_pow_lzc_wk = corr(log10(mean(delta_sp))',wk_scl_lzc');
r_pow_lzc_sp = corr(log10(mean(delta_sp))',sp_scl_lzc');
rpb_lzc_wk = mes(repmat(allgenotype',1,size(wk_scl_lzc,1)),wk_scl_lzc','requiv');
rpb_lzc_sp = mes(repmat(allgenotype',1,size(sp_scl_lzc,1)),sp_scl_lzc','requiv');
r_age_lzc_wk = corr(all_ages',wk_scl_lzc');
r_age_lzc_sp = corr(all_ages',sp_scl_lzc');
r_pow_cog_wk = corr(cog_score',wk_scl_lzc','rows','complete');
r_pow_cog_sp = corr(cog_score',sp_scl_lzc','rows','complete');
r_pow_lng_wk = corr(lng_score',wk_scl_lzc','rows','complete');
r_pow_lng_sp = corr(lng_score',sp_scl_lzc','rows','complete');

figure('Color','w')
hold on
plot(log2(nyqfoi),r_pow_lzc_wk,'r','linewidth',2)
plot(log2(nyqfoi),r_pow_lzc_sp,'r:','linewidth',2)
plot(log2(nyqfoi),rpb_lzc_wk.requiv,'c','linewidth',2)
plot(log2(nyqfoi),rpb_lzc_sp.requiv,'c:','linewidth',2)
plot(log2(nyqfoi),r_age_lzc_wk,'b','linewidth',2)
plot(log2(nyqfoi),r_age_lzc_sp,'b:','linewidth',2)
plot(log2(nyqfoi),r_pow_cog_wk,'m','linewidth',2)
plot(log2(nyqfoi),r_pow_cog_sp,'m:','linewidth',2)
plot(log2(nyqfoi),r_pow_lng_wk,'k','linewidth',2)
plot(log2(nyqfoi),r_pow_lng_sp,'k:','linewidth',2)
ylim([-1 1])
legend({'\delta power awake','\delta power sleep','Genotype awake','Genotype sleep','Age awake','Age asleep','Cognition awake', ...
    'Cognition sleep','Language awake','Langauge sleep'},'location','best')
legend boxoff
ylabel('r')
xlabel('Nyquist frequency (Hz)')
makefigpretty
xlim([log2(5) log2(100)])
ticks = 2:1:7;
xticks(ticks)
ticklabels = cell(1,length(ticks));
for i = 1:length(ticklabels)
    ticklabels{i} = num2str(2^ticks(i));
end
xticklabels(ticklabels)
xlabel('Nyquist frequency (Hz)')
ylabel('Correlation (r)')
makefigpretty
% ADD TOP X-AXIS FOR TAU
ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation',...
    'right', 'Color','none');
makefigpretty
line(log2(nyqfoi),mean(wk_scl_lzc,2),'Parent',ax2,'Color','none')
xlim([log2(5) log2(100)])
pltscl = fliplr([1 2 5 10 20]);
ticks = log2(nyqfoi(pltscl));
xticks(ticks)
ticklabels = cell(1,length(ticks));
for i = 1:length(ticklabels)
    ticklabels{i} = num2str(pltscl(i));
end
xticklabels(ticklabels)
xlabel('time scale')
set(gca,'ytick',[]);
yAX = get(gca,'YAxis');
yAX.Visible = 'off';
% SAVE FIGURE
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'LZC_correlations'))
print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,'LZC_correlations'))

% %% Mean MSE vs Delta power, ASLEEP
% figure('Color','w')
% scatter(mean(sp_scl),log10(mean(delta_sp)),'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 1])
% mylsline
% ylabel('Power log_1_0(\muV^{2}/Hz)')
% xlabel(sprintf('MSE avg(%i-%i)',scales(1),scales(end)))
% [r,p] = corrcoef(log10(mean(delta_sp)),mean(sp_scl));
% %axis([1.2 1.8 2.25 4.25])
% title(sprintf('Sleep delta power vs MSE (scales %i-%i), r = %1.2f, p = %1.3f', scales(1),scales(end),r(2,1), p(2,1)),'fontsize',18)
% figname = sprintf('sleep_delta_power_vs_MSE_%i_%i',scales(1),scales(end));
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% 
% %% Mean MSE vs Delta power, AWAKE
% 
% figure('Color','w')
% scatter(mean(wk_scl),log10(mean(delta_wk)),'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 1])
% mylsline
% ylabel('Power log_1_0(\muV^{2}/Hz)')
% xlabel(sprintf('MSE avg(%i-%i)',scales(1),scales(end)))
% [r,p] = corrcoef(log10(mean(delta_wk)),mean(wk_scl));
% %axis([1.2 1.8 2.25 4.25])
% title(sprintf('Awake delta power vs MSE (scales %i-%i), r = %1.2f, p = %1.3f',scales(1),scales(end), r(2,1), p(2,1)),'fontsize',18)
% figname = sprintf('wake_delta_power_vs_MSE_%i_%i',scales(1),scales(end));
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% %% Mean LZC vs Delta power, ASLEEP
% figure('Color','w')
% scatter(mean(sp_scl_lzc),log10(mean(delta_sp)),'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 1])
% mylsline
% ylabel('Power log_1_0(\muV^{2}/Hz)')
% xlabel(sprintf('LZC avg(%i-%i)',scales(1),scales(end)))
% [r,p] = corrcoef(log10(mean(delta_sp)),mean(sp_scl_lzc));
% %axis([0.15 0.4 2.25 4.25])
% title(sprintf('Sleep delta power vs LZC (scales %i-%i), r = %1.2f, p = %1.3f', scales(1),scales(end),r(2,1), p(2,1)),'fontsize',18)
% figname = sprintf('sleep_delta_power_vs_LZC_%i_%i',scales(1),scales(end));
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% 
% %% Mean LZC vs Delta power, AWAKE
% 
% figure('Color','w')
% scatter(mean(wk_scl_lzc),log10(mean(delta_wk)),'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 1])
% mylsline
% ylabel('Power log_1_0(\muV^{2}/Hz)')
% xlabel(sprintf('LZC avg(%i-%i)',scales(1),scales(end)))
% [r,p] = corrcoef(log10(mean(delta_wk)),mean(wk_scl_lzc));
% %axis([0.15 0.4 2.25 4.25])
% title(sprintf('Awake delta power vs LZC (scales %i-%i), r = %1.2f, p = %1.3f',scales(1),scales(end), r(2,1), p(2,1)),'fontsize',18)
% figname = sprintf('wake_delta_power_vs_LZC_%i_%i',scales(1),scales(end));
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% %% Mean MSE vs relative Delta power, ASLEEP
% figure('Color','w')
% scatter(mean(sp_scl),log10(mean(delta_sp)),'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 1])
% mylsline
% ylabel('Power log_1_0(\muV^{2}/Hz)')
% xlabel(sprintf('MSE avg(%i-%i)',scales(1),scales(end)))
% [r,p] = corrcoef(log10(mean(delta_sp)),mean(sp_scl));
% %axis([1.2 1.8 2.25 4.25])
% title(sprintf('Sleep delta relative power vs MSE (scales %i-%i), r = %1.2f, p = %1.3f', scales(1),scales(end),r(2,1), p(2,1)),'fontsize',18)
% figname = sprintf('sleep_delta_rel_power_vs_MSE_%i_%i',scales(1),scales(end));
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% 
% %% Mean MSE vs relative Delta power, AWAKE
% 
% figure('Color','w')
% scatter(mean(wk_scl),log10(mean(delta_wk)),'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 1])
% mylsline
% ylabel('Power log_1_0(\muV^{2}/Hz)')
% xlabel(sprintf('MSE avg(%i-%i)',scales(1),scales(end)))
% [r,p] = corrcoef(log10(mean(delta_wk)),mean(wk_scl));
% %axis([1.2 1.8 2.25 4.25])
% title(sprintf('Awake delta relative power vs MSE (scales %i-%i), r = %1.2f, p = %1.3f',scales(1),scales(end), r(2,1), p(2,1)),'fontsize',18)
% figname = sprintf('wake_delta_rel_power_vs_MSE_%i_%i',scales(1),scales(end));
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

%% Regress out effects of delta power, age, and genotype on MSE

B_wk = nan(4,size(wk_all,1),size(wk_all,2));
R_wk = nan(size(wk_all));
B_sp = nan(4,size(sp_all,1),size(sp_all,2));
R_sp = nan(size(sp_all));

for i = 1:size(wk_all,1)
    for j = 1:size(wk_all,2)
        [b,~,r] = regress(squeeze(wk_all(i,j,:)),[log10(delta_wk(i,:))' ...
            all_ages' allgenotype' ones(size(delta_wk,2),1)]);
        B_wk(:,i,j) = b;
        R_wk(i,j,:) = r;
    end
end

for i = 1:size(sp_all,1)
    for j = 1:size(sp_all,2)
        [b,~,r] = regress(squeeze(sp_all(i,j,:)),[log10(delta_sp(i,:))' ...
            all_ages' allgenotype' ones(size(delta_sp,2),1)]);
        B_sp(:,i,j) = b;
        R_sp(i,j,:) = r;
    end
end

% Use residuals to examine MSE with delta power effects regressed out

regout_wk = R_wk+repmat(squeeze(B_wk(4,:,:)),1,1,size(R_wk,3)); % regress out effect of delta power (residuals + y-intercept)
regout_sp = R_sp+repmat(squeeze(B_sp(4,:,:)),1,1,size(R_sp,3)); % regress out effect of delta power (residuals + y-intercept)

%% Also regress out Bayley score effects on MSE for subjects who have them

NDX = ~isnan(lng_score) & ~isnan(cog_score);

B_wk = nan(3,size(regout_wk,1),size(regout_wk,2));
R_wk = nan(size(regout_wk(:,:,NDX)));
B_sp = nan(3,size(regout_sp,1),size(regout_sp,2));
R_sp = nan(size(regout_sp(:,:,NDX)));

for i = 1:size(regout_wk,1)
    for j = 1:size(regout_wk,2)
        [b,~,r] = regress(squeeze(regout_wk(i,j,NDX)),[cog_score(NDX)' lng_score(NDX)' ones(length(cog_score(NDX)),1)]);
        B_wk(:,i,j) = b;
        R_wk(i,j,:) = r;
    end
end

for i = 1:size(regout_sp,1)
    for j = 1:size(regout_sp,2)
        [b,~,r] = regress(squeeze(regout_sp(i,j,NDX)),[cog_score(NDX)' lng_score(NDX)' ones(length(cog_score(NDX)),1)]);
        B_sp(:,i,j) = b;
        R_sp(i,j,:) = r;
    end
end

regout_wk(:,:,NDX) = R_wk+repmat(squeeze(B_wk(3,:,:)),1,1,size(R_wk,3)); % regress out effect of delta power (residuals + y-intercept)
regout_sp(:,:,NDX) = R_sp+repmat(squeeze(B_sp(3,:,:)),1,1,size(R_sp,3)); % regress out effect of delta power (residuals + y-intercept)

%% Mediation analysis (MSE)

sleep = cat(1,zeros(size(wk_scl,2),1),ones(size(wk_scl,2),1));
options.alpha = 0.05;
options.verbose = false;
options.display = false;
medP = nan(1,size(wk_scl,1));

for iscl = 1:size(wk_scl,1)
    out = mediationAnalysis0(sleep,[log10(mean(delta_wk)) ...
        log10(mean(delta_sp))]',[wk_scl(iscl,:) sp_scl(iscl,:)]',options);
    medP(iscl) = out.sobel.p;
end

medQ = mafdr(medP,'BHFDR',true);

for iscl = 1:size(wk_scl,1)
    if medQ(iscl) < 0.05
        fprintf('\n%s for MSE @ scale #%i, p = %1.2f FDR corrected\n','Full mediation',iscl,medQ(iscl))
    else
        fprintf('\n%s for MSE @ scale #%i, p = %1.2f FDR corrected\n','No mediation',iscl,medQ(iscl))
    end
end
        
figure('Color','w')
hold on
plot(log2(nyqfoi),-log10(medQ),'linewidth',2)
plot(log2(nyqfoi),ones(1,length(medQ)).*-log10(0.05),'k--')
ylim([0 3])
xlim([log2(5) log2(100)])
legend({'corrected p-value (Sobel test)','alpha level'},'location','best','fontsize',14)
legend boxoff
ticks = 2:1:7;
xticks(ticks)
ticklabels = cell(1,length(ticks));
for i = 1:length(ticklabels)
    ticklabels{i} = num2str(2^ticks(i));
end
xticklabels(ticklabels)
xlabel('Nyquist frequency (Hz)')
ylabel('-log_{10}(p)')
makefigpretty
ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation',...
    'right', 'Color','none');
makefigpretty
line(log2(nyqfoi),wsdelta,'Parent',ax2,'Color','none')
xlim([log2(5) log2(100)])
pltscl = fliplr([1 2 5 10 20]);
ticks = log2(nyqfoi(pltscl));
xticks(ticks)
ticklabels = cell(1,length(ticks));
for i = 1:length(ticklabels)
    ticklabels{i} = num2str(pltscl(i));
end
xticklabels(ticklabels)
xlabel('time scale')
set(gca,'ytick',[]);
yAX = get(gca,'YAxis');
yAX.Visible = 'off';
title(sprintf('MSE mediation analysis p-values'),'fontsize',18)
figname = sprintf('MSE mediation analysis p-values');
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,figname))

%% plot average MSE scales after regression

figure('Color','w')
hold on
[~,~,ci] = ttest(squeeze(mean(regout_wk))');
plot(log2(nyqfoi),mean(squeeze(mean(regout_wk)),2),'r','linewidth',4)
patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
    1.0,'EdgeColor','none')
[~,~,ci] = ttest(squeeze(mean(regout_sp))');
plot(log2(nyqfoi),mean(squeeze(mean(regout_sp)),2),'b','linewidth',4)
patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [ci(1,:) fliplr(ci(2,:))], 'b','facealpha', ...
    1.0,'EdgeColor','none')
legend('wake','wake 95% CI','sleep','95% sleep','AutoUpdate','off',...
    'FontSize',14,'location','northeast')
legend boxoff
ylim([0 2])
%title('MSE values over all time scales (no regression)','fontsize',18)
xlim([log2(5) log2(100)])
ticks = 2:1:7;
xticks(ticks)
ticklabels = cell(1,length(ticks));
for i = 1:length(ticklabels)
    ticklabels{i} = num2str(2^ticks(i));
end
xticklabels(ticklabels)
xlabel('Nyquist frequency (Hz)')
ylabel('sample entropy')
makefigpretty
% ADD TOP X-AXIS FOR TAU
ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation',...
    'right', 'Color','none');
makefigpretty
line(log2(nyqfoi),mean(squeeze(mean(regout_wk)),2),'Parent',ax2,'Color','none')
xlim([log2(5) log2(100)])
pltscl = fliplr([1 2 5 10 20]);
ticks = log2(nyqfoi(pltscl));
xticks(ticks)
ticklabels = cell(1,length(ticks));
for i = 1:length(ticklabels)
    ticklabels{i} = num2str(pltscl(i));
end
xticklabels(ticklabels)
xlabel('time scale')
set(gca,'ytick',[]);
yAX = get(gca,'YAxis');
yAX.Visible = 'off';
% SAVE FIGURE
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'MSE_scales_with_reg'))
print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,'MSE_scales_with_reg'))

%% MSE percent change from sleep after regressing out covariates

figure('Color','w')
hold on
wsdelta = (mean(squeeze(mean(regout_wk)),2) - mean(squeeze(mean(regout_sp)),2)) ... 
    ./mean(squeeze(mean(regout_sp)),2).*100;
[~,~,ci] = ttest(squeeze(mean(regout_wk)),squeeze(mean(regout_sp)),'dim',2); ci = ci';
plot(log2(nyqfoi),wsdelta ,'m','linewidth',4)
patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [ci(1,:)./mean(squeeze(mean(regout_sp)),2)'.*100 ... 
    fliplr(ci(2,:)./mean(squeeze(mean(regout_sp)),2)'.*100)], 'm','facealpha', 1.0,'EdgeColor','none')
legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',14,'location','southeast')
legend boxoff
plot(1:length(wsdelta),zeros(1,length(wsdelta)),'-.k')
xlim([log2(5) log2(100)])
ticks = 2:1:7;
xticks(ticks)
ticklabels = cell(1,length(ticks));
for i = 1:length(ticklabels)
    ticklabels{i} = num2str(2^ticks(i));
end
xticklabels(ticklabels)
%title('MSE difference over all time scales (no regression)','fontsize',18)
xlabel('Nyquist frequency (Hz)')
ylabel('SampEn change from sleep (%)')
makefigpretty
% ADD TOP X-AXIS FOR TAU
ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation',...
    'right', 'Color','none');
makefigpretty
line(log2(nyqfoi),wsdelta,'Parent',ax2,'Color','none')
xlim([log2(5) log2(100)])
pltscl = fliplr([1 2 5 10 20]);
ticks = log2(nyqfoi(pltscl));
xticks(ticks)
ticklabels = cell(1,length(ticks));
for i = 1:length(ticklabels)
    ticklabels{i} = num2str(pltscl(i));
end
xticklabels(ticklabels)
xlabel('time scale')
set(gca,'ytick',[]);
yAX = get(gca,'YAxis');
yAX.Visible = 'off';
% SAVE FIGURE
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'MSE_scales_delta_with_reg'))
print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,'MSE_scales_delta_with_reg'))


%% MSE permutation clustering after regression

clst_rgt = RMpermclusttest(regout_wk,regout_sp,Nperm);

sgcl_pos = [];
sgcl_neg = [];

if sum(clst_rgt.n_pos > 0) > 0 
    sgcl_pos = clst_rgt.P_val_pos < 0.05;
end

if sum(clst_rgt.n_neg > 0) > 0 
    sgcl_neg = clst_rgt.P_val_neg < 0.05;
end
cnt = 0;
for icls = 1:length(sgcl_pos)
    if sgcl_pos(icls)
        cnt = cnt +1;
        figure('Color','w')
        plot(log2(nyqfoi),sum(clst_rgt.lbl_pos==icls,2),'linewidth',2)
        ylim([0 20])
        xlim([log2(5) log2(100)])
        ticks = 2:1:7;
        xticks(ticks)
        ticklabels = cell(1,length(ticks));
        for i = 1:length(ticklabels)
            ticklabels{i} = num2str(2^ticks(i));
        end
        xticklabels(ticklabels)
        xlabel('Nyquist frequency (Hz)')
        ylabel('Number of channels')
        makefigpretty
        ax1 = gca; % current axes
        ax1_pos = ax1.Position; % position of first axes
        ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation',...
            'right', 'Color','none');
        makefigpretty
        line(log2(nyqfoi),wsdelta,'Parent',ax2,'Color','none')
        xlim([log2(5) log2(100)])
        pltscl = fliplr([1 2 5 10 20]);
        ticks = log2(nyqfoi(pltscl));
        xticks(ticks)
        ticklabels = cell(1,length(ticks));
        for i = 1:length(ticklabels)
            ticklabels{i} = num2str(pltscl(i));
        end
        xticklabels(ticklabels)
        xlabel('time scale')
        set(gca,'ytick',[]);
        yAX = get(gca,'YAxis');
        yAX.Visible = 'off';
        title(sprintf('MSE with reg postivie cluster, p = %1.3f',clst_rgt.P_val_pos(icls)),'fontsize',18)
        figname = sprintf('MSE with reg pos clust spect #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,figname))
        
        figure('Color','w')
        plot_topo_AS(sum(clst_rgt.lbl_pos==icls)')
        title(sprintf('MSE with reg pos clust, p = %1.3f',clst_rgt.P_val_pos(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(wk_lzc_all,2)])
        mycolorbar
        figname = sprintf('MSE with reg pos clust toop #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,figname))
    end
end
cnt = 0;
for icls = 1:length(sgcl_neg)
    if sgcl_neg(icls)
        cnt = cnt + 1;
        figure('Color','w')
        plot(log2(nyqfoi),sum(clst_rgt.lbl_neg==icls,2),'linewidth',2)
        ylim([0 20])
        xlim([log2(5) log2(100)])
        ticks = 2:1:7;
        xticks(ticks)
        ticklabels = cell(1,length(ticks));
        for i = 1:length(ticklabels)
            ticklabels{i} = num2str(2^ticks(i));
        end
        xticklabels(ticklabels)
        xlabel('Nyquist frequency (Hz)')
        ylabel('Number of channels')
        makefigpretty
        ax1 = gca; % current axes
        ax1_pos = ax1.Position; % position of first axes
        ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation',...
            'right', 'Color','none');
        makefigpretty
        line(log2(nyqfoi),wsdelta,'Parent',ax2,'Color','none')
        xlim([log2(5) log2(100)])
        pltscl = fliplr([1 2 5 10 20]);
        ticks = log2(nyqfoi(pltscl));
        xticks(ticks)
        ticklabels = cell(1,length(ticks));
        for i = 1:length(ticklabels)
            ticklabels{i} = num2str(pltscl(i));
        end
        xticklabels(ticklabels)
        xlabel('time scale')
        set(gca,'ytick',[]);
        yAX = get(gca,'YAxis');
        yAX.Visible = 'off';
        title(sprintf('MSE with reg  negative cluster, p = %1.3f',clst_rgt.P_val_neg(icls)),'fontsize',18)
        figname = sprintf('MSE neg clust spect #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,figname))
         
        figure('Color','w')
        plot_topo_AS(sum(clst_rgt.lbl_neg==icls)')
        title(sprintf('MSE neg clust, p = %1.3f',clst_rgt.P_val_neg(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(wk_lzc_all,2)])
        mycolorbar
        figname = sprintf('MSE with reg negs clust topo #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,figname))
    end
end

%% Regress out effects of delta power, age, and genotype on gMLZ

B_wk = nan(4,size(wk_lzc_all,1),size(wk_lzc_all,2));
R_wk = nan(size(wk_lzc_all));
B_sp = nan(4,size(sp_lzc_all,1),size(sp_lzc_all,2));
R_sp = nan(size(sp_lzc_all));

% del_wk_avg = log10(mean(delta_wk))';
% del_sp_avg = log10(mean(delta_sp))';

for i = 1:size(wk_lzc_all,1)
    for j = 1:size(wk_lzc_all,2)
        [b,~,r] = regress(squeeze(wk_lzc_all(i,j,:)),[log10(delta_wk(i,:))' ...
            all_ages' allgenotype' ones(size(delta_wk,2),1)]);
        B_wk(:,i,j) = b;
        R_wk(i,j,:) = r;
    end
end

for i = 1:size(sp_lzc_all,1)
    for j = 1:size(sp_lzc_all,2)
        [b,~,r] = regress(squeeze(sp_lzc_all(i,j,:)),[log10(delta_wk(i,:))' ...
            all_ages' allgenotype' ones(size(delta_wk,2),1)]);
        B_sp(:,i,j) = b;
        R_sp(i,j,:) = r;
    end
end

% Use residuals to examine gMLZ with delta power effects regressed out

reg_lzc_wk = R_wk+repmat(squeeze(B_wk(4,:,:)),1,1,size(R_wk,3)); % regress out effect of delta power (residuals + y-intercept)
reg_lzc_sp = R_sp+repmat(squeeze(B_sp(4,:,:)),1,1,size(R_sp,3)); % regress out effect of delta power (residuals + y-intercept)


%% Also regress out Bayley scores for subjects who have them

NDX = ~isnan(lng_score) & ~isnan(cog_score);

B_wk = nan(3,size(reg_lzc_wk,1),size(reg_lzc_wk,2));
R_wk = nan(size(reg_lzc_wk(:,:,NDX)));
B_sp = nan(3,size(reg_lzc_sp,1),size(reg_lzc_sp,2));
R_sp = nan(size(reg_lzc_sp(:,:,NDX)));

for i = 1:size(reg_lzc_wk,1)
    for j = 1:size(reg_lzc_wk,2)
        [b,~,r] = regress(squeeze(reg_lzc_wk(i,j,NDX)),[cog_score(NDX)' lng_score(NDX)' ones(length(cog_score(NDX)),1)]);
        B_wk(:,i,j) = b;
        R_wk(i,j,:) = r;
    end
end

for i = 1:size(reg_lzc_sp,1)
    for j = 1:size(reg_lzc_sp,2)
        [b,~,r] = regress(squeeze(reg_lzc_sp(i,j,NDX)),[cog_score(NDX)' lng_score(NDX)' ones(length(cog_score(NDX)),1)]);
        B_sp(:,i,j) = b;
        R_sp(i,j,:) = r;
    end
end

reg_lzc_wk(:,:,NDX) = R_wk+repmat(squeeze(B_wk(3,:,:)),1,1,size(R_wk,3)); % regress out effect of delta power (residuals + y-intercept)
reg_lzc_sp(:,:,NDX) = R_sp+repmat(squeeze(B_sp(3,:,:)),1,1,size(R_sp,3)); % regress out effect of delta power (residuals + y-intercept)

%% Mediation analysis (gMLZ)

sleep = cat(1,zeros(size(wk_scl_lzc,2),1),ones(size(wk_scl_lzc,2),1));
options.alpha = 0.05;
options.verbose = false;
options.display = false;
medP = nan(1,size(wk_scl_lzc,1));

for iscl = 1:size(wk_scl_lzc,1)
    out = mediationAnalysis0(sleep,[log10(mean(delta_wk)) ...
        log10(mean(delta_sp))]',[wk_scl_lzc(iscl,:) sp_scl_lzc(iscl,:)]',options);
    medP(iscl) = out.sobel.p;
end

medQ = mafdr(medP,'BHFDR',true);

for iscl = 1:size(wk_scl_lzc,1)
    if medQ(iscl) < 0.05
        fprintf('\n%s for gMLZ @ scale #%i, p = %1.2f FDR corrected\n','Full mediation',iscl,medQ(iscl))
    else
        fprintf('\n%s for gMLZ @ scale #%i, p = %1.2f FDR corrected\n','No mediation',iscl,medQ(iscl))
    end
end
    
%%
figure('Color','w')
h=subplot(1,1,1)
hold on
plot(log2(LZC_foi),-log10(medQ),'linewidth',2)
plot(log2(LZC_foi),ones(1,length(medQ)).*-log10(0.05),'k--')
makefigpretty
ylim([0 3])
xlim([0 log2(LZC_foi(end))])
legend({'corrected p-value (Sobel test)','alpha level'},'location','best','fontsize',14)
legend boxoff
xlabel('Center frequency (Hz)')
ylabel('-log10(p)')
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
title(sprintf('gMLZ mediation analysis p-values'),'fontsize',18)
figname = sprintf('gMLZ mediation analysis p-values');
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,figname))

%% plot average gMLZ scales after regression

figure('Color','w')
h=subplot(1,1,1);
hold on
[~,~,ci] = ttest(squeeze(mean(reg_lzc_wk))');
plot(log2(LZC_foi(scales)),mean(squeeze(mean(reg_lzc_wk)),2),'r','linewidth',4)
patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
    1.0,'EdgeColor','none')
[~,~,ci] = ttest(squeeze(mean(reg_lzc_sp))');
plot(log2(LZC_foi(scales)),mean(squeeze(mean(reg_lzc_sp)),2),'b','linewidth',4)
patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], 'b','facealpha', ...
    1.0,'EdgeColor','none')
legend('wake','wake 95% CI','sleep','95% sleep','AutoUpdate','off','FontSize',14,'location','southeast')
legend boxoff
xlim([1 length(scales)])
title('gMLZ values covary nuisance over all time scales ','fontsize',18)
xlabel('Center frequency (Hz)')
ylabel('LZ complexity')
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
xlim([0 5])
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_scales_no_reg'))

%% gMLZ percent change from sleep

figure('Color','w')
h=subplot(1,1,1);
hold on
wsdelta = (mean(squeeze(mean(reg_lzc_wk)),2) - mean(squeeze(mean(reg_lzc_sp)),2)) ... 
    ./mean(squeeze(mean(reg_lzc_sp)),2).*100;
[~,~,ci] = ttest(squeeze(mean(reg_lzc_wk)),squeeze(mean(reg_lzc_sp)),'dim',2); ci = ci';
plot(log2(LZC_foi(scales)),wsdelta ,'m','linewidth',4)
patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:)./mean(squeeze(mean(reg_lzc_sp)),2)'.*100 ...
    fliplr(ci(2,:)./mean(squeeze(mean(reg_lzc_sp)),2)'.*100)], 'r','facealpha', 1.0,'EdgeColor','none')
legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',14,'location','southeast')
legend boxoff
plot(log2(LZC_foi(scales)),zeros(1,length(wsdelta)),'-.k')
xlim([0 5])
%ylim([-0.04 0.04])
title('gMLZ covary nuisance difference over all time scales ','fontsize',18)
xlabel('Center frequency (Hz)')
ylabel('gLZ change from sleep (%)')
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_scales_delta_no_reg'))


%% gMLZ permutation clustering after regression

clst_LZC_rgt = RMpermclusttest(reg_lzc_wk,reg_lzc_sp,Nperm);

sgcl_pos = [];
sgcl_neg = [];

if sum(clst_LZC_rgt.n_pos > 0) > 0 
    sgcl_pos = clst_LZC_rgt.P_val_pos < 0.05;
end

if sum(clst_LZC_rgt.n_neg > 0) > 0 
    sgcl_neg = clst_LZC_rgt.P_val_neg < 0.05;
end
cnt = 0;

for icls = 1:length(sgcl_pos)
    if sgcl_pos(icls)
        cnt = cnt + 1;
        figure('Color','w')
        h=subplot(1,1,1);
        plot(log2(LZC_foi(scales)),sum(clst_LZC_rgt.lbl_pos==icls,2),'linewidth',2)
        ylim([0 20])
        xlim([0 5])
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        xlim([0 5])
        xlabel('Center frequency (Hz)')
        ylabel('Number of channels')
        makefigpretty
        title(sprintf('gMLZ postivie cluster covary nuisance, p = %1.3f',clst_LZC_rgt.P_val_pos(icls)),'fontsize',18)
        figname = sprintf('LZC pos clust spect #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,figname))
        
        figure('Color','w')
        plot_topo_AS(sum(clst_LZC_rgt.lbl_pos==icls)')
        title(sprintf('gMLZ pos clust, p = %1.3f',clst_LZC_rgt.P_val_pos(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(reg_lzc_wk,2)])
        mycolorbar
        figname = sprintf('LZC pos clust topo #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,figname))
    end
end

cnt = 0;
for icls = 1:length(sgcl_neg)
    if sgcl_neg(icls)
        cnt = cnt+1;
        figure('Color','w')
        h=subplot(1,1,1);
        plot(log2(LZC_foi(scales)),sum(clst_LZC_rgt.lbl_neg==icls,2),'linewidth',2)
        ylim([0 20])
        xlim([0 5])
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        xlim([0 5])
        xlabel('Center frequency (Hz)')
        ylabel('Number of channels')
        makefigpretty
        title(sprintf('gMLZ negative cluster covary nuisance, p = %1.3f',clst_LZC_rgt.P_val_neg(icls)),'fontsize',18)
        figname = sprintf('gMLZ neg clust spect #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,figname))
        
        figure('Color','w')
        plot_topo_AS(sum(clst_LZC_rgt.lbl_neg==icls)')
        title(sprintf('gMLZ neg clust, p = %1.3f',clst_LZC_rgt.P_val_neg(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(WP,2)])
        mycolorbar
        figname = sprintf('gMLZ neg clust topo #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-depsc','-r600','-painters',sprintf('%s%s.eps',DIRFIGURE,figname))
    end
end
    


% %% frequency vs -log10(p-val) corrected
% 
% figure('Color','w')
% h=subplot(1,1,1);
% foi_hd = 2.^[0:0.01:5]; % interpolation resolution for plotting
% lbl1 = '-log10(p) FDR corrected';
% lbl2 = 'p = 0.05';
% plot(log2(foi_hd),interp1(foi,-log10(Q)',foi_hd,'spline'),'r','LineWidth',2)
% hold on
% plot(log2(foi_hd),ones(1,length(foi_hd)).*-log10(0.05),'k--')
% legend({lbl1,lbl2},'AutoUpdate','off','FontSize',14)
% xlim([0 5])
% xlabel('Frequency (Hz)','FontSize',12)
% ylabel('-log_1_0(p-value)','FontSize',14)
% title('Delta power, sleep vs wake','FontSize',20)
% set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
% box off
% %axis square
% legend boxoff 
% ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
% annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
%     'LineStyle','none','FontSize',14);
% makefigpretty
% figname = sprintf('freq_-log10p');
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% 
% %% Deletion genotype PSD plots absolute power
% 
% figure('Color','w')
% h=subplot(1,1,1);
% [~,~,wkci] = ttest( log10(squeeze(nanmean(WP(:,:,logical(allgenotype)))))' );
% [~,~,spci] = ttest( log10(squeeze(nanmean(SP(:,:,logical(allgenotype)))))' );
% foi_hd = 2.^[0:0.01:5]; % interpolation resolution for plotting
% wk = log10(squeeze(nanmean(WP(:,:,logical(allgenotype)))))';
% sp = log10(squeeze(nanmean(SP(:,:,logical(allgenotype)))))';
% wkclr = [1 0 0];
% spclr = [0 0 1];
% lbl1 = 'Wake';
% lbl2 = 'Sleep';
% plot(log2(foi_hd),interp1(foi,nanmean(wk)',foi_hd,'spline'),'Color',wkclr,'LineWidth',4)
% hold on
% patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,wkci(1,:),foi_hd,'spline') ...
%     fliplr(interp1(foi,wkci(2,:),foi_hd,'spline'))],wkclr,'facealpha',1.0,'EdgeColor','none')
% plot(log2(foi_hd),interp1(foi,nanmean(sp)',foi_hd,'spline'),'Color',spclr,'LineWidth',4)
% patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,spci(1,:),foi_hd,'spline') ...
%     fliplr(interp1(foi,spci(2,:),foi_hd,'spline'))],spclr,'facealpha',1.0,'EdgeColor','none')
% legend({lbl1,sprintf('%s%s',lbl1,' 95% CI'),lbl2,sprintf('%s%s',lbl2,' 95% CI')},...
%     'AutoUpdate','off','FontSize',14)
% xlim([0 5])
% xlabel('Frequency (Hz)','FontSize',12)
% ylabel('Power log_1_0(\muV^{2}/Hz)','FontSize',14)
% title('Deletion Absolute Power','FontSize',20)
% set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
% box off
% %axis square
% legend boxoff 
% ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
% annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
%     'LineStyle','none','FontSize',14);
% set(gca,'linewidth',3)
% xAX = get(gca,'XAxis');
% set(xAX,'FontSize', 22)
% yAX = get(gca,'YAxis');
% set(yAX,'FontSize', 22)
% set(gca, 'TickDir', 'out')
% figname = sprintf('Deletion_PSDs_abs_power_sp_wk');
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% clear wkci spci wk sp
% 
% %% Non-deletion genotype PSD plots absolute power
% 
% figure('Color','w')
% h=subplot(1,1,1);
% [~,~,wkci] = ttest( log10(squeeze(nanmean(WP(:,:,logical(~allgenotype)))))' );
% [~,~,spci] = ttest( log10(squeeze(nanmean(SP(:,:,logical(~allgenotype)))))' );
% foi_hd = 2.^[0:0.01:5]; % interpolation resolution for plotting
% wk = log10(squeeze(nanmean(WP(:,:,logical(~allgenotype)))))';
% sp = log10(squeeze(nanmean(SP(:,:,logical(~allgenotype)))))';
% wkclr = [1 0 0];
% spclr = [0 0 1];
% lbl1 = 'Wake';
% lbl2 = 'Sleep';
% plot(log2(foi_hd),interp1(foi,nanmean(wk)',foi_hd,'spline'),'Color',wkclr,'LineWidth',4)
% hold on
% patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,wkci(1,:),foi_hd,'spline') ...
%     fliplr(interp1(foi,wkci(2,:),foi_hd,'spline'))],wkclr,'facealpha',1.0,'EdgeColor','none')
% plot(log2(foi_hd),interp1(foi,nanmean(sp)',foi_hd,'spline'),'Color',spclr,'LineWidth',4)
% patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,spci(1,:),foi_hd,'spline') ...
%     fliplr(interp1(foi,spci(2,:),foi_hd,'spline'))],spclr,'facealpha',1.0,'EdgeColor','none')
% legend({lbl1,sprintf('%s%s',lbl1,' 95% CI'),lbl2,sprintf('%s%s',lbl2,' 95% CI')},...
%     'AutoUpdate','off','FontSize',14)
% xlim([0 5])
% %ylim([-1 4])
% xlabel('Frequency (Hz)','FontSize',12)
% ylabel('Power log_1_0(\muV^{2}/Hz)','FontSize',14)
% title('Non-Deletion Absolute Power','FontSize',20)
% set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
% box off
% %axis square
% legend boxoff 
% ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
% annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
%     'LineStyle','none','FontSize',14);
% set(gca,'linewidth',3)
% xAX = get(gca,'XAxis');
% set(xAX,'FontSize', 22)
% yAX = get(gca,'YAxis');
% set(yAX,'FontSize', 22)
% set(gca, 'TickDir', 'out')
% figname = sprintf('Non-Deletion_PSDs_abs_power_sp_wk');
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% clear wkci spci wk sp
% 
% %% PSD plots relative power
% 
% figure('Color','w')
% h=subplot(1,1,1);
% [~,~,wkci] = ttest( log10(squeeze(nanmean(WPr)))' );
% [~,~,spci] = ttest( log10(squeeze(nanmean(SPr)))' );
% foi_hd = 2.^[0:0.01:5]; % interpolation resolution for plotting
% wk = log10(squeeze(nanmean(WPr)))';
% sp = log10(squeeze(nanmean(SPr)))';
% wkclr = [1 0 0];
% spclr = [0 0 1];
% lbl1 = 'Wake';
% lbl2 = 'Sleep';
% plot(log2(foi_hd),interp1(foi,nanmean(wk)',foi_hd,'spline'),'Color',wkclr,'LineWidth',4)
% hold on
% patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,wkci(1,:),foi_hd,'spline') ...
%     fliplr(interp1(foi,wkci(2,:),foi_hd,'spline'))],wkclr,'facealpha',1.0,'EdgeColor','none')
% plot(log2(foi_hd),interp1(foi,nanmean(sp)',foi_hd,'spline'),'Color',spclr,'LineWidth',4)
% patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,spci(1,:),foi_hd,'spline') ...
%     fliplr(interp1(foi,spci(2,:),foi_hd,'spline'))],spclr,'facealpha',1.0,'EdgeColor','none')
% legend({lbl1,sprintf('%s%s',lbl1,' 95% CI'),lbl2,sprintf('%s%s',lbl2,' 95% CI')},...
%     'AutoUpdate','off','FontSize',14)
% xlim([0 5])
% xlabel('Frequency (Hz)','FontSize',12)
% ylabel('Power log_1_0(1/log2(Hz))','FontSize',14)
% title('Relative Power','FontSize',20)
% set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
% box off
% %axis square
% legend boxoff 
% ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
% annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
%     'LineStyle','none','FontSize',14);
% set(gca,'linewidth',3)
% xAX = get(gca,'XAxis');
% set(xAX,'FontSize', 22)
% yAX = get(gca,'YAxis');
% set(yAX,'FontSize', 22)
% set(gca, 'TickDir', 'out')
% figname = sprintf('PSDs_rel_power_sp_wk');
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% %% PSD plots power variance
% 
% figure('Color','w')
% h=subplot(1,1,1);
% [~,~,wkci] = ttest( log10(squeeze(nanmean(WPv)))' );
% [~,~,spci] = ttest( log10(squeeze(nanmean(SPv)))' );
% foi_hd = 2.^[0:0.01:5]; % interpolation resolution for plotting
% wk = log10(squeeze(nanmean(WPv)))';
% sp = log10(squeeze(nanmean(SPv)))';
% wkclr = [1 0 0];
% spclr = [0 0 1];
% lbl1 = 'Wake';
% lbl2 = 'Sleep';
% plot(log2(foi_hd),interp1(foi,nanmean(wk)',foi_hd,'spline'),'Color',wkclr,'LineWidth',4)
% hold on
% patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,wkci(1,:),foi_hd,'spline') ...
%     fliplr(interp1(foi,wkci(2,:),foi_hd,'spline'))],wkclr,'facealpha',1.0,'EdgeColor','none')
% plot(log2(foi_hd),interp1(foi,nanmean(sp)',foi_hd,'spline'),'Color',spclr,'LineWidth',4)
% patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,spci(1,:),foi_hd,'spline') ...
%     fliplr(interp1(foi,spci(2,:),foi_hd,'spline'))],spclr,'facealpha',1.0,'EdgeColor','none')
% legend({lbl1,sprintf('%s%s',lbl1,' 95% CI'),lbl2,sprintf('%s%s',lbl2,' 95% CI')},...
%     'AutoUpdate','off','FontSize',14)
% xlim([0 5])
% xlabel('Frequency (Hz)','FontSize',12)
% ylabel('Power log_1_0(1/log2(Hz))','FontSize',14)
% title('Power variance','FontSize',20)
% set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
% box off
% %axis square
% legend boxoff 
% ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
% annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
%     'LineStyle','none','FontSize',14);
% set(gca,'linewidth',3)
% xAX = get(gca,'XAxis');
% set(xAX,'FontSize', 22)
% yAX = get(gca,'YAxis');
% set(yAX,'FontSize', 22)
% set(gca, 'TickDir', 'out')
% figname = sprintf('PSDs_power_var_sp_wk');
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))


%% Bayley cognition scores

% figure('Color','w')
% scatter(mean(wk_avg),cog_score,'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 1])
% [r,p] = corrcoef(mean(wk_avg),cog_score,'rows','complete');
% xlabel('awake MSE')
% ylabel('Bayley cognitive score')
% mylsline
% xlim([min(mean(wk_avg))-0.1 max(mean(wk_avg))+0.1])
% title(sprintf('r = %1.3f, p = %1.3f, scales = %i-%i',r(2,1),p(2,1),scales(1),scales(end)),'fontsize',18)
% makefigpretty
% figname = sprintf('cog_score_wk_mse_%i_%i',scales(1),scales(end));
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% figure('Color','w')
% scatter(mean(sp_avg),cog_score,'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 1])
% [r,p] = corrcoef(mean(sp_avg),cog_score,'rows','complete');
% xlabel('sleep MSE')
% ylabel('Bayley cognitive score')
% mylsline
% xlim([min(mean(sp_avg)) - 0.1 max(mean(sp_avg))+0.1])
% title(sprintf('r = %1.3f, p = %1.3f, scales = %i-%i',r(2,1),p(2,1),scales(1),scales(end)),'fontsize',18)
% makefigpretty
% figname = sprintf('cog_score_wk_mse_%i_%i',scales(1),scales(end));
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))


% %% plot R^2 as a function of time scale, Bayley cognitive scores vs MSE
% Rwk = nan(1,size(wk_scl,1));
% Rsp = nan(1,size(sp_scl,1));
% for iscale = 1:size(wk_scl,1)
%     rwk = corrcoef(wk_scl(iscale,:),cog_score,'rows','complete');
%     rsp = corrcoef(sp_scl(iscale,:),cog_score,'rows','complete');
%     Rwk(iscale) = rwk(2,1)^2;
%     Rsp(iscale) = rsp(2,1)^2;
% end
% 
% figure('Color','w')
% hold on
% plot(log2(nyqfoi),Rwk,'linewidth',3,'color','r')
% plot(log2(nyqfoi),Rsp,'linewidth',3,'color','b')
% legend({'Wake','Sleep'},'location','best','fontsize',14)
% legend boxoff
% xlim([log2(5) log2(100)])
% ticks = 2:1:7;
% xticks(ticks)
% ticklabels = cell(1,length(ticks));
% for i = 1:length(ticklabels)
%     ticklabels{i} = num2str(2^ticks(i));
% end
% xticklabels(ticklabels)
% xlabel('Nyquist frequency (Hz)')
% ylabel('R^{2}')
% makefigpretty
% % ADD TOP X-AXIS FOR TAU
% ax1 = gca; % current axes
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation',...
%     'right', 'Color','none');
% makefigpretty
% line(log2(nyqfoi),mean(wk_scl,2),'Parent',ax2,'Color','none')
% xlim([log2(5) log2(100)])
% pltscl = fliplr([1 2 5 10 20]);
% ticks = log2(nyqfoi(pltscl));
% xticks(ticks)
% ticklabels = cell(1,length(ticks));
% for i = 1:length(ticklabels)
%     ticklabels{i} = num2str(pltscl(i));
% end
% xticklabels(ticklabels)
% xlabel('time scale')
% set(gca,'ytick',[]);
% yAX = get(gca,'YAxis');
% yAX.Visible = 'off';
% % SAVE FIGURE
% xlabel('Time scale')
% ylabel('R^{2}')
% title('Variance of cognitive scores explained by MSE')
% makefigpretty
% figname = 'R2MSEcognition';
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% %% plot R^2 as a function of time scale, Bayley language scores vs MSE
% Rwk = nan(1,size(wk_scl,1));
% Rsp = nan(1,size(sp_scl,1));
% for iscale = 1:size(wk_scl,1)
%     rwk = corrcoef(wk_scl(iscale,:),lng_score,'rows','complete');
%     rsp = corrcoef(sp_scl(iscale,:),lng_score,'rows','complete');
%     Rwk(iscale) = rwk(2,1)^2;
%     Rsp(iscale) = rsp(2,1)^2;
% end
% 
% figure('Color','w')
% hold on
% plot(log2(nyqfoi),Rwk,'linewidth',3,'color','r')
% plot(log2(nyqfoi),Rsp,'linewidth',3,'color','b')
% legend({'Wake','Sleep'},'location','best','fontsize',14)
% legend boxoff
% xlim([log2(5) log2(100)])
% ticks = 2:1:7;
% xticks(ticks)
% ticklabels = cell(1,length(ticks));
% for i = 1:length(ticklabels)
%     ticklabels{i} = num2str(2^ticks(i));
% end
% xticklabels(ticklabels)
% xlabel('Nyquist frequency (Hz)')
% ylabel('R^{2}')
% makefigpretty
% % ADD TOP X-AXIS FOR TAU
% ax1 = gca; % current axes
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation',...
%     'right', 'Color','none');
% makefigpretty
% line(log2(nyqfoi),mean(wk_scl,2),'Parent',ax2,'Color','none')
% xlim([log2(5) log2(100)])
% pltscl = fliplr([1 2 5 10 20]);
% ticks = log2(nyqfoi(pltscl));
% xticks(ticks)
% ticklabels = cell(1,length(ticks));
% for i = 1:length(ticklabels)
%     ticklabels{i} = num2str(pltscl(i));
% end
% xticklabels(ticklabels)
% xlabel('time scale')
% set(gca,'ytick',[]);
% yAX = get(gca,'YAxis');
% yAX.Visible = 'off';
% % SAVE FIGURE
% xlabel('Time scale')
% ylabel('R^{2}')
% title('Variance of language scores explained by MSE')
% makefigpretty
% figname = 'R2MSElanguage';
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% 
% %% plot R^2 as a function of time scale, Bayley cognitive scores vs gMLZ
% Rwk = nan(1,size(wk_scl_lzc,1));
% Rsp = nan(1,size(sp_scl_lzc,1));
% for iscale = 1:size(wk_scl_lzc,1)
%     rwk = corrcoef(wk_scl_lzc(iscale,:),cog_score,'rows','complete');
%     rsp = corrcoef(sp_scl_lzc(iscale,:),cog_score,'rows','complete');
%     Rwk(iscale) = rwk(2,1)^2;
%     Rsp(iscale) = rsp(2,1)^2;
% end
% 
% figure('Color','w')
% hold on
% plot(Rwk,'linewidth',3,'color','r')
% plot(Rsp,'linewidth',3,'color','b')
% legend({'Wake','Sleep'},'location','southeast','fontsize',14)
% legend boxoff
% xlabel('Time scale')
% ylabel('R^{2}')
% legend('Wake','Sleep')
% legend boxoff
% title('Variance of cognitive scores explained by gMLZ')
% makefigpretty
% figname = 'R2LZCcognition';
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% %% plot R^2 as a function of time scale, Bayley language scores vs gMLZ
% Rwk = nan(1,size(wk_scl_lzc,1));
% Rsp = nan(1,size(sp_scl_lzc,1));
% for iscale = 1:size(wk_scl_lzc,1)
%     rwk = corrcoef(wk_scl_lzc(iscale,:),lng_score,'rows','complete');
%     rsp = corrcoef(sp_scl_lzc(iscale,:),lng_score,'rows','complete');
%     Rwk(iscale) = rwk(2,1)^2;
%     Rsp(iscale) = rsp(2,1)^2;
% end
% 
% figure('Color','w')
% hold on
% plot(Rwk,'linewidth',3,'color','r')
% plot(Rsp,'linewidth',3,'color','b')
% legend({'Wake','Sleep'},'location','southeast','fontsize',14)
% legend boxoff
% xlabel('Time scale')
% ylabel('R^{2}')
% legend('Wake','Sleep')
% legend boxoff
% title('Variance of language scores explained by gMLZ')
% makefigpretty
% figname = 'R2LZClanguage';
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 


%%

% Report stats from all time scales after averaging across channels

[h,pval3,ci,stat] = ttest(wk_scl',sp_scl'); %,'tail','right');
Q = mafdr(pval3,'BHFDR',true); % Bejamini Hochberg FDR
ngood = sum(Q < 0.05);


figure('Color','w')
plot_topo_AS(-log10(Q)')
title(sprintf('-log10(p), FDR corrected, scales %i -%i',scales(1),scales(end)),'fontsize',18)
colormap jet
caxis([0 8])
mycolorbar
figname = sprintf('sp_wk_log10p_%i_%i',scales(1),scales(end));
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
%%
figure('Color','w')
plot_topo_AS(tvals')
title(sprintf('t-stats, scales %i -%i',scales(1),scales(end)),'fontsize',18)
colormap jet
caxis([0 2])
mycolorbar
figname = sprintf('so_wk_tstats_%i_%i',scales(1),scales(end));
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

figure('Color','w')
plot_topo_AS(mean(wk_avg,2))
title(sprintf('Mean MSE values, wake condition, scales %i -%i',scales(1),scales(end)),'fontsize',18)
caxis([0 2.5])
colormap jet
mycolorbar
figname = sprintf('wk_mse_%i_%i',scales(1),scales(end));
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

figure('Color','w')
plot_topo_AS(mean(sp_avg,2))
caxis([0 2.5])
mycolorbar
title(sprintf('Mean MSE values, sleep condition, scales %i -%i',scales(1),scales(end)),'fontsize',18)
figname = sprintf('sp_mse_%i_%i',scales(1),scales(end));
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

figure('Color','w')
delta = mean(wk_avg,2) - mean(sp_avg,2);
plot_topo_AS(delta)
title(sprintf('Mean MSE difference, wake - sleep, scales %i -%i',scales(1),scales(end)),'fontsize',18)
colormap jet
caxis([-0.15 0.15])
mycolorbar
figname = sprintf('delta_mse_%i_%i',scales(1),scales(end));
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

% % analyze by age
%
% all_ages = nan(1,length(allage));
% for icell = 1:length(allage)
%     all_ages(icell) = mean(allage{icell});
% end
% all_ages = all_ages(~isnan(all_ages)); % purge NaNs

% sanity check
assert(length(wk_avg)==length(sp_avg),'Sleep and wake not balanced!')


fprintf('\n%i time scales are significant  after averaging across channels and correcting for multiple comparisons (FDR)\n',ngood)

%% Demographics figure

fprintf('n = %i\n',length(wk_avg))
fprintf('Age = %2.2f +/- %2.2f (yr)\n',mean(all_ages)/12,std(all_ages)/12)
fprintf('Using data from these subjects: \n')
% sanity check, make sure same subjects are used for sleep and wake

assert(length(SID_wk) == length(SID_sp),'Sleep and wake aren''t balanced!')
for isb = 1:length(SID_wk)
    assert(SID_wk{isb} == SID_sp{isb}) % sanity check
    fprintf('Paired MSE wake and sleep: %i %i\n',SID_wk{isb},SID_sp{isb})
end

fprintf('\n')

assert(length(SID_wk_lzc) == length(SID_sp_lzc),'Sleep and wake aren''t balanced!')
for isb = 1:length(SID_wk_lzc)
    assert(SID_wk_lzc{isb} == SID_sp_lzc{isb}) % sanity check
    fprintf('Paired LZC wake and sleep: %i %i\n',SID_wk_lzc{isb},SID_sp_lzc{isb})
end

fprintf('\n')

assert(length(SID_wk) == length(SID_wk),'Sleep and wake aren''t balanced!')
for isb = 1:length(SID_wk_lzc)
    assert(SID_wk_lzc{isb} == SID_wk{isb}) % sanity check
    fprintf('Paired wake data, MSE and LZC: %i %i\n',SID_wk{isb},SID_wk_lzc{isb})
end

figure('Color','w')
histogram(all_ages./12,1:11)
xlabel('age (years)')
ylabel('count (n)')
title('Angelman age distribution','fontsize',18)
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'AS_demographics'))

malecnt = 0;
delcnt = 0;

for i = 1:length(allgender)
    if allgender(i) == 1
        malecnt = malecnt + 1;
    end
    if allgenotype(i) == 1
        delcnt = delcnt + 1;
    end
end

fprintf('n = %i male \n',malecnt)
fprintf('n = %i deletion \n',delcnt)

figure('Color','w')
hold on
delta = mean(wk_avg) - mean(sp_avg);
scatter(all_ages,delta,'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 1])
[r,p] = corrcoef(all_ages,delta);

%% UNCOMMENT BELOW IF YOU WANT THE JACKKNIEFED CORRELATIONS (E.G., IF THERE
% ARE OUTLIERS)
% % jack-kniefed correlation coefficient (Spearman)
% jackrho = jackknife(@corr,all_ages',delta','type','spearman');
% meanrho = mean(jackrho);

plot(1:200,zeros(1,200),'-.r')
xlabel('age (months)')
ylabel('deltaMSE')
mylsline
title(sprintf('r = %1.3f, p = %1.3f, scales = %i-%i',r(2,1),p(2,1),scales(1),scales(end)),'fontsize',18)
makefigpretty
figname = sprintf('age_delta_%i_%i',scales(1),scales(end));
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% 
% switch match
%     case true
%         save MSE_scores_matched wk_scl sp_scl SID_wk SID_sp 
%     case false
%         save MSE_scores_all wk_scl sp_scl SID_wk SID_sp 
% end

% %% Correlate with LZC
% 
% figure('Color','w')
% scatter(mean(sp_scl),mean(sp_scl_lzc))
% mylsline
% ylabel(sprintf('LZ complexity avg(%i-%i)',scales(1),scales(end)))
% xlabel(sprintf('MSE avg(%i-%i)',scales(1),scales(end)))
% [r,p] = corrcoef(mean(sp_scl),mean(sp_scl_lzc));
% %axis([1.2 1.8 0.15 0.3])
% title(sprintf('Sleep MSE vs LZC (scales %i-%i), r = %1.2f, p = %1.2f', scales(1),scales(end),r(2,1), p(2,1)),'fontsize',18)
% figname = sprintf('sleep_MSE_vs_LZC_%i_%i',scales(1),scales(end));
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% figure('Color','w')
% scatter(mean(wk_scl),mean(wk_scl_lzc))
% mylsline
% ylabel(sprintf('LZ complexity avg(%i-%i)',scales(1),scales(end)))
% xlabel(sprintf('MSE avg(%i-%i)',scales(1),scales(end)))
% [r,p] = corrcoef(mean(wk_scl),mean(wk_scl_lzc));
% %axis([1.2 1.8 0.15 0.3])
% title(sprintf('Wake MSE vs LZC (scales %i-%i), r = %1.2f, p = %1.2f', scales(1),scales(end),r(2,1), p(2,1)),'fontsize',18)
% figname = sprintf('wake_MSE_vs_LZC_%i_%i',scales(1),scales(end));
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% 
% %% Do all MSE t-tests and correct for multiple testing
% 
% Pval_mse = nan(size(wk_all,1),(size(wk_all,2)));
% Tval_mse = nan(size(wk_all,1),(size(wk_all,2)));
% 
% for ichan = 1:size(wk_all,1)
%     for jscale = 1:size(wk_all,2)
%         [~,p,~,stat] = ttest(wk_all(ichan,jscale,:),sp_all(ichan,jscale,:));
%         Pval_mse(ichan,jscale) = p;
%         Tval_mse(ichan,jscale) = stat.tstat;
%     end
% end
% 
% alpha = 0.05/(size(wk_all,1)*size(wk_all,2));
% Bonf_mse = Pval_mse < alpha;
% Q_mse = mafdr(Pval_mse(1:end),'BHFDR',true); % Bejamini Hochberg FDR
% Q_mse = reshape(Q_mse,size(wk_all,1),size(wk_all,2));
% 
% %% Do all LZC t-stats and correct for multiple testing 
% 
% Pval_lzc = nan(size(wk_lzc_all,1),(size(wk_lzc_all,2)));
% Tval_lzc = nan(size(wk_lzc_all,1),(size(wk_lzc_all,2)));
% 
% for ichan = 1:size(wk_lzc_all,1)
%     for jscale = 1:size(wk_lzc_all,2)
%         [~,p,~,stat] = ttest(wk_lzc_all(ichan,jscale,:),sp_lzc_all(ichan,jscale,:));
%         Pval_lzc(ichan,jscale) = p;
%         Tval_lzc(ichan,jscale) = stat.tstat;
%     end
% end
% 
% 
% alpha = 0.05/(size(wk_lzc_all,1)*size(wk_lzc_all,2));
% Bonf_lzc = Pval_lzc < alpha;
% Q_lzc = mafdr(Pval_lzc(1:end),'BHFDR',true); % Bejamini Hochberg FDR
% Q_lzc = reshape(Q_lzc,size(wk_lzc_all,1),size(wk_lzc_all,2));
% 
% %% FDR correction accounting for all comparisons across MSE and gMLZ
% 
% Pval_all = [Pval_mse(1:end) Pval_lzc(1:end)];
% Tval_all = [Tval_mse(1:end) Tval_lzc(1:end)];
% Q_all   = mafdr(Pval_all,'BHFDR','true');
% Q_wake  = Q_all;
% Q_sleep = Q_all;
% Q_wake(Tval_all < 0) = NaN; % just those tests that found greater complexity in the awake state
% Q_sleep(Tval_all > 0) = NaN; % just those tests that found greater complexity in the asleep state
% 
% fprintf('\n %i MSE and %i gMLZ comparisons are statistically significant\n', ...
%     sum(Q_all(1:size(wk_all,1)*size(wk_all,2)) < 0.05),sum(Q_all(size(wk_all,1)*size(wk_all,2)+1:end) < 0.05));
% 
% fprintf('\n %i MSE and %i gMLZ comparisons are greater in the awake state\n', ...
%     sum(Q_wake(1:size(wk_all,1)*size(wk_all,2)) < 0.05),sum(Q_wake(size(wk_all,1)*size(wk_all,2)+1:end) < 0.05));
% 
% fprintf('\n %i MSE and %i gMLZ comparisons are greater in the asleep state\n', ...
%     sum(Q_sleep(1:size(wk_all,1)*size(wk_all,2)) < 0.05),sum(Q_sleep(size(wk_all,1)*size(wk_all,2)+1:end) < 0.05));
% 
% % build table 
% 
% load joel_lay.mat lay
% channels = lay.label';
% scalevec = ceil(linspace(10e-5,20,380))';
% 
% stat_table = table([repmat({'MSE'},380,1); repmat({'gMLZ'},380,1)],...
%     repmat(channels,40,1),repmat(scalevec,2,1),(Q_all < 0.05)',Pval_all',Q_all',Tval_all');
% 
% writetable(stat_table,'Stats_table.csv')


%% Export subject table

if match
    sex = cell(1,1);
    for i = 1:length(allgender)
        if allgender(i)
            sex{i} = 'male';
        else
            sex{i} = 'female';
        end
    end
    % sanity check 
    subtab = table(all_ages',sex',logical(allgenotype)',round(dur_wk,1)',round(dur_sp,1)',round(trg_wk,1)',...
        round(trg_wk_seg,1)',round(trg_sp,1)',round(trg_sp_seg,1)');
    writetable(subtab,'Table1.csv')
end

%% Save variables for analysis in R 
% 
% MSE = cat(4,wk_all,sp_all);
% LZC = cat(4,wk_lzc_all,sp_lzc_all);
% 
% assert(length(all_ages)==size(MSE,3),'Wrong number of subject ages')
% assert(length(allgenotype)==size(MSE,3),'Wrong number of subject genotypes')
% assert(length(allgender)==size(MSE,3),'Wrong number of subject genders')
% 
% age = [];
% genotype = [];
% sex = [];
% 
% for isb = 1:size(MSE,3)
%     age = [age; repmat(all_ages(isb),size(MSE,1)*size(MSE,2)*size(MSE,4),1)];
%     genotype = [genotype; repmat(allgenotype(isb),size(MSE,1)*size(MSE,2)*size(MSE,4),1)];
%     sex = [sex; repmat(allgender(isb),size(MSE,1)*size(MSE,2)*size(MSE,4),1)];
% end
% 
% % mean across channels and scales
% mse_avg = [mean(wk_scl); mean(sp_scl)];
% lzc_avg = [mean(wk_scl_lzc); mean(sp_scl_lzc)];
% 
% save('MSE','MSE')
% save('LZC','LZC')
% save('MSEAVG','mse_avg')
% save('LZCAVG','lzc_avg')
% save('AS_AGE','age')
% save('AS_GENOTYPE','genotype')
% save('AS_SEX','sex')
% %%
% channel = ones(size(MSE));
% for irow = 1:size(channel,1)
%     channel(irow,:,:,:) = channel(irow,:,:,:).*irow;
% end
% 
% scale = ones(size(MSE));
% for icol = 1:size(scale,2)
%     scale(:,icol,:,:) = scale(:,icol,:,:).*icol;
% end
% 
% state = ones(size(MSE));
% for ily = 1:size(state,4)
%     state(:,:,:,ily) = state(:,:,:,ily).*ily;
% end
% 
% group = {channel(1:end) scale(1:end) state(1:end)};
% 
% p = anovan(MSE(1:end),group,'varnames',char('Channel', 'Scale','Sleep'));
% 

%statmat = [wk_avg sp_avg];


%
% % two-way anova
% %[p,tbl,stats] = anova2(statmat,size(wk_avg,2));
%
% chankern = 1:size(statmat,1); chankern = chankern';
% chanmat = repmat(chankern,1,size(statmat,2),size(statmat,3));
% spconmat = [zeros(size(statmat,1),size(statmat,2)/2) ...
%     ones(size(statmat,1),size(statmat,2)/2)];
% spconmat = repmat(spconmat,1,1,size(statmat,3));
% subjects = 1:size(wk_avg,2);
% submat = repmat(subjects,19,2);
%
% statvec = statmat(1:end);
% chanvec = chanmat(1:end);
% spconvec = spconmat(1:end);
% subvec = submat(1:end);
%
% % repeated measures two-way anova
% stats = rm_anova2(statvec,subvec,chanvec,spconvec,{'channel','sleep'})
%
% %c = multcompare(stats)
%
% % one-way anova on sleep (channel factor)
% sleep_table = table();
% for ich = 1:size(sp_avg,1)
%     sleep_table = setfield(sleep_table,strcat('ch',num2str(ich)),sp_avg(ich,:));
% end
%
% T = table(spconmat(1,:)',statmat(1,:)',statmat(1,:)',statmat(1,:)',statmat(1,:)',...
%     statmat(1,:)',statmat(1,:)',statmat(1,:)',statmat(1,:)',statmat(1,:)',...
%     statmat(1,:)',statmat(1,:)',statmat(1,:)',statmat(1,:)',statmat(1,:)',...
%     statmat(1,:)',statmat(1,:)',statmat(1,:)',statmat(1,:)',statmat(1,:)',...
%     'VariableNames',{'Sleep','Ch1','Ch2','Ch3','Ch4','Ch5','Ch6','Ch7',...
%     'Ch8','Ch9','Ch10','Ch11','Ch12','Ch13','Ch14','Ch15','Ch16','Ch17',...
%     'Ch18','Ch19',});
%
% vars = table( [ones(1,19); ones(1,19).*2]', [1:19]','VariableNames',{'SleepCondition','Channels'})
% %chans = 1:19';
% rm = fitrm(T,'Ch1-Ch19 ~ Sleep','WithinDesign',vars)
%
%
% %
% % [P,ANOVATAB,STATS] = anova1(sp_avg');
% %
% %
% %
% % % one-way anova on sleep
% %
% % [P,ANOVATAB,STATS] = anova1(wk_avg');
