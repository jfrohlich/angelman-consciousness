%function[wk_avg,sp_avg,pvals] = ana_AS_stats_cluster_short(scales,p_thresh,match)

% check if we are in the right directory
if ~contains(pwd,'Monti')
    warning('It looks like you might be in the wrong directory!')
end

% input: scales are the time scales we want to average
clearvars
rng(345800);  % For reproducibility (jack-kniefed correlations)

%if ~exist('nargin','var') || nargin < 1 % If we are executing as a script OR if there are no inpu args
%     scales = 1:20;
%     p_thresh = 0.01; %0.0005; % p-value treshold for cluster permutation stats
%     match = true; 
%end

scales = 1:20;
p_thresh = 0.01; %0.0005; % p-value treshold for cluster permutation stats
match = true; 

fsample = 200; % all files should be downsampled to this freq before computing complexity measures
nyqfoi = (fsample/2)./scales; % MSE
strict = true; % do we filter out clinical EEG?
questionable = {'AS_100007_20180706','AS_104972_20170515'};

close all
dbstop if error
viewts = false;
Nperm = 10^4; % number of permutations
surrtest = false; % surrogate data testing?
method = 'gMLZ'; % Lempel-Ziv method
myload 1020_labels
if ~match
    myload select_files
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

% add subdirectory for p-value threshold

DIRFIGURE = sprintf('%sp_=_%0.4f/',DIRFIGURE,p_thresh);

if ~exist(DIRFIGURE,'dir'), mkdir(DIRFIGURE), end

switch match
    case true
        pth = './MSE_March/match/Xie/r=0.15/dynr/';
        wake = mydir(sprintf('%swake/*.mat',pth));
        sleep = mydir(sprintf('%ssleep/*.mat',pth));
        %frqpth = './freq_out/March/nan=0.2/match/';
        frqpth = './freq_out/July/match/';
        frqwake  = mydir(sprintf('%swake/*.mat',frqpth));
        frqsleep = mydir(sprintf('%ssleep/*.mat',frqpth));
        pthLZC = './gLZC_March/match/';
    case false
        pth = './MSE_March/Xie/r=0.15/dynr/';
        wake = mydir(sprintf('%swake/*.mat',pth));
        sleep = mydir(sprintf('%ssleep/*.mat',pth));
        %frqpth = './freq_out/March/nan=0.2/';
        frqpth = './freq_out/July/';
        frqwake  = mydir(sprintf('%swake/*.mat',frqpth));
        frqsleep = mydir(sprintf('%ssleep/*.mat',frqpth));
        pthLZC = './gLZC_March/';
end
pthTD = './MSE_TD/dynr/';
pthTD_LZC = './gLZC_TD/';
wakeTD= mydir(sprintf('%s/*.mat',pthTD));
wakeTD_LZC= mydir(sprintf('%s/*.mat',pthTD_LZC));

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
alldate = cell(1,1);
allgender = [];
allgenotype = [];
allsite = [];
allcog = cell(1,1);
alllng = cell(1,1);
allVCog = cell(1,1);
allNVCog = cell(1,1);
allSTN  = cell(1,1);
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

T  = readtable('AS_table_2019_JF.csv','Delimiter',','); % exported from the master spreadsheet
B  = readtable('AS_Bayley.xlsx');
SS = readtable('SID_sites.xlsx');

OS = computer; % detect operating system
switch OS
    case 'PCWIN64'
        WindowsProof % script to fix file strings
end
% 
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
%     % MSE
%     
%     % get SID and age
%     underscore = strfind(wakeTD(ifile).name,'_');
%     stop = underscore(3)-1; % everything up to the third underscore 
%     fstr = wakeTD(ifile).name(1:stop);
%     SID = fstr(underscore(1)+1:underscore(2)-1);
%     age = fstr(underscore(2)+1:strfind(fstr,'m')-1);
%     
%     % load empirical values
%     W_TD = myload(sprintf('%s/%s',pthTD,wakeTD(ifile).name));
%     fprintf('\nSuccessfully loaded MSE from TD %s\n',SID)
%     
%     % find "buried" data structures (these got buried due to a glitch in
%     % the script use gather() to turn gpuArray back to regular array
%    
%     % recursively find the real structure
%     while isfield(W_TD.MSEout,'MSEout')
%         W_TD.MSEout = W_TD_Slzc.MSEout.MSEout;
%     end
%     
%     min_valid_mse = 100; % Grandy et al., 2016
%     wk_gdwin = W_TD.MSEout.n_valid >= min_valid_mse;
%     if sum(wk_gdwin) == 0 % at least one good window!
%         fprintf('Skipping this file, no good windows: %s\n', wakeTD(ifile).name)
%         continue
%     end
%     
%     % LZC
%     
%     wk_lzc_scale = [];
%     wk_lzc_avg = [];
%     
%     % load empirical values
%     W_TD_lzc = myload(sprintf('%s/%s',pthTD_LZC,wakeTD_LZC(ifile).name));
%     fprintf('\nSuccessfully loaded LZC from TD %s\n',SID)
%     
%     % find "buried" data structures (these got buried due to a glitch in
%     % the script use gather() to turn gpuArray back to regular array
%    
%     % recursively find the real structure
%     while isfield(W_TD_lzc.LZCout,'LZCout')
%         W_TD_lzc.LZCout = W_TD_lzc.LZCout.LZCout;
%     end
%     
%     min_valid_lzc = 2000; % Grandy et al., 2016
%     
%     % find good windows for multiscale (gMLZ or dLZC)
%     switch method
%         case 'gMLZ'
%             wk_lzc_gdwin = W_TD_lzc.LZCout.n_valid_gLZC >= min_valid_lzc;
%             %sp_lzc_gdwin = Slzc.n_valid_gLZC >= min_valid_lzc;
%             wk_lzc_win_avg = nanmean(W_TD_lzc.LZCout.gLZC(wk_lzc_gdwin),3);
%             %sp_lzc_win_avg = nanmean(Slzc.gLZC(sp_lzc_gdwin),3);
%         case 'dLZC'
%             wk_lzc_gdwin = W_TD_lzc.LZCout.n_valid_dLZC(:,1) >= min_valid_lzc;
%             %sp_lzc_gdwin = Slzc.n_valid_dLZC(:,1) >= min_valid_lzc;
%             wk_lzc_win_avg = nanmean(W_TD_lzc.LZCout.dLZC(wk_lzc_gdwin),3);
%             %sp_lzc_win_avg = nanmean(Slzc.dLZC(sp_lzc_gdwin),3);
%         otherwise
%             error('Method unknown')
%     end
%     
%     if sum(wk_lzc_gdwin(:)) == 0 % at least one good window!
%         fprintf('Skipping this file, no good windows: %i\n', SlzcID)
%         continue
%     end
%     
%     % load into variables 
%     
%     TD_lzc_scl(:,ifile) = nanmean(nanmean(W_TD_lzc.LZCout.gLZC,3));
%     TD_lzc_avg(:,ifile)   = nanmean(nanmean(W_TD_lzc.LZCout.gLZC,3),2)';
%     
%     wk_win_avg = nanmean(W_TD.MSEout.mse(:,:,wk_gdwin),3);
%     
%     mse_td{ifile} = wk_win_avg;
%     win_td(ifile) = sum(W_TD.MSEout.n_valid(wk_gdwin));
%     dur_td_all(ifile) = W_TD.MSEout.dur_used;
%     SID_td(ifile) = str2num(SID);
%     age_td(ifile) = str2num(age);
%     
%     td_avg(:,ifile) = mean(mse_td{ifile}(:,scales),2);
%     td_scl(:,ifile) = mean(mse_td{ifile}(:,scales),1);
%     td_all      = cat(3,td_all,mse_td{ifile});
% 
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
    
    % Load sleep data
    St = readtable('SleepData.xlsx');
    % Number of nights per week that child sleeps through the night
    STN = St.sleeps_through_night_num(St.ReferenceID == SID); 
    
    RSP = B.RegularSleepingPattern(B.ID == SID & round(B.Age) == age);
    if isempty(RSP) % if empty, take the closest age
        [~,idx] = min(abs(B.Age(B.ID==SID) - age));
        RSPall = B.RegularSleepingPattern(B.ID == SID);
        RSP = RSPall(idx);
    end
    
    if strcmpi(RSP,'Yes')
        RSP = 1; % use dummy variables (T/F doesn't work with NaNs) 
    elseif strcmpi(RSP,'No')
        RSP = 0; % use dummy variables (T/F doesn't work with NaNs) 
    else 
        RSP = NaN;
    end
    
    site_str = unique(SS.Site(SS.SID == SID));
    if strcmp(site_str,'SD')
        site = 1; % dummy variable
    elseif strcmp(site_str,'Bos')
        site = 0; % dummy variable
    else
        error('Site not recognized')
    end
    assert(length(site)==1,'Not one unique site found')
    
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
    
    %% load power (do this first so we can skip EEGs without valid power output)
    
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
            fprintf('     loading power...\n')
            myload(sprintf('%s%s',frqpth,'sleep/',frqsleep(jfile).name),'pow','pow_full','pow_var','foi','n','cfg')
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
        %fprintf('\n Bad sleep power, only %i good windows, skipping this one ...\n',n(1))
        %keyboard
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
            fprintf('     loading power...\n')
            myload(sprintf('%s%s',frqpth,'wake/',frqwake(jfile).name),'pow','pow_full','pow_var','foi','n','cfg')
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
        %keyboard
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
    
    %% load MSE
    
    % load empirical values
    We = myload(sprintf('%swake/%s%s',pth,fstr,'_wake_MSE'));
    
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
            Ws = myload(sprintf('%swake/Surrogate/%s',pth,sprintf('%s%s',...
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
    
    fprintf('\nloading awake file %s \n',wake(ifile).name)
    try % attempt to load the corresponding sleep file
        sleepfile = strcat(fstr,'_sleep_MSE');
        
        % load empirical values
        Se = myload(sprintf('%ssleep/%s.mat',pth,sleepfile));
        
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
            Ss = myload(sprintf('%ssleep/Surrogate/%s.mat',pth,sprintf('%s%s',...
                sleepfile,'_SURROGATE')));
        catch
            %             Ss = myload(sprintf('%ssleep/Surrogate/%s.mat',pth,sprintf('%s%s',...
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
    
    %% load LZC
    
    % load empirical values
    %try
    Welzc = myload(sprintf('%swake/%s%s',pthLZC,fstr,'_wake_LZC'));
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
        Wslzc = myload(sprintf('%swake/Surrogate/%s',pthLZC,sprintf('%s%s',...
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
    fprintf('\nloading awake file %s \n',fstr)
    
    try % attempt to load the corresponding sleep file
        sleepfile = strcat(fstr,'_sleep_LZC');
        
        % load empirical values
        Selzc = myload(strcat(sprintf('%ssleep/',pthLZC),sleepfile));
        
        fprintf('\n loading sleep file %s \n',sleepfile)
    catch
        fprintf('\n No sleep data for this subject ....\n')
        continue
    end
    
    if isfield(Selzc.LZCout,'surrogate') % not all files will have this field
        assert(~Selzc.LZCout.surrogate,'This empirical file is actually a surrogate')
    end
    assert(Selzc.LZCout.cfg.new_srate == fsample,'File wasn''t downsampled?')
    
    if surrtest
        % load surrogate values
        Sslzc = myload(sprintf('%ssleep/Surrogate/%s%s.mat',pthLZC,...
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
        alldate{idx} = cat(2,alldate{idx},date);
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
        alldate{idx} = date;
        % Gender, genotype, and site -- you only need to add these once,
        % they don't change
        allgender(idx) = gender;
        allgenotype(idx) = genotype;
        allsite(idx) = site;
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
    if idx > length(allSTN) || isempty(allSTN(idx)) % if this is the first visit for this participant
        allSTN{idx} = STN; % fill in the entry for this visit
        if ~isempty(STN) 
            allSTN{idx} = STN;
        else
            allSTN{idx} = NaN;
        end
    else % if we already have at least one visit for this participant
        idx2 = length(allSTN{idx})+1; % fill in the next visit
        if ~isempty(STN)
            allSTN{idx}(idx2) = STN;
        else
            allSTN{idx}(idx2) = NaN;
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
allsite(badidx) = [];
allgender(badidx) = [];
allage(badidx) = [];
alldate(badidx) = [];
allNVCog(badidx) = [];
allVCog(badidx) = [];
allcog(badidx) = [];
alllng(badidx) = [];
allSTN(badidx) = [];
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
            all_dates(i) = alldate{i}(gdx);
            cog_score(i) = allcog{i}(gdx);
            lng_score(i) = alllng{i}(gdx);
            VC_score(i) = allVCog{i}(gdx);
            NVC_score(i) = allNVCog{i}(gdx);
            STN_TF(i)    = allSTN{i}(gdx);
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
            all_dates(i) = alldate{i};
            cog_score(i) = allcog{i};
            lng_score(i) = alllng{i};
            STN_TF(i)    = allSTN{i};
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

% % Find the appropriate age match
% keep = age_td > mean(all_ages)-std(all_ages)*2 & age_td < ...
%     mean(all_ages)+std(all_ages)*2 & age_td > 12 & age_td < 12*18;
% 
% td_scl = td_scl(:,keep);
% TD_lzc_scl = TD_lzc_scl(:,keep);
% TD_lzc_avg = TD_lzc_avg(:,keep);
% td_avg = td_avg(:,keep);


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

%myfigure
myfigure
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
'FontSize',18,'location','northeast')
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
ylabel('mSampEn'); ylim([0 2]);
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
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'MSE_scales_no_reg'))

% %% plot average awake MSE scales w/TD
% 
% %myfigure
% myfigure
% hold on
% [~,~,ci] = ttest(wk_scl');
% plot(log2(nyqfoi),mean(wk_scl,2),'r','linewidth',2)
% patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
%     0.3,'EdgeColor','none')
% 
% [~,~,ci] = ttest(td_scl');
% plot(log2(nyqfoi),mean(td_scl,2),'color',[0 0.5 0.5],'linewidth',2)
% patch( [log2(nyqfoi) fliplr(log2(nyqfoi))], [ci(1,:) fliplr(ci(2,:))],[0 0.5 0.5],'facealpha', ...
%     0.3,'EdgeColor','none')
% 
% % [~,~,ci] = ttest(sp_scl');
% % plot(log2(nyqfoi),mean(sp_scl,2),'b','linewidth',4)
% % patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [ci(1,:) fliplr(ci(2,:))], 'b','facealpha', ...
% %     1.0,'EdgeColor','none')
% legend('AS wake','AS wake 95% CI','TD wake','TD wake 95% CI','AutoUpdate','off',...
%     'FontSize',18,'location','northeast')
% legend boxoff
% %title('MSE values over all time scales (no regression)','fontsize',18)
% xlim([log2(5) log2(100)])
% ticks = 2:1:7;
% xticks(ticks)
% ticklabels = cell(1,length(ticks));
% for i = 1:length(ticklabels)
%     ticklabels{i} = num2str(2^ticks(i));
% end
% xticklabels(ticklabels)
% xlabel('Nyquist frequency (Hz)')
% ylabel('mSampEn'); ylim([0 2]);
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
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'MSE_and_TD_scales_no_reg'))
% print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'MSE_and_TD_scales_no_reg'))



%% plot MSE percent change

myfigure
hold on
wsdelta = ((mean(wk_scl,2) - mean(sp_scl,2))./mean(sp_scl,2)).*100;
[~,~,ci] = ttest(wk_scl,sp_scl,'dim',2); ci = ci';
plot(log2(nyqfoi),wsdelta ,'m','linewidth',4)
patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [(ci(1,:)./mean(sp_scl,2)').*100 ...
    fliplr((ci(2,:)./mean(sp_scl,2)').*100)], 'm','facealpha', 1.0,'EdgeColor','none')
legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',18,'location','southeast')
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
ylabel('mSampEn change from sleep (%)'); ylim([0 40]);
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
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'MSE_scales_delta_no_reg'))

%% plot MSE raw change

myfigure
hold on
wsdelta = mean(wk_scl,2) - mean(sp_scl,2);
[~,~,ci] = ttest(wk_scl,sp_scl,'dim',2); ci = ci';
plot(log2(nyqfoi),wsdelta ,'m','linewidth',4)
patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [ci(1,:) fliplr(ci(2,:))], 'm','facealpha', 1.0,'EdgeColor','none')
legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',18,'location','southeast')
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
ylabel('mSampEn change from sleep'); 
ylim([-0.2 0.2]);
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
ylim([-0.1 0.1])
% SAVE FIGURE
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'MSE_scales_delta_no_reg_RAW'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'MSE_scales_delta_no_reg_RAW'))


%% MSE permutation clustering

clst = RMpermclusttest(wk_all,sp_all,Nperm,p_thresh);

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
        myfigure
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
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
        
        myfigure
        plot_topo_AS(sum(clst.lbl_pos==icls)')
        title(sprintf('MSE pos clust, p = %1.3f',clst.P_val_pos(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(wk_lzc_all,2)])
        mycolorbar
        figname = sprintf('MSE pos clust toop #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
    end
    
    % dot plots for all clusters (even if p > 0.05)
    if any(~isnan(clst.lbl_pos(:)))   
        dot_slope_plot(wk_all,sp_all,clst.lbl_pos,icls,all_ages,allgenotype,p_thresh)
        title(sprintf('MSE, d_median = %1.2f',clst.d_pos(icls)),'fontsize',30)
        ylabel('SampEn')
        %ylim([0.7 1.7])
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s%i.svg',DIRFIGURE,'MSE_pos_clust_dot_plot_#',icls))
        print(gcf,'-dpng','-r300','-painters',sprintf('%s%s%i.png',DIRFIGURE,'MSE_pos_clust_dot_plot_#',icls))
    end
end

cnt = 0;
for icls = 1:length(sgcl_neg)
    if sgcl_neg(icls)
        cnt = cnt + 1;
        myfigure
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
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
        
        myfigure
        plot_topo_AS(sum(clst.lbl_neg==icls)')
        title(sprintf('MSE neg clust, p = %1.3f',clst.P_val_neg(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(wk_lzc_all,2)])
        mycolorbar
        figname = sprintf('MSE negs clust topo #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
    end
    
    % dot plots for all clusters (even if p > 0.05)
    if any(~isnan(clst.lbl_neg(:)))    
        dot_slope_plot(wk_all,sp_all,clst.lbl_neg,icls,all_ages,allgenotype,p_thresh)
        title(sprintf('MSE, d_median = %1.2f',clst.d_neg(icls)),'fontsize',30)
        ylabel('SampEn')
        %ylim([0.7 1.7])
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s%i.svg',DIRFIGURE,'MSE_neg_clust_dot_plot_#',icls))
        print(gcf,'-dpng','-r300','-painters',sprintf('%s%s%i.png',DIRFIGURE,'MSE_neg_clust_dot_plot_#',icls))
    end
end      

%% plot average gMLZ scales

myfigure
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
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
legend('wake','wake 95% CI','sleep','95% sleep','AutoUpdate','off','FontSize',18,'location','northwest')
legend boxoff
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
xlim([0 5])
% ticks = 0:1:5;
% xticks(ticks)
% ticklabels = cell(1,length(ticks));
% for i = 1:length(ticklabels)
%     ticklabels{i} = num2str(2^ticks(i));
% end
% xticklabels(ticklabels)
%title('MSE difference over all time scales (no regression)','fontsize',18)
xlabel('Center frequency (Hz)')
ylabel('LZ complexity'); ylim([0 0.8]); 
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_scales_no_reg'))
print(gcf,'-depsc',sprintf('%s%s.svg',DIRFIGURE,'gMLZ_scales_no_reg'))

% %% plot average gMLZ scales w/TD
% 
% myfigure
% %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% h=subplot(1,1,1);
% hold on
% [~,~,ci] = ttest(wk_scl_lzc');
% plot(log2(LZC_foi(scales)),mean(wk_scl_lzc,2),'r','linewidth',2)
% patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
%     0.3,'EdgeColor','none')
% [~,~,ci] = ttest(TD_lzc_scl');
% plot(log2(LZC_foi(scales)),mean(TD_lzc_scl,2),'color',[0 0.5 0.5],'linewidth',2)
% patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))],[0 0.5 0.5],'facealpha', ...
%     0.3,'EdgeColor','none')
% legend('AS wake','AS wake 95% CI','TD wake','TD 95% wake CI','AutoUpdate','off','FontSize',18,'location','northwest')
% legend boxoff
% set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
% xlim([0 5])
% % ticks = 0:1:5;
% % xticks(ticks)
% % ticklabels = cell(1,length(ticks));
% % for i = 1:length(ticklabels)
% %     ticklabels{i} = num2str(2^ticks(i));
% % end
% % xticklabels(ticklabels)
% % title('MSE difference over all time scales (no regression)','fontsize',18)
% xlabel('Center frequency (Hz)')
% ylabel('LZ complexity'); ylim([0 0.8]); 
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_scales_AS_TD'))
% print(gcf,'-depsc',sprintf('%s%s.svg',DIRFIGURE,'gMLZ_scales_AS_TD'))

%% plot average gMLZ scales ALT X-AXIS

myfigure
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
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
legend('wake','wake 95% CI','sleep','95% sleep','AutoUpdate','off','FontSize',18,'location','northwest')
legend boxoff
pltscl = [1 5 10 15 20];
ticks = log2(LZC_foi(pltscl));
xticks(ticks)
ticklabels = cell(1,length(ticks));
for i = 1:length(ticklabels)
    ticklabels{i} = num2str(pltscl(i));
end
xticklabels(ticklabels)
xlabel('time scale')
xlabel('Center frequency (Hz)')
ylabel('LZ complexity'); ylim([0 0.8]); 
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_scales_no_reg_ALT_XAXIS'))
print(gcf,'-depsc',sprintf('%s%s.svg',DIRFIGURE,'gMLZ_scales_no_reg_ALT_XAXIS'))

%% gMLZ percent change from sleep

myfigure
h=subplot(1,1,1);
hold on
wsdelta = (mean(wk_scl_lzc,2) - mean(sp_scl_lzc,2))./mean(sp_scl_lzc,2).*100;
[~,~,ci] = ttest(wk_scl_lzc,sp_scl_lzc,'dim',2); ci = ci';
plot(log2(LZC_foi(scales)),wsdelta ,'m','linewidth',4)
patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , ...
    [ci(1,:)./mean(sp_scl_lzc,2)'.*100 fliplr(ci(2,:)./mean(sp_scl_lzc,2)'.*100)], 'm','facealpha', ...
    1.0,'EdgeColor','none')
legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',18,'location','northwest')
legend boxoff
plot(log2(LZC_foi(scales)),zeros(1,length(wsdelta)),'-.k')
xlim([0 5])
ylim([0 40])
title('gMLZ difference over all time scales (no regression)','fontsize',18)
xlabel('Center frequency (Hz)')
ylabel('gLZ change from sleep (%)'); ones(1,length(foi));
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_scales_delta_no_reg'))
print(gcf,'-depsc',sprintf('%s%s.svg',DIRFIGURE,'gMLZ_scales_delta_no_reg'))

%% gMLZ raw change from sleep

myfigure
h=subplot(1,1,1);
hold on
wsdelta = mean(wk_scl_lzc,2) - mean(sp_scl_lzc,2);
[~,~,ci] = ttest(wk_scl_lzc,sp_scl_lzc,'dim',2); ci = ci';
plot(log2(LZC_foi(scales)),wsdelta ,'m','linewidth',4)
patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , ...
    [ci(1,:) fliplr(ci(2,:))], 'm','facealpha', ...
    1.0,'EdgeColor','none')
legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',18,'location','northwest')
legend boxoff
plot(log2(LZC_foi(scales)),zeros(1,length(wsdelta)),'-.k')
xlim([0 5])
ylim([-0.2 0.2])
title('gMLZ difference over all time scales (no regression)','fontsize',18)
xlabel('Center frequency (Hz)')
ylabel('gLZ change from sleep'); ones(1,length(foi));
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_scales_delta_no_reg_RAW'))
print(gcf,'-depsc',sprintf('%s%s.svg',DIRFIGURE,'gMLZ_scales_delta_no_reg_RAW'))




%% gMLZ permutation clustering

clst_LZC = RMpermclusttest(wk_lzc_all,sp_lzc_all,Nperm,p_thresh);

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
        myfigure
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
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
        
        myfigure
        plot_topo_AS(sum(clst_LZC.lbl_pos==icls)')
        title(sprintf('gMLZ pos clust, p = %1.3f',clst_LZC.P_val_pos(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(wk_lzc_all,2)])
        mycolorbar
        figname = sprintf('LZC pos clust topo #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
    end
    
    % dot plots for all clusters (even if p > 0.05)
    if any(~isnan(clst_LZC.lbl_pos(:)))    
        dot_slope_plot(wk_lzc_all,sp_lzc_all,clst_LZC.lbl_pos,icls,all_ages,allgenotype,p_thresh)
        title(sprintf('gLZC, d_median = %1.2f',clst_LZC.d_pos(icls)))
        ylabel('gMLZ')
        %ylim([0.14 0.36])
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s%i.svg',DIRFIGURE,'LZC_pos_clust_dot_plot_#',icls))
        print(gcf,'-dpng','-r300','-painters',sprintf('%s%s%i.png',DIRFIGURE,'LZC_pos_clust_dot_plot_#',icls))
    end
end

cnt = 0;
for icls = 1:length(sgcl_neg)
    if sgcl_neg(icls)
        cnt = cnt+1;
        myfigure
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
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
        
        myfigure
        plot_topo_AS(sum(clst_LZC.lbl_neg==icls)')
        title(sprintf('gMLZ neg clust, p = %1.3f',clst_LZC.P_val_neg(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(WP,2)])
        mycolorbar
        figname = sprintf('gMLZ neg clust topo #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
    end
    
    % dot plots for all clusters (even if p > 0.05)
    if any(~isnan(clst_LZC.lbl_neg(:)))            
        dot_slope_plot(wk_lzc_all,sp_lzc_all,clst_LZC.lbl_neg,icls,all_ages,allgenotype,p_thresh)
        ylabel('gMLZ')
        title(sprintf('gLZC, d_median = %1.2f',clst_LZC.d_neg(icls)))
        %ylim([0.14 0.36])
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s%i.svg',DIRFIGURE,'LZC_neg_clust_dot_plot_#',icls))
        print(gcf,'-dpng','-r300','-painters',sprintf('%s%s%i.png',DIRFIGURE,'LZC_neg_clust_dot_plot_#',icls))
    end
end

%% Do scatter plots of age

myfigure
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

myfigure
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

myfigure
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

myfigure
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

myfigure
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

myfigure
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

myfigure
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

myfigure
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

myfigure
plot([Rwk; Rsp]','linewidth',2)
xlabel('Time scale')
ylabel('R^{2}')
legend('Wake','Sleep')
legend boxoff
title('Variance of LZC explained by age')
makefigpretty
figname = 'R2LZCage';
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

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

myfigure
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
    'AutoUpdate','off','FontSize',18)
xlim([0 5])
xlabel('Frequency (Hz)','FontSize',12)
ylabel('Power log_1_0(\muV^{2}/Hz)','FontSize',18)
title('Absolute Power','FontSize',20)
set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
box off
%axis square
legend boxoff
%ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
%annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
    %%'LineStyle','none','FontSize',18);
set(gca,'linewidth',3)
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 22)
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 22)
set(gca, 'TickDir', 'out')
figname = sprintf('PSDs_abs_power_sp_wk');
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
clear wkci spci wk sp

% %% PSD difference plot
% 
% % Doesn't quite match after we convert dB to percent change ... don't use
% % this one 
% 
% myfigure
% h=subplot(1,1,1);
% foi_hd = 2.^[0:0.01:5]; % interpolation resolution for plotting
% wk = 10*log10(squeeze(nanmean(WP)))';
% sp = 10*log10(squeeze(nanmean(SP)))';
% wsd_dB = mean(wk)-mean(sp);
% wsd = (10.^(wsd_dB/10).*100)-100; % convert dB to percent change
% [~,~,ci] = ttest(wk,sp);
% ci = (10.^(ci/10).*100)-100;% convert CI to percentage 
% lbl1 = 'Wake - Sleep';
% plot(log2(foi_hd),interp1(foi,wsd,foi_hd,'spline'),'m','LineWidth',4)
% hold on
% patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ci(1,:),foi_hd,'spline') ...
%     fliplr(interp1(foi,ci(2,:),foi_hd,'spline'))],'m','facealpha',1.0,'EdgeColor','none')
% legend({lbl1,sprintf('%s%s',lbl1,' 95% CI')},'AutoUpdate','off','FontSize',18)
% plot(log2(foi_hd),zeros(1,length(foi_hd)),'k--')
% xlim([0 5])
% xlabel('Frequency (Hz)','FontSize',12)
% ylabel('Power change from sleep (%)'); ylim([-25 25])
% title('DON''T USE THIS ONE','FontSize',20)
% set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
% legend boxoff
% %ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
% %annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
%     %'LineStyle','none','FontSize',18);
% figname = sprintf('PSDs_abs_power_DIFFERNCE');
% makefigpretty
% ylim([-100 100])
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))

%% PSD plots absolute power (linear)

myfigure
h=subplot(1,1,1);
[~,~,wkci] = ttest( squeeze(nanmean(WP))' );
[~,~,spci] = ttest( squeeze(nanmean(SP))' );
foi_hd = 2.^[0:0.01:5]; % interpolation resolution for plotting
wk = squeeze(nanmean(WP))';
sp = squeeze(nanmean(SP))';
wkclr = [1 0 0];
spclr = [0 0 1];
lbl1 = 'Wake';
lbl2 = 'Sleep';
plot(log2(foi_hd),interp1(foi,nanmean(wk)',foi_hd,'spline'),'Color',wkclr,'LineWidth',4)
hold on
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,wkci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,wkci(2,:),foi_hd,'spline'))],wkclr,'facealpha',0.3,'EdgeColor','none')
plot(log2(foi_hd),interp1(foi,nanmean(sp)',foi_hd,'spline'),'Color',spclr,'LineWidth',4)
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,spci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,spci(2,:),foi_hd,'spline'))],spclr,'facealpha',0.3,'EdgeColor','none')
legend({lbl1,sprintf('%s%s',lbl1,' 95% CI'),lbl2,sprintf('%s%s',lbl2,' 95% CI')},...
    'AutoUpdate','off','FontSize',18)
xlim([0 5])
xlabel('Frequency (Hz)','FontSize',12)
ylabel('Power log_1_0(\muV^{2}/Hz)','FontSize',18)
title('Absolute Power','FontSize',20)
set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
box off
%axis square
legend boxoff
%ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
%annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
    %%'LineStyle','none','FontSize',18);
set(gca,'linewidth',3)
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 22)
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 22)
set(gca, 'TickDir', 'out')
figname = sprintf('PSDs_abs_power_sp_wk_LINEAR');
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
clear wkci spci wk sp

%% PSD difference plot (linear power ref'ed to sleep)

% Zero-crossings match with PSD plot (checked 12.03.19)

myfigure
h=subplot(1,1,1);
foi_hd = 2.^[0:0.01:5]; % interpolation resolution for plotting
wk = squeeze(nanmean(WP))';
sp = squeeze(nanmean(SP))';
wsd = (wk-sp)./sp.*100; % percent change
[~,~,ci] = ttest(wsd);
lbl1 = 'Wake - Sleep';
plot(log2(foi_hd),interp1(foi,mean(wsd),foi_hd,'spline'),'m','LineWidth',4)
hold on
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,ci(2,:),foi_hd,'spline'))],'m','facealpha',0.3,'EdgeColor','none')
legend({lbl1,sprintf('%s%s',lbl1,' 95% CI')},'AutoUpdate','off','FontSize',18)
plot(log2(foi_hd),zeros(1,length(foi_hd)),'k--')
xlim([0 5])
xlabel('Frequency (Hz)','FontSize',12)
ylabel('Power change from sleep (%)'); 
title('Absolute Power','FontSize',20)
set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
legend boxoff
ylim([-100 100])
%ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
%annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
    %'LineStyle','none','FontSize',18);
figname = sprintf('PSDs_abs_power_LINEAR_DIFFERNCE');
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))

%% PSD difference plot (linear power ref'ed to sleep) -- RAW VERSION 

% Zero-crossings match with PSD plot (checked 12.03.19)

myfigure
h=subplot(1,1,1);
foi_hd = 2.^[0:0.01:5]; % interpolation resolution for plotting
wk = squeeze(nanmean(WP))';
sp = squeeze(nanmean(SP))';
wsd = log10(wk)-log10(sp); % percent change
[~,~,ci] = ttest(wsd);
lbl1 = 'Wake - Sleep';
plot(log2(foi_hd),interp1(foi,mean(wsd),foi_hd,'spline'),'m','LineWidth',4)
hold on
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,ci(2,:),foi_hd,'spline'))],'m','facealpha',1.0,'EdgeColor','none')
legend({lbl1,sprintf('%s%s',lbl1,' 95% CI')},'AutoUpdate','off','FontSize',18,'location','southeast')
plot(log2(foi_hd),zeros(1,length(foi_hd)),'k--')
xlim([0 5])
xlabel('Frequency (Hz)','FontSize',12)
ylabel('Power change from sleep (log10 (\muV^{2}/log_{2}(Hz)))'); 
title('Absolute Power','FontSize',20)
set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
legend boxoff
ylim([-1 1])
%ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
%annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
    %'LineStyle','none','FontSize',18);
figname = sprintf('PSDs_abs_power_LINEAR_DIFFERNCE_RAW');
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))


% %% PSD difference plot (linear power ref'ed to sleep)
% 
% myfigure
% h=subplot(1,1,1);
% foi_hd = 2.^[0:0.01:5]; % interpolation resolution for plotting
% wk = squeeze(nanmean(WP))';
% sp = squeeze(nanmean(SP))';
% wsd = wk-sp; % percent change
% [~,~,ci] = ttest(wsd);
% lbl1 = 'Wake - Sleep';
% plot(log2(foi_hd),interp1(foi,mean(wsd),foi_hd,'spline'),'m','LineWidth',4)
% hold on
% patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ci(1,:),foi_hd,'spline') ...
%     fliplr(interp1(foi,ci(2,:),foi_hd,'spline'))],'m','facealpha',1.0,'EdgeColor','none')
% legend({lbl1,sprintf('%s%s',lbl1,' 95% CI')},'AutoUpdate','off','FontSize',18,'location','southeast')
% plot(log2(foi_hd),zeros(1,length(foi_hd)),'k--')
% xlim([0 5])
% xlabel('Frequency (Hz)','FontSize',12)
% ylabel('Power change from sleep'); 
% title('Absolute Power','FontSize',20)
% set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
% legend boxoff
% ylim([-1250 100])
% %ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
% %annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
%     %'LineStyle','none','FontSize',18);
% figname = sprintf('PSDs_abs_power_LINEAR_DIFFERNCE_RAW');
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))

%% PSD difference plot (linear power ref'ed to sleep)

% myfigure
% h=subplot(1,1,1);
% foi_hd = 2.^[0:0.01:5]; % interpolation resolution for plotting
% wk = squeeze(nanmean(WP))';
% sp = squeeze(nanmean(SP))';
% wsd = wk-sp; % percent change
% [~,~,ci] = ttest(wsd);
% lbl1 = 'Wake - Sleep';
% plot(log2(foi_hd),interp1(foi,mean(wsd),foi_hd,'spline'),'m','LineWidth',4)
% hold on
% patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ci(1,:),foi_hd,'spline') ...
%     fliplr(interp1(foi,ci(2,:),foi_hd,'spline'))],'m','facealpha',1.0,'EdgeColor','none')
% legend({lbl1,sprintf('%s%s',lbl1,' 95% CI')},'AutoUpdate','off','FontSize',18,'location','southeast')
% plot(log2(foi_hd),zeros(1,length(foi_hd)),'k--')
% xlim([0 5])
% xlabel('Frequency (Hz)','FontSize',12)
% ylabel('Power change from sleep'); 
% title('Absolute Power','FontSize',20)
% set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
% legend boxoff
% ylim([-1250 100])
% %ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
% %annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
%     %'LineStyle','none','FontSize',18);
% figname = sprintf('PSDs_abs_power_LINEAR_DIFFERNCE_RAW');
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))

%% PSD difference plot (linear power referenced to awake)

%  WARNING: Zero-crossings DON'T match with PSD plot (checked 12.03.19)

myfigure
h=subplot(1,1,1);
foi_hd = 2.^[0:0.01:5]; % interpolation resolution for plotting
wk = squeeze(nanmean(WP))';
sp = squeeze(nanmean(SP))';
swd = (sp-wk)./wk.*100; % percent change
swd2 = mean(sp-wk)./wk.*100; % alternative formula (zero crossings will match but now CIs look buggy)
[~,~,ci] = ttest(swd);
lbl1 = 'Sleep - Wake';
plot(log2(foi_hd),interp1(foi,mean(swd),foi_hd,'spline'),'m','LineWidth',4)
hold on
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,ci(2,:),foi_hd,'spline'))],'m','facealpha',0.3,'EdgeColor','none')
legend({lbl1,sprintf('%s%s',lbl1,' 95% CI')},'AutoUpdate','off','FontSize',18)
plot(log2(foi_hd),zeros(1,length(foi_hd)),'k--')
xlim([0 5])
xlabel('Frequency (Hz)','FontSize',12)
ylabel('Power change from wake (%)'); 
title('Absolute Power','FontSize',20)
set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
legend boxoff
ylim([-50 300])
%ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
%annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
    %'LineStyle','none','FontSize',18);
figname = sprintf('PSDs_abs_power_LINEAR_WK_REF');
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))


%% PSD difference plot (linear power referenced to awake) -- RAW VERSION

% Zero-crossings match with PSD plot (checked 12.03.19)

myfigure
h=subplot(1,1,1);
foi_hd = 2.^[0:0.01:5]; % interpolation resolution for plotting
wk = squeeze(nanmean(WP))';
sp = squeeze(nanmean(SP))';
swd = log10(sp)-log10(wk); % percent change
[~,~,ci] = ttest(swd);
lbl1 = 'Sleep - Wake';
plot(log2(foi_hd),interp1(foi,mean(swd),foi_hd,'spline'),'m','LineWidth',4)
hold on
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,ci(2,:),foi_hd,'spline'))],'m','facealpha',0.3,'EdgeColor','none')
legend({lbl1,sprintf('%s%s',lbl1,' 95% CI')},'AutoUpdate','off','FontSize',18)
plot(log2(foi_hd),zeros(1,length(foi_hd)),'k--')
xlim([0 5])
xlabel('Frequency (Hz)','FontSize',12)
ylabel('Power change from wake (log10(\muV^{2}/log_{2}(Hz)))'); 
title('Absolute Power','FontSize',20)
set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
legend boxoff
ylim([-1 1])
%ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
%annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
    %'LineStyle','none','FontSize',18);
figname = sprintf('PSDs_abs_power_LINEAR_WK_REF_RAW');
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))

%% power permutation clustering

clst_pow = RMpermclusttest(log10(WP),log10(SP),Nperm,p_thresh);

sgcl_pos = [];
sgcl_neg = [];

if sum(clst_pow.n_pos > 0) > 0
    sgcl_pos = clst_pow.P_val_pos < 0.05;
end

if sum(clst_pow.n_neg > 0) > 0
    sgcl_neg = clst_pow.P_val_neg < 0.05;
end

cnt = 0;
for icls = 1:length(sgcl_pos)
    if sgcl_pos(icls)
        cnt = cnt + 1;
        myfigure
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
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
        
        myfigure
        plot_topo_AS(sum(clst_pow.lbl_pos==icls)')
        title(sprintf('power pos clust, p = %1.3f',clst_pow.P_val_pos(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(WP,2)])
        mycolorbar
        figname = sprintf('Power pos clust toop #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
    end
    
    % dot plots for all clusters (even if p > 0.05)
    if any(~isnan(clst_pow.lbl_pos(:)))            
        dot_slope_plot(WP,SP,clst_pow.lbl_pos,icls,all_ages,allgenotype,p_thresh)
        title(sprintf('Power, d_median = %1.2f',clst_pow.d_pos(icls)))
        ylabel('Power log_10(\muV^{2}/log(Hz))')
        %ylim([-1 5])
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s%i.svg',DIRFIGURE,'Power_pos_clust_dot_plot_#',icls))
        print(gcf,'-dpng','-r300','-painters',sprintf('%s%s%i.png',DIRFIGURE,'Power_pos_clust_dot_plot_#',icls))
    end
end

cnt = 0;
for icls = 1:length(sgcl_neg)
    if sgcl_neg(icls)
        cnt = cnt + 1;
        myfigure
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
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
        
        myfigure
        plot_topo_AS(sum(clst_pow.lbl_neg==icls)')
        title(sprintf('power neg clust, p = %1.3f',clst_pow.P_val_neg(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(WP,2)])
        mycolorbar
        figname = sprintf('Power neg clust toop #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
    end
    
    % dot plots for all clusters (even if p > 0.05)
    if any(~isnan(clst_pow.lbl_neg(:)))            
        dot_slope_plot(WP,SP,clst_pow.lbl_neg,icls,all_ages,allgenotype,p_thresh)
        title(sprintf('Power, d_median = %1.2f',clst_pow.d_neg(icls)))
        ylabel('Power log_10(\muV^{2}/log(Hz))')
        %ylim([-1 5])
        print(gcf,'-dsvg','-r600','-painters',sprintf('%s%s%i.svg',DIRFIGURE,'Power_neg_clust_dot_plot_#',icls))
        print(gcf,'-dpng','-r600','-painters',sprintf('%s%s%i.png',DIRFIGURE,'Power_neg_clust_dot_plot_#',icls))
    end
end

%% Power topoplots for reviewer

delta = foi >= 1 & foi < 4;
theta = foi >= 4 & foi < 8;
alpha = foi >= 8 & foi < 12;
beta  = foi >= 12 & foi < 25;
gamma = foi >= 25 & foi < 32;

dband_wk = mean(log10(WP(:,delta,:)),3); % average participants after log-scaling
plot_dw = log10(trapz(foi(delta), 10.^dband_wk'))'; % back transform, take integral, log-scale again

tband_wk = mean(log10(WP(:,theta,:)),3); % average participants after log-scaling
plot_tw = log10(trapz(foi(theta), 10.^tband_wk'))'; % back transform, take integral, log-scale again

aband_wk = mean(log10(WP(:,alpha,:)),3); % average participants after log-scaling
plot_aw = log10(trapz(foi(alpha), 10.^aband_wk'))'; % back transform, take integral, log-scale again

bband_wk = mean(log10(WP(:,beta,:)),3); % average participants after log-scaling
plot_bw = log10(trapz(foi(beta), 10.^bband_wk'))'; % back transform, take integral, log-scale again

gband_wk = mean(log10(WP(:,gamma,:)),3); % average participants after log-scaling
plot_gw = log10(trapz(foi(gamma), 10.^gband_wk'))'; % back transform, take integral, log-scale again

dband_sp = mean(log10(SP(:,delta,:)),3); % average participants after log-scaling
plot_ds = log10(trapz(foi(delta), 10.^dband_sp'))'; % back transform, take integral, log-scale again

tband_sp = mean(log10(SP(:,theta,:)),3); % average participants after log-scaling
plot_ts = log10(trapz(foi(theta), 10.^tband_sp'))'; % back transform, take integral, log-scale again

aband_sp = mean(log10(SP(:,alpha,:)),3); % average participants after log-scaling
plot_as = log10(trapz(foi(alpha), 10.^aband_sp'))'; % back transform, take integral, log-scale again

bband_sp = mean(log10(SP(:,beta,:)),3); % average participants after log-scaling
plot_bs = log10(trapz(foi(beta), 10.^bband_sp'))'; % back transform, take integral, log-scale again

gband_sp = mean(log10(SP(:,gamma,:)),3); % average participants after log-scaling
plot_gs = log10(trapz(foi(gamma), 10.^gband_sp'))'; % back transform, take integral, log-scale again


d1 = round(min([plot_dw; plot_ds]),1); 
d2 = round(max([plot_dw; plot_ds]),1);

myfigure
plot_topo_AS(plot_dw)
colormap jet
mycolorbar
caxis([d1 d2])
title('Delta power (1 - 4 Hz)','fontsize',20)
figname = 'AWAKE_delta_power_topo';
print(gcf,'-dsvg','-r300',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300',sprintf('%s%s.png',DIRFIGURE,figname))

myfigure
plot_topo_AS(plot_ds)
colormap jet
mycolorbar
caxis([d1 d2])
title('Delta power (1 - 4 Hz)','fontsize',20)
figname = 'SLEEP_delta_power_topo';
print(gcf,'-dsvg','-r300',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300',sprintf('%s%s.png',DIRFIGURE,figname))

t1 = round(min([plot_tw; plot_ts]),1); 
t2 = round(max([plot_tw; plot_ts]),1);

myfigure
plot_topo_AS(plot_tw)
colormap jet
mycolorbar
caxis([t1 t2])
title('Theta power (1 - 4 Hz)','fontsize',20)
figname = 'AWAKE_theta_power_topo';
print(gcf,'-dsvg','-r300',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300',sprintf('%s%s.png',DIRFIGURE,figname))

myfigure
plot_topo_AS(plot_ts)
colormap jet
mycolorbar
caxis([t1 t2])
title('Theta power (4 - 8 Hz)','fontsize',20)
figname = 'SLEEP_theta_power_topo';
print(gcf,'-dsvg','-r300',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300',sprintf('%s%s.png',DIRFIGURE,figname))


a1 = round(min([plot_aw; plot_as]),1); 
a2 = round(max([plot_aw; plot_as]),1);

myfigure
plot_topo_AS(plot_aw)
colormap jet
mycolorbar
caxis([a1 a2])
title('Alpha power (8 - 12 Hz)','fontsize',20)
figname = 'AWAKE_alpha_power_topo';
print(gcf,'-dsvg','-r300',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300',sprintf('%s%s.png',DIRFIGURE,figname))

myfigure
plot_topo_AS(plot_as)
colormap jet
mycolorbar
caxis([a1 a2])
title('Alpha power (8 - 12 Hz)','fontsize',20)
figname = 'SLEEP_alpha_power_topo';
print(gcf,'-dsvg','-r300',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300',sprintf('%s%s.png',DIRFIGURE,figname))


b1 = round(min([plot_bw; plot_bs]),1); 
b2 = round(max([plot_bw; plot_bs]),1);

myfigure
plot_topo_AS(plot_bw)
colormap jet
mycolorbar
caxis([b1 b2])
title('beta power (12 - 25 Hz)','fontsize',20)
figname = 'AWAKE_beta_power_topo';
print(gcf,'-dsvg','-r300',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300',sprintf('%s%s.png',DIRFIGURE,figname))

myfigure
plot_topo_AS(plot_bs)
colormap jet
mycolorbar
caxis([b1 b2])
title('beta power (12 - 25 Hz)','fontsize',20)
figname = 'SLEEP_beta_power_topo';
print(gcf,'-dsvg','-r300',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300',sprintf('%s%s.png',DIRFIGURE,figname))

g1 = round(min([plot_gw; plot_gs]),1); 
g2 = round(max([plot_gw; plot_gs]),1);

myfigure
plot_topo_AS(plot_gw)
colormap jet
mycolorbar
caxis([g1 g2])
title('Gamma power (1 - 4 Hz)','fontsize',20)
figname = 'AWAKE_gamma_power_topo';
print(gcf,'-dsvg','-r300',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300',sprintf('%s%s.png',DIRFIGURE,figname))

myfigure
plot_topo_AS(plot_gs)
colormap jet
mycolorbar
caxis([g1 g2])
title('Gamma power (1 - 4 Hz)','fontsize',20)
figname = 'SLEEP_gamma_power_topo';
print(gcf,'-dsvg','-r300',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300',sprintf('%s%s.png',DIRFIGURE,figname))

fprintf('Delta min - max, %2.1f - %2.1f log10(uV^2)\n',d1,d2)
fprintf('Theta min - max, %2.1f - %2.1f log10(uV^2)\n',t1,t2)
fprintf('Alpha min - max, %2.1f - %2.1f log10(uV^2)\n',a1,a2)
fprintf('Beta  min - max, %2.1f - %2.1f log10(uV^2)\n',b1,b2)
fprintf('Gamma min - max, %2.1f - %2.1f log10(uV^2)\n',g1,g2)

%% Topoplots for wake - sleep

myfigure
plot_topo_AS(plot_dw - plot_ds)
colormap jet
mycolorbar
caxis([-0.2 0.2])
%caxis([0 d2])
title('Delta power (1 - 4 Hz)','fontsize',20)
figname = 'wake_minus_sleep_delta_power_topo';
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300','-painters',sprintf('%s%s.png',DIRFIGURE,figname))

myfigure
plot_topo_AS(plot_tw - plot_ts)
colormap jet
mycolorbar
caxis([-0.2 0.2])
%caxis([0 t2])
title('Theta power (4 - 8 Hz)','fontsize',20)
figname = 'wake_minus_sleep_theta_power_topo';
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300','-painters',sprintf('%s%s.png',DIRFIGURE,figname))

myfigure
plot_topo_AS(plot_aw - plot_as)
colormap jet
mycolorbar
caxis([-0.2 0.2])
%caxis([0 a2])
title('Alpha power (8 - 12 Hz)','fontsize',20)
figname = 'wake_minus_sleep_alpha_power_topo';
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300','-painters',sprintf('%s%s.png',DIRFIGURE,figname))

myfigure
plot_topo_AS(plot_bw - plot_bs)
colormap jet
mycolorbar
caxis([-0.2 0.2])
%caxis([0 b2])
title('Beta power (12 - 25 Hz)','fontsize',20)
figname = 'wake_minus_sleep_beta_power_topo';
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300','-painters',sprintf('%s%s.png',DIRFIGURE,figname))

myfigure
plot_topo_AS(plot_gw - plot_gs)
colormap jet
mycolorbar
caxis([-0.2 0.2])
%caxis([round(min(plotme),1) round(max(plotme),1)])
title('Gamma power (25 - 32 Hz)','fontsize',20)
figname = 'wake_minus_sleep_gamma_power_topo';
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300','-painters',sprintf('%s%s.png',DIRFIGURE,figname))

%% Now do again for t-stats

dband_wk = WP(:,delta,:);
dband_sp = SP(:,delta,:);
dw_tot = nan(size(WP,1),size(WP,3));
ds_tot = nan(size(SP,1),size(SP,3));

for isb = 1:size(WP,3)
    dw_tot(:,isb) = trapz(foi(delta), squeeze(dband_wk(:,:,isb))');
    ds_tot(:,isb) = trapz(foi(delta), squeeze(dband_sp(:,:,isb))');
end

[~,~,~,stat] = ttest(log10(dw_tot),log10(ds_tot),'dim',2);
myfigure
plot_topo_AS(stat.tstat)
colormap jet
mycolorbar
caxis([-4 4])
title('Delta_t_stats','fontsize',20)
figname = 'Delta_t_stats';
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300','-painters',sprintf('%s%s.png',DIRFIGURE,figname))

tband_wk = WP(:,theta,:);
tband_sp = SP(:,theta,:);
tw_tot = nan(size(WP,1),size(WP,3));
ts_tot = nan(size(SP,1),size(SP,3));

for isb = 1:size(WP,3)
    tw_tot(:,isb) = trapz(foi(theta), squeeze(tband_wk(:,:,isb))');
    ts_tot(:,isb) = trapz(foi(theta), squeeze(tband_sp(:,:,isb))');
end

[~,~,~,stat] = ttest(log10(tw_tot),log10(ts_tot),'dim',2);
myfigure
plot_topo_AS(stat.tstat)
colormap jet
mycolorbar
caxis([-4 4])
title('theta_t_stats','fontsize',20)
figname = 'theta_t_stats';
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300','-painters',sprintf('%s%s.png',DIRFIGURE,figname))

aband_wk = WP(:,alpha,:);
aband_sp = SP(:,alpha,:);
aw_tot = nan(size(WP,1),size(WP,3));
as_tot = nan(size(SP,1),size(SP,3));

for isb = 1:size(WP,3)
    aw_tot(:,isb) = trapz(foi(alpha), squeeze(aband_wk(:,:,isb))');
    as_tot(:,isb) = trapz(foi(alpha), squeeze(aband_sp(:,:,isb))');
end

[~,~,~,stat] = ttest(log10(aw_tot),log10(as_tot),'dim',2);
myfigure
plot_topo_AS(stat.tstat)
colormap jet
mycolorbar
caxis([-4 4])
title('alpha_t_stats','fontsize',20)
figname = 'alpha_t_stats';
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300','-painters',sprintf('%s%s.png',DIRFIGURE,figname))

bband_wk = WP(:,beta,:);
bband_sp = SP(:,beta,:);
bw_tot = nan(size(WP,1),size(WP,3));
bs_tot = nan(size(SP,1),size(SP,3));

for isb = 1:size(WP,3)
    bw_tot(:,isb) = trapz(foi(beta), squeeze(bband_wk(:,:,isb))');
    bs_tot(:,isb) = trapz(foi(beta), squeeze(bband_sp(:,:,isb))');
end

[~,~,~,stat] = ttest(log10(bw_tot),log10(bs_tot),'dim',2);
myfigure
plot_topo_AS(stat.tstat)
colormap jet
mycolorbar
caxis([-4 4])
title('beta_t_stats','fontsize',20)
figname = 'beta_t_stats';
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300','-painters',sprintf('%s%s.png',DIRFIGURE,figname))


gband_wk = WP(:,gamma,:);
gband_sp = SP(:,gamma,:);
gw_tot = nan(size(WP,1),size(WP,3));
gs_tot = nan(size(SP,1),size(SP,3));

for isb = 1:size(WP,3)
    gw_tot(:,isb) = trapz(foi(gamma), squeeze(gband_wk(:,:,isb))');
    gs_tot(:,isb) = trapz(foi(gamma), squeeze(gband_sp(:,:,isb))');
end

[~,~,~,stat] = ttest(log10(gw_tot),log10(gs_tot),'dim',2);
myfigure
plot_topo_AS(stat.tstat)
colormap jet
mycolorbar
caxis([-4 4])
title('gamma_t_stats','fontsize',20)
figname = 'gamma_t_stats';
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
title('')
%colorbar off
print(gcf,'-dpng','-r300','-painters',sprintf('%s%s.png',DIRFIGURE,figname))



%%
    
myfigure
for ifoi = 1:1:length(foi)
    subplot(7,6,ifoi)
    plot_topo_AS_classic( mean(log10(WP(:,ifoi,:)),3) )
    caxis([-1.5 3])
    colormap jet
    title(sprintf('f = %1.2f Hz',foi(ifoi)))
end

figname = 'AWAKE_all_freq_topos';
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
print(gcf,'-dpng','-r300','-painters',sprintf('%s%s.png',DIRFIGURE,figname))

myfigure
for ifoi = 1:length(foi)
    subplot(7,6,ifoi)
    plot_topo_AS_classic( mean(log10(SP(:,ifoi,:)),3) )
    caxis([-1.5 3])
    colormap jet
    title(sprintf('f = %1.2f Hz',foi(ifoi)))
end

figname = 'SLEEP_all_freq_topos';
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
print(gcf,'-dpng','-r300','-painters',sprintf('%s%s.png',DIRFIGURE,figname))

% Another way of doing it -- maybe this is more what the reviewer wanted?

WPavg = squeeze(mean(WP,3));
NDX = clst_pow.lbl_pos';
NDX(isnan(NDX)) = 0;
clstavg = nan(size(WPavg,1),1);
for ich = 1:size(clstavg,1)
    if sum(NDX(ich,:)) > 0
        clstavg(ich) = mean(WPavg(ich,logical(NDX(ich,:))));
    else
        ich
        clstavg(ich) = 0;
    end
end

myfigure
plot_topo_AS(clstavg)


%% power Correlations

logWP = log10(squeeze(mean(WP)));
logSP = log10(squeeze(mean(SP)));

% categorical variables

r_geno_wk = mes(logWP(:,logical(allgenotype))',logWP(:,~logical(allgenotype))','requiv');
r_geno_sp = mes(logSP(:,logical(allgenotype))',logSP(:,~logical(allgenotype))','requiv');

r_site_wk = mes(logWP(:,logical(allsite))',logWP(:,~logical(allsite))','requiv');
r_site_sp = mes(logSP(:,logical(allsite))',logSP(:,~logical(allsite))','requiv');

% Take out subjects without sleep pattern data
logWP2 = logWP;
logWP2(:,isnan(STN_TF)) = [];
logSP2 = logSP;
logSP2(:,isnan(STN_TF)) = [];

STN_TF2 = STN_TF;
STN_TF2(isnan(STN_TF2)) = [];

r_STN_wk = mes(logWP2(:,logical(STN_TF2))',logWP2(:,~logical(STN_TF2))','requiv','missVal','listwise');
r_STN_sp = mes(logSP2(:,logical(STN_TF2))',logSP2(:,~logical(STN_TF2))','requiv','missVal','listwise');

% continuous variables 
[r_age_power_wk,p_age_power_wk] = corr(log2(all_ages)',logWP');
[r_age_power_sp,p_age_power_sp] = corr(log2(all_ages)',logSP');

[r_pow_cog_wk,p_pow_cog_wk] = corr(cog_score',logWP','rows','complete');
[r_pow_cog_sp,p_pow_cog_sp] = corr(cog_score',logSP','rows','complete');

[r_pow_lng_wk,p_pow_lng_wk] = corr(lng_score',logWP','rows','complete');
[r_pow_lng_sp,p_pow_lng_sp] = corr(lng_score',logSP','rows','complete');

% significance

PVALS = [r_geno_wk.t.p; r_geno_sp.t.p; r_site_wk.t.p; ...
    r_site_sp.t.p; r_STN_wk.t.p; r_STN_sp.t.p; p_age_power_wk; p_age_power_sp; ...
    p_pow_cog_wk; p_pow_cog_sp; p_pow_lng_wk; p_pow_lng_sp];

RVALS = [r_geno_wk.requiv; r_geno_sp.requiv; r_site_wk.requiv; ...
    r_site_sp.requiv; r_STN_wk.requiv; r_STN_sp.requiv; r_age_power_wk; r_age_power_sp; ...
    r_pow_cog_wk; r_pow_cog_sp; r_pow_lng_wk; r_pow_lng_sp];

% FDR correction

QVALS = mafdr(PVALS(:),'BHFDR',true);
QVALS = reshape(QVALS,size(PVALS,1),size(PVALS,2));

% color order for plotting
co = [0 0.7 0.7;
      0 0.7 0;
      0.850 0.325 0.098;
      0 0 1;
      1 0 1;
      0 0 0];
set(groot,'defaultAxesColorOrder',co)

myfigure
h=subplot(1,1,1);
hold on
% awake
plot(log2(foi),-log10(QVALS(1:2:end,:)),'linewidth',2)
legend({'Genotype','Site','STN','log_{2}(Age)','Cognition', ...
    'Language'},'location','eastoutside','autoupdate','off')
% asleep
plot(log2(foi),-log10(QVALS(2:2:end,:)),':','linewidth',2)
plot(log2(foi),ones(1,length(foi)).*-log10(0.05),'--','color',[0.3 0.3 0.3],'linewidth',2)
legend boxoff
xlabel('Frequency (Hz)')
ylabel('-log_{10}(p)'); ylim([0 0.8])
XTICK = [0:1:5];
set(h,'XTick',XTICK,'XTickLabel',2.^XTICK)
xlim([0 5])
ylim([0 4])
makefigpretty
% SAVE FIGURE
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'power_corr_pval'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'power_corr_pval'))

% R^2 x sign(r)

myfigure
h=subplot(1,1,1);
hold on
% awake
plot(log2(foi),RVALS(1:2:end,:).^2.*sign(RVALS(1:2:end,:)),'linewidth',2)
legend({'Genotype','Site','STN','log_{2}(Age)','Cognition', ...
    'Language'},'location','eastoutside','autoupdate','off')
% asleep
plot(log2(foi),RVALS(2:2:end,:).^2.*sign(RVALS(2:2:end,:)),':','linewidth',2)
%plot(log2(foi),zeros(1,20),'color',[0.3 0.3 0.3],'linewidth',2)
legend boxoff
ylabel('R^{2} x sign(r)')
xlabel('Frequency (Hz)')
ylabel('R^{2} x sign(r)')
XTICK = [0:1:5];
set(h,'XTick',XTICK,'XTickLabel',2.^XTICK)
xlim([0 5])
ylim([-1 1])
makefigpretty
% SAVE FIGURE
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'power_correlations'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'power_correlations'))

%% Figure for reviewers (POWER correlations with age, sleep, genotype)

% significance

PVALS = [r_geno_wk.t.p; r_geno_sp.t.p; r_STN_wk.t.p; r_STN_sp.t.p; p_age_power_wk; p_age_power_sp];
RVALS = [r_geno_wk.requiv; r_geno_sp.requiv; r_STN_wk.requiv; r_STN_sp.requiv; r_age_power_wk; r_age_power_sp];

% FDR correction

QVALS = mafdr(PVALS(:),'BHFDR',true);
QVALS = reshape(QVALS,size(PVALS,1),size(PVALS,2));

% color order for plotting
co = [0 0 1;
      1 0 1;
      0 0 0];
set(groot,'defaultAxesColorOrder',co)

myfigure
h=subplot(1,1,1);
hold on
% awake
plot(log2(foi),-log10(QVALS(1:2:end,:)),'linewidth',2)
legend({'Genotype','STN','log_{2}(Age)'},'location','eastoutside','autoupdate','off')
% asleep
plot(log2(foi),-log10(QVALS(2:2:end,:)),':','linewidth',2)
plot(log2(foi),ones(1,length(foi)).*-log10(0.05),'--','color',[0.3 0.3 0.3],'linewidth',2)
legend boxoff
xlabel('Frequency (Hz)')
ylabel('-log_{10}(p)'); ylim([0 0.8])
XTICK = [0:1:5];
set(h,'XTick',XTICK,'XTickLabel',2.^XTICK)
xlim([0 5])
ylim([0 4])
makefigpretty
% SAVE FIGURE
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'power_corr_pval'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'power_corr_pval'))

% R^2 x sign(r)

myfigure
h=subplot(1,1,1);
hold on
% awake
plot(log2(foi),RVALS(1:2:end,:).^2.*sign(RVALS(1:2:end,:)),'linewidth',2)
legend({'Genotype','STN','log_{2}(Age)'},'location','eastoutside','autoupdate','off')
% asleep
plot(log2(foi),RVALS(2:2:end,:).^2.*sign(RVALS(2:2:end,:)),':','linewidth',2)
%plot(log2(foi),zeros(1,20),'color',[0.3 0.3 0.3],'linewidth',2)
legend boxoff
ylabel('R^{2} x sign(r)')
xlabel('Frequency (Hz)')
ylabel('R^{2} x sign(r)')
XTICK = [0:1:5];
set(h,'XTick',XTICK,'XTickLabel',2.^XTICK)
xlim([0 5])
ylim([-1 1])
makefigpretty
% SAVE FIGURE
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'power_correlations'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'power_correlations'))

%% figure for reviewers (mMSE correlations)

mse_wk_chavg = squeeze(mean(wk_all));
mse_sp_chavg  = squeeze(mean(sp_all));

% categorical variables

r_geno_wk = mes(mse_wk_chavg(:,logical(allgenotype))',mse_wk_chavg(:,~logical(allgenotype))','requiv');
r_geno_sp = mes(mse_sp_chavg(:,logical(allgenotype))',mse_sp_chavg(:,~logical(allgenotype))','requiv');

r_site_wk = mes(mse_wk_chavg(:,logical(allsite))',mse_wk_chavg(:,~logical(allsite))','requiv');
r_site_sp = mes(mse_sp_chavg(:,logical(allsite))',mse_sp_chavg(:,~logical(allsite))','requiv');

% Take out subjects without regular sleep pattern data
mse_wk_chavg2 = mse_wk_chavg;
mse_wk_chavg2(:,isnan(STN_TF)) = [];
mse_sp_chavg2 = mse_sp_chavg;
mse_sp_chavg2(:,isnan(STN_TF)) = [];

STN_TF2 = STN_TF;
STN_TF2(isnan(STN_TF2)) = [];

r_STN_wk = mes(mse_wk_chavg2(:,logical(STN_TF2))',mse_wk_chavg2(:,~logical(STN_TF2))','requiv','missVal','listwise');
r_STN_sp = mes(mse_sp_chavg2(:,logical(STN_TF2))',mse_sp_chavg2(:,~logical(STN_TF2))','requiv','missVal','listwise');

% continuous variables 
[r_age_power_wk,p_age_power_wk] = corr(log2(all_ages)',mse_wk_chavg');
[r_age_power_sp,p_age_power_sp] = corr(log2(all_ages)',mse_sp_chavg');

[r_pow_cog_wk,p_pow_cog_wk] = corr(cog_score',mse_wk_chavg','rows','complete');
[r_pow_cog_sp,p_pow_cog_sp] = corr(cog_score',mse_sp_chavg','rows','complete');

[r_pow_lng_wk,p_pow_lng_wk] = corr(lng_score',mse_wk_chavg','rows','complete');
[r_pow_lng_sp,p_pow_lng_sp] = corr(lng_score',mse_sp_chavg','rows','complete');

% significance

PVALS = [r_geno_wk.t.p; r_geno_sp.t.p; r_STN_wk.t.p; r_STN_sp.t.p; p_age_power_wk; p_age_power_sp];
RVALS = [r_geno_wk.requiv; r_geno_sp.requiv; r_STN_wk.requiv; r_STN_sp.requiv; r_age_power_wk; r_age_power_sp];

% FDR correction

QVALS = mafdr(PVALS(:),'BHFDR',true);
QVALS = reshape(QVALS,size(PVALS,1),size(PVALS,2));

% color order for plotting
co = [0 0 1;
      1 0 1;
      0 0 0];
set(groot,'defaultAxesColorOrder',co)

% -log10(p-value)

myfigure
h=subplot(1,1,1);
hold on
% awake
plot(log2(nyqfoi),-log10(QVALS(1:2:end,:)),'linewidth',2)
legend({'Genotype','STN','log_{2}(Age)'},'location','eastoutside','autoupdate','off')
% asleep
plot(log2(nyqfoi),-log10(QVALS(2:2:end,:)),':','linewidth',2)
plot(log2(nyqfoi),ones(1,length(nyqfoi)).*-log10(0.05),'--','color',[0.3 0.3 0.3],'linewidth',2)
legend boxoff
xlabel('Frequency (Hz)')
ylabel('-log_{10}(p)'); ylim([0 0.8])
%XTICK = [0:1:5];
%set(h,'XTick',XTICK,'XTickLabel',2.^XTICK)
ylim([0 4])
makefigpretty
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
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'MSE_corr_pval'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'MSE_corr_pval'))

%% R^2 x sign(r)

myfigure
h=subplot(1,1,1);
hold on
% awake
plot(log2(nyqfoi),RVALS(1:2:end,:).^2.*sign(RVALS(1:2:end,:)),'linewidth',2)
legend({'Genotype','STN','log_{2}(Age)'},'location','eastoutside','autoupdate','off')
% asleep
plot(log2(nyqfoi),RVALS(2:2:end,:).^2.*sign(RVALS(2:2:end,:)),':','linewidth',2)
plot(log2(nyqfoi),zeros(1,20),'color',[0.3 0.3 0.3],'linewidth',2)
legend boxoff
ylabel('R^{2} x sign(r)')
xlabel('Frequency (Hz)')
ylabel('R^{2} x sign(r)')
%set(h,'XTick',XTICK,'XTickLabel',2.^XTICK)
ylim([-1 1])
makefigpretty
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
ylabel('R^{2} x sign(r)')
makefigpretty
ylim([-1 1])
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
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'MSE_correlations'))


%% figure for reviewers (gLZC correlations)

lzc_wk_chavg = squeeze(mean(wk_lzc_all));
lzc_sp_chavg  = squeeze(mean(sp_lzc_all));

% categorical variables

r_geno_wk = mes(lzc_wk_chavg(:,logical(allgenotype))',lzc_wk_chavg(:,~logical(allgenotype))','requiv');
r_geno_sp = mes(lzc_sp_chavg(:,logical(allgenotype))',lzc_sp_chavg(:,~logical(allgenotype))','requiv');

r_site_wk = mes(lzc_wk_chavg(:,logical(allsite))',lzc_wk_chavg(:,~logical(allsite))','requiv');
r_site_sp = mes(lzc_sp_chavg(:,logical(allsite))',lzc_sp_chavg(:,~logical(allsite))','requiv');

% Take out subjects without regular sleep pattern data
lzc_wk_chavg2 = lzc_wk_chavg;
lzc_wk_chavg2(:,isnan(STN_TF)) = [];
lzc_sp_chavg2 = lzc_sp_chavg;
lzc_sp_chavg2(:,isnan(STN_TF)) = [];

STN_TF2 = STN_TF;
STN_TF2(isnan(STN_TF2)) = [];

r_STN_wk = mes(lzc_wk_chavg2(:,logical(STN_TF2))',lzc_wk_chavg2(:,~logical(STN_TF2))','requiv','missVal','listwise');
r_STN_sp = mes(lzc_sp_chavg2(:,logical(STN_TF2))',lzc_sp_chavg2(:,~logical(STN_TF2))','requiv','missVal','listwise');

% continuous variables 
[r_age_power_wk,p_age_power_wk] = corr(log2(all_ages)',lzc_wk_chavg');
[r_age_power_sp,p_age_power_sp] = corr(log2(all_ages)',lzc_sp_chavg');

[r_pow_cog_wk,p_pow_cog_wk] = corr(cog_score',lzc_wk_chavg','rows','complete');
[r_pow_cog_sp,p_pow_cog_sp] = corr(cog_score',lzc_sp_chavg','rows','complete');

[r_pow_lng_wk,p_pow_lng_wk] = corr(lng_score',lzc_wk_chavg','rows','complete');
[r_pow_lng_sp,p_pow_lng_sp] = corr(lng_score',lzc_sp_chavg','rows','complete');

% significance

PVALS = [r_geno_wk.t.p; r_geno_sp.t.p; r_STN_wk.t.p; r_STN_sp.t.p; p_age_power_wk; p_age_power_sp];
RVALS = [r_geno_wk.requiv; r_geno_sp.requiv; r_STN_wk.requiv; r_STN_sp.requiv; r_age_power_wk; r_age_power_sp];

% FDR correction

QVALS = mafdr(PVALS(:),'BHFDR',true);
QVALS = reshape(QVALS,size(PVALS,1),size(PVALS,2));

% color order for plotting
co = [0 0 1;
      1 0 1;
      0 0 0];
set(groot,'defaultAxesColorOrder',co)

% -log10(p-value)

myfigure
h=subplot(1,1,1);
hold on
% awake
plot(log2(LZC_foi),-log10(QVALS(1:2:end,:)),'linewidth',2)
legend({'Genotype','STN','log_{2}(Age)'},'location','eastoutside','autoupdate','off')
% asleep
plot(log2(LZC_foi),-log10(QVALS(2:2:end,:)),':','linewidth',2)
plot(log2(LZC_foi),ones(1,length(LZC_foi)).*-log10(0.05),'--','color',[0.3 0.3 0.3],'linewidth',2)
legend boxoff
xlabel('Frequency (Hz)')
ylabel('-log_{10}(p)'); ylim([0 0.8])
%XTICK = [0:1:5];
%set(h,'XTick',XTICK,'XTickLabel',2.^XTICK)
ylim([0 4])
makefigpretty
xlabel('Center frequency (Hz)')
makefigpretty
xlim([0 5])
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))

% SAVE FIGURE
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'LZC_corr_pval'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'LZC_corr_pval'))

% R^2 x sign(r)

myfigure
h=subplot(1,1,1);
hold on
% awake
plot(log2(LZC_foi),RVALS(1:2:end,:).^2.*sign(RVALS(1:2:end,:)),'linewidth',2)
legend({'Genotype','STN','log_{2}(Age)'},'location','eastoutside','autoupdate','off')
% asleep
plot(log2(LZC_foi),RVALS(2:2:end,:).^2.*sign(RVALS(2:2:end,:)),':','linewidth',2)
plot(log2(LZC_foi),zeros(1,20),'color',[0.3 0.3 0.3],'linewidth',2)
legend boxoff
ylabel('R^{2} x sign(r)')
xlabel('Frequency (Hz)')
ylabel('R^{2} x sign(r)')
%set(h,'XTick',XTICK,'XTickLabel',2.^XTICK)
ylim([-1 1])
makefigpretty
xlabel('Center frequency (Hz)')
makefigpretty
xlim([0 5])
ylim([-1 1])
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))



% SAVE FIGURE
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'LZC_correlations'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'LZC_correlations'))

%% Regress out effects of age, and genotype on power

% Note: the regression MUST be done on the log-scaled values, otherwise we
% end up with negative power (nonsensical, no physical meaning) which we
% cannot be log-scaled later. Solution -- log-scale --> regress -->
% back-transform 

B_wk = nan(4,size(WP,1),size(WP,2));
R_wk = nan(size(WP));
B_sp = nan(4,size(SP,1),size(SP,2));
R_sp = nan(size(SP));

% del_wk_avg = log10(mean(delta_wk))';
% del_sp_avg = log10(mean(delta_sp))';

for i = 1:size(WP,1)
    for j = 1:size(WP,2)
        [b,~,r] = regress(squeeze(log10(WP(i,j,:))), [log2(all_ages)' allsite' allgenotype' ones(size(delta_wk,2),1)]);
        B_wk(:,i,j) = b;
        R_wk(i,j,:) = r;
    end
end

for i = 1:size(SP,1)
    for j = 1:size(SP,2)
        [b,~,r] = regress(squeeze(log10(SP(i,j,:))), [log2(all_ages)' allsite' allgenotype' ones(size(delta_sp,2),1)]);
        B_sp(:,i,j) = b;
        R_sp(i,j,:) = r;
    end
end

% Use residuals and back-transform to examine power with delta power effects regressed out

reg_WP = 10.^(R_wk+repmat(squeeze(B_wk(4,:,:)),1,1,size(R_wk,3))); % regress out effect of delta power (residuals + y-intercept)
reg_SP = 10.^(R_sp+repmat(squeeze(B_sp(4,:,:)),1,1,size(R_sp,3))); % regress out effect of delta power (residuals + y-intercept)

keyboard

% 
% 
% %% Also regress out Bayley scores for subjects who have them
% 
% % Note: the regression MUST be done on the log-scaled values, otherwise we
% % end up with negative power (nonsensical, no physical meaning) which we
% % cannot be log-scaled later. Solution -- log-scale --> regress -->
% % back-transform 
% 
% NDX = ~isnan(lng_score) & ~isnan(cog_score);
% 
% B_wk = nan(3,size(reg_WP,1),size(reg_WP,2));
% R_wk = nan(size(reg_WP(:,:,NDX)));
% B_sp = nan(3,size(reg_SP,1),size(reg_SP,2));
% R_sp = nan(size(reg_SP(:,:,NDX)));
% 
% for i = 1:size(reg_WP,1)
%     for j = 1:size(reg_WP,2)
%         [b,~,r] = regress(squeeze(log10(reg_WP(i,j,NDX))),[cog_score(NDX)' lng_score(NDX)' ones(length(cog_score(NDX)),1)]);
%         B_wk(:,i,j) = b;
%         R_wk(i,j,:) = r;
%     end
% end
% 
% for i = 1:size(reg_SP,1)
%     for j = 1:size(reg_SP,2)
%         [b,~,r] = regress(squeeze(log10(reg_SP(i,j,NDX))),[cog_score(NDX)' lng_score(NDX)' ones(length(cog_score(NDX)),1)]);
%         B_sp(:,i,j) = b;
%         R_sp(i,j,:) = r;
%     end
% end
% 
% reg_WP(:,:,NDX) = 10.^(R_wk+repmat(squeeze(B_wk(3,:,:)),1,1,size(R_wk,3))); % regress out effect of delta power (residuals + y-intercept)
% reg_SP(:,:,NDX) = 10.^(R_sp+repmat(squeeze(B_sp(3,:,:)),1,1,size(R_sp,3))); % regress out effect of delta power (residuals + y-intercept)
% 
% %% Also regress out STN scores for subjects who have them
% 
% % Note: the regression MUST be done on the log-scaled values, otherwise we
% % end up with negative power (nonsensical, no physical meaning) which we
% % cannot be log-scaled later. Solution -- log-scale --> regress -->
% % back-transform 
% 
% NDX = ~isnan(STN_TF);
% 
% B_wk = nan(2,size(reg_WP,1),size(reg_WP,2));
% R_wk = nan(size(reg_WP(:,:,NDX)));
% B_sp = nan(2,size(reg_SP,1),size(reg_SP,2));
% R_sp = nan(size(reg_SP(:,:,NDX)));
% 
% for i = 1:size(reg_WP,1)
%     for j = 1:size(reg_WP,2)
%         [b,~,r] = regress(squeeze(log10(reg_WP(i,j,NDX))),[STN_TF(NDX)' ones(length(STN_TF(NDX)),1)]);
%         B_wk(:,i,j) = b;
%         R_wk(i,j,:) = r;
%     end
% end
% 
% for i = 1:size(reg_SP,1)
%     for j = 1:size(reg_SP,2)
%         [b,~,r] = regress(squeeze(log10(reg_SP(i,j,NDX))),[STN_TF(NDX)' ones(length(STN_TF(NDX)),1)]);
%         B_sp(:,i,j) = b;
%         R_sp(i,j,:) = r;
%     end
% end
% 
% reg_WP(:,:,NDX) = 10.^(R_wk+repmat(squeeze(B_wk(2,:,:)),1,1,size(R_wk,3))); % regress out effect of delta power (residuals + y-intercept)
% reg_SP(:,:,NDX) = 10.^(R_sp+repmat(squeeze(B_sp(2,:,:)),1,1,size(R_sp,3))); % regress out effect of delta power (residuals + y-intercept)
% 
% %% PSD plots absolute power after regression
% 
% myfigure
% h=subplot(1,1,1);
% [~,~,wkci] = ttest( log10(squeeze(nanmean(reg_WP)))' );
% [~,~,spci] = ttest( log10(squeeze(nanmean(reg_SP)))' );
% foi_hd = 2.^[0:0.01:5]; % interpolation resolution for plotting
% wk = log10(squeeze(nanmean(reg_WP)))';
% sp = log10(squeeze(nanmean(reg_SP)))';
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
%     'AutoUpdate','off','FontSize',18)
% xlim([0 5])
% xlabel('Frequency (Hz)','FontSize',12)
% ylabel('Power log_1_0(\muV^{2}/Hz)','FontSize',18)
% title('Absolute Power','FontSize',20)
% set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
% box off
% %axis square
% legend boxoff
% %ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
% %annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
%     %%'LineStyle','none','FontSize',18);
% set(gca,'linewidth',3)
% xAX = get(gca,'XAxis');
% set(xAX,'FontSize', 22)
% yAX = get(gca,'YAxis');
% set(yAX,'FontSize', 22)
% set(gca, 'TickDir', 'out')
% figname = sprintf('PSDs_abs_power_regression_sp_wk');
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
% 
% clear wkci spci wk sp
% 
% %% PSD difference plot after regression
% 
% myfigure
% h=subplot(1,1,1);
% foi_hd = 2.^[0:0.01:5]; % interpolation resolution for plotting
% wk = log10(squeeze(nanmean(reg_WP)))';
% sp = log10(squeeze(nanmean(reg_SP)))';
% wsd_dB = mean(wk)-mean(sp);
% wsd = (10.^(wsd_dB/10).*100)-100; % convert dB to percent change
% [~,~,ci] = ttest(wk,sp);
% ci = (10.^(ci/10).*100)-100;% convert CI to percentage 
% lbl1 = 'Wake - Sleep';
% plot(log2(foi_hd),interp1(foi,wsd',foi_hd,'spline'),'m','LineWidth',4)
% hold on
% patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ci(1,:),foi_hd,'spline') ...
%     fliplr(interp1(foi,ci(2,:),foi_hd,'spline'))],'m','facealpha',0.3,'EdgeColor','none')
% legend({lbl1,sprintf('%s%s',lbl1,' 95% CI')},'AutoUpdate','off','FontSize',18)
% plot(log2(foi_hd),zeros(1,length(foi_hd)),'k--')
% xlim([0 5])
% xlabel('Frequency (Hz)','FontSize',12)
% %ylabel('Power log_1_0(\muV^{2}/Hz)','FontSize',18)
% ylabel('Power change from sleep (%)'); ylim([-25 25])
% title('Absolute Power','FontSize',20)
% set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
% box off
% %axis square
% legend boxoff
% %ano_str = sprintf('%s%i\n',' n = ',sum(~isnan(sum(wk+sp,2))));
% %annotation('textbox',[.65 0 0.3 .7],'String',ano_str,'FitBoxToText','on',...
%     %%'LineStyle','none','FontSize',18);
% set(gca,'linewidth',3)
% xAX = get(gca,'XAxis');
% set(xAX,'FontSize', 22)
% yAX = get(gca,'YAxis');
% set(yAX,'FontSize', 22)
% set(gca, 'TickDir', 'out')
% figname = sprintf('PSDs_after_reg_abs_power_DIFFERNCE');
% makefigpretty
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
% 
% %% power permutation clustering after regression
% 
% clst_pow_reg = RMpermclusttest(log10(reg_WP),log10(reg_SP),Nperm,p_thresh);
% 
% sgcl_pos = [];
% sgcl_neg = [];
% 
% if sum(clst_pow_reg.n_pos > 0) > 0
%     sgcl_pos = clst_pow_reg.P_val_pos < 0.05;
% end
% 
% if sum(clst_pow_reg.n_neg > 0) > 0
%     sgcl_neg = clst_pow_reg.P_val_neg < 0.05;
% end
% 
% cnt = 0;
% for icls = 1:length(sgcl_pos)
%     if sgcl_pos(icls)
%         cnt = cnt + 1;
%         myfigure
%         h=subplot(1,1,1);
%         plot(log2(foi),sum(clst_pow_reg.lbl_pos==icls,2),'linewidth',2)
%         ylim([0 20])
%         xlim([0 5])
%         set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
%         xlim([0 5])
%         xlabel('Frequency (Hz)')
%         ylabel('Number of channels')
%         makefigpretty
%         title(sprintf('power after reg postivie cluster, p = %1.3f',clst_pow_reg.P_val_pos(icls)),'fontsize',18)
%         figname = sprintf('Power after reg pos clust spect #%i',cnt);
%         print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
%         print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
%         
%         myfigure
%         plot_topo_AS(sum(clst_pow_reg.lbl_pos==icls)')
%         title(sprintf('power after reg pos clust, p = %1.3f',clst_pow_reg.P_val_pos(icls)),'fontsize',18)
%         colormap jet
%         caxis([0 size(reg_WP,2)])
%         mycolorbar
%         figname = sprintf('Power after reg pos clust toop #%i',cnt);
%         print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
%         print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
%     end
% end
% 
% cnt = 0;
% for icls = 1:length(sgcl_neg)
%     if sgcl_neg(icls)
%         cnt = cnt + 1;
%         myfigure
%         h=subplot(1,1,1);
%         plot(log2(foi),sum(clst_pow_reg.lbl_neg==icls,2),'linewidth',2)
%         ylim([0 20])
%         xlim([0 5])
%         set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
%         xlim([0 5])
%         xlabel('Frequency (Hz)')
%         ylabel('Number of channels')
%         makefigpretty
%         title(sprintf('power after reg negative cluster, p = %1.3f',clst_pow_reg.P_val_neg(icls)),'fontsize',18)
%         figname = sprintf('Power after reg neg clust spectrum #%i',cnt);
%         print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
%         print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
%         
%         myfigure
%         plot_topo_AS(sum(clst_pow_reg.lbl_neg==icls)')
%         title(sprintf('power after reg neg clust, p = %1.3f',clst_pow_reg.P_val_neg(icls)),'fontsize',18)
%         colormap jet
%         caxis([0 size(reg_WP,2)])
%         mycolorbar
%         figname = sprintf('Power after reg neg clust toop #%i',cnt);
%         print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
%         print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
%     end
% end
% 

% %% Individual subject changes (more figures for reviewers)
% 
% % change the color order default
% 
% chrome = hsv;
% chrome = chrome(1:1.5:end,:)./1.333;
% set(groot,'defaultAxesColorOrder',chrome)
% scat_chrome = [chrome(1:length(W_POW),:); chrome(1:length(W_POW),:)];
% 
% % Power
% 
% W_POW = nan(1,size(WP,3));
% S_POW = nan(1,size(SP,3));
% 
% W_tmp = nan(size(WP,1),size(WP,3));
% S_tmp = nan(size(SP,1),size(SP,3));
% for ich = 1:size(WP,1)
%     for isb = 1:size(WP,3)
%         W_tmp(ich,isb) = trapz(foi,WP(ich,:,isb));
%         S_tmp(ich,isb) = trapz(foi,SP(ich,:,isb));
%     end
% end
% 
% 
% W_POW = log10(squeeze(mean(W_tmp)))'; % average across channels and frequency
% S_POW = log10(squeeze(mean(S_tmp)))';
% 
% myfigure
% h=subplot(1,1,1);
% plot([W_POW S_POW]')
% scatter([ones(length(W_POW),1); ones(length(W_POW),1).*2],[W_POW; S_POW],67,scat_chrome,'filled')
% makefigpretty
% xlim([0.8 2.2])
% xticks([1 2])
% xticklabels({'Awake','Asleep'})
% ax1 = gca;                   % gca = get current axis
% ax1.XAxis.Visible = 'off';   % remove x-axis
% title('Power','fontsize',25')
% print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'Power_dot_plot'))
% 
% 
% % MSE
% 
% W_MSE = squeeze(mean(mean(wk_all))); % average across channels and time-scales
% S_MSE = squeeze(mean(mean(sp_all)));
% 
% myfigure
% h=subplot(1,1,1);
% plot([W_MSE S_MSE]')
% scatter([ones(length(W_MSE),1); ones(length(W_MSE),1).*2],[W_MSE; S_MSE],67,scat_chrome,'filled')
% makefigpretty
% xlim([0.8 2.2])
% xticks([1 2])
% xticklabels({'Awake','Asleep'})
% ax1 = gca;                   % gca = get current axis
% ax1.XAxis.Visible = 'off';   % remove x-axis
% title('MSE','fontsize',25')
% print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'MSE_dot_plot'))
% 
% % LZC
% 
% W_LZC = squeeze(mean(mean(wk_lzc_all))); % average across channels and time-scales
% S_LZC = squeeze(mean(mean(sp_lzc_all)));
% 
% myfigure
% h=subplot(1,1,1);
% plot([W_LZC S_LZC]')
% scatter([ones(length(W_LZC),1); ones(length(W_LZC),1).*2],[W_LZC; S_LZC],67,scat_chrome,'filled')
% makefigpretty
% xlim([0.8 2.2])
% xticks([1 2])
% xticklabels({'Awake','Asleep'})
% ax1 = gca;                   % gca = get current axis
% ax1.XAxis.Visible = 'off';   % remove x-axis
% title('LZC','fontsize',25')
% print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'LZC_dot_plot'))

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
    'FontSize',18,'location','southeast')
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

myfigure
plot_topo_AS(cd_pow_top')
title('Spectral power effect sizes \delta frequencies','fontsize',18)
mycolorbar
caxis([-1 1])
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'CohenD_pow_topo'))

myfigure
plot_topo_AS(cd_lzc_top')
title('gMLZ effect sizes \delta scales','fontsize',18)
mycolorbar
caxis([-1 1])
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'CohenD_lzc_topo'))


%% MSE Correlation

% allocation
r_pow_MSE_wk = nan(2,2,20);
r_pow_MSE_sp = nan(2,2,20);
r_pow_MSE_wk_lo = nan(2,2,20);
r_pow_MSE_sp_lo = nan(2,2,20);
r_pow_MSE_wk_hi = nan(2,2,20);
r_pow_MSE_sp_hi = nan(2,2,20);
p_pow_MSE_wk = nan(2,2,20);
p_pow_MSE_sp = nan(2,2,20);

for i = 1:size(wk_scl,1)
    [r_pow_MSE_wk(:,:,i),p_pow_MSE_wk(:,:,i),r_pow_MSE_wk_lo(:,:,i),r_pow_MSE_wk_hi(:,:,i)] ...
        = corrcoef(log10(mean(delta_wk)),wk_scl(i,:));
    [r_pow_MSE_sp(:,:,i),p_pow_MSE_sp(:,:,i),r_pow_MSE_sp_lo(:,:,i),r_pow_MSE_sp_hi(:,:,i)] ...
        = corrcoef(log10(mean(delta_sp)),sp_scl(i,:));
end

% discard extraneous elements, square, and multiply by sign

r_pow_MSE_wk    = squeeze( r_pow_MSE_wk(1,2,:).^2.*sign( r_pow_MSE_wk(1,2,:)) )';
r_pow_MSE_sp    = squeeze(r_pow_MSE_sp(1,2,:).^2.*sign( r_pow_MSE_sp(1,2,:)))';
r_pow_MSE_wk_lo = squeeze(r_pow_MSE_wk_lo(1,2,:).^2.*sign( r_pow_MSE_wk_lo(1,2,:)))';
r_pow_MSE_sp_lo = squeeze(r_pow_MSE_sp_lo(1,2,:).^2.*sign( r_pow_MSE_sp_lo(1,2,:)))';
r_pow_MSE_wk_hi = squeeze(r_pow_MSE_wk_hi(1,2,:).^2.*sign( r_pow_MSE_wk_hi(1,2,:)))';
r_pow_MSE_sp_hi = squeeze(r_pow_MSE_sp_hi(1,2,:).^2.*sign( r_pow_MSE_sp_hi(1,2,:)))';
p_pow_MSE_wk    = squeeze(p_pow_MSE_wk(1,2,:))';
p_pow_MSE_sp    = squeeze(p_pow_MSE_sp(1,2,:))';

% R^2 x sign(r)

myfigure
h=subplot(1,1,1)
hold on
plot(log2(nyqfoi),r_pow_MSE_wk,'r','linewidth',2)
patch([log2(nyqfoi(scales)) fliplr(log2(nyqfoi(scales)))] , [r_pow_MSE_wk_lo fliplr(r_pow_MSE_wk_hi)], 'r','facealpha', ...
    1.0,'EdgeColor','none')

plot(log2(nyqfoi),r_pow_MSE_sp,'b','linewidth',2)
patch([log2(nyqfoi(scales)) fliplr(log2(nyqfoi(scales)))] , [r_pow_MSE_sp_lo fliplr(r_pow_MSE_sp_hi)], 'b','facealpha', ...
    1.0,'EdgeColor','none')
legend({'Awake R^{2} x sign(r)','Awake 95% CI','Asleep R^{2} x sign(r)','Asleep 95% CI'} ,...
    'location','northeast','fontsize',18,'autoupdate','off')
legend boxoff
plot(zeros(1,length(nyqfoi)),'-.','color','k','linewidth',2)
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
ylabel('R^{2} x sign(r)')
makefigpretty
ylim([-1 1])
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
% SAVE FIGURE
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'MSE_correlations'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'MSE_correlations'))


%% gMLZ Correlation

% allocation
r_pow_gMLZ_wk = nan(2,2,20);
r_pow_gMLZ_sp = nan(2,2,20);
r_pow_gMLZ_wk_lo = nan(2,2,20);
r_pow_gMLZ_sp_lo = nan(2,2,20);
r_pow_gMLZ_wk_hi = nan(2,2,20);
r_pow_gMLZ_sp_hi = nan(2,2,20);
p_pow_gMLZ_wk = nan(2,2,20);
p_pow_gMLZ_sp = nan(2,2,20);

for i = 1:size(wk_scl_lzc,1)
    [r_pow_gMLZ_wk(:,:,i),p_pow_gMLZ_wk(:,:,i),r_pow_gMLZ_wk_lo(:,:,i),r_pow_gMLZ_wk_hi(:,:,i)] ...
        = corrcoef(log10(mean(delta_wk)),wk_scl_lzc(i,:));
    [r_pow_gMLZ_sp(:,:,i),p_pow_gMLZ_sp(:,:,i),r_pow_gMLZ_sp_lo(:,:,i),r_pow_gMLZ_sp_hi(:,:,i)] ...
        = corrcoef(log10(mean(delta_sp)),sp_scl_lzc(i,:));
end

% discard extraneous elements

r_pow_gMLZ_wk    = squeeze( r_pow_gMLZ_wk(1,2,:).^2.*sign( r_pow_gMLZ_wk(1,2,:)) )';
r_pow_gMLZ_sp    = squeeze(r_pow_gMLZ_sp(1,2,:).^2.*sign( r_pow_gMLZ_sp(1,2,:)))';
r_pow_gMLZ_wk_lo = squeeze(r_pow_gMLZ_wk_lo(1,2,:).^2.*sign( r_pow_gMLZ_wk_lo(1,2,:)))';
r_pow_gMLZ_sp_lo = squeeze(r_pow_gMLZ_sp_lo(1,2,:).^2.*sign( r_pow_gMLZ_sp_lo(1,2,:)))';
r_pow_gMLZ_wk_hi = squeeze(r_pow_gMLZ_wk_hi(1,2,:).^2.*sign( r_pow_gMLZ_wk_hi(1,2,:)))';
r_pow_gMLZ_sp_hi = squeeze(r_pow_gMLZ_sp_hi(1,2,:).^2.*sign( r_pow_gMLZ_sp_hi(1,2,:)))';
p_pow_gMLZ_wk    = squeeze(p_pow_gMLZ_wk(1,2,:))';
p_pow_gMLZ_sp    = squeeze(p_pow_gMLZ_sp(1,2,:))';

% R^2 x sign(r)

myfigure
h=subplot(1,1,1)
hold on
plot(log2(LZC_foi),r_pow_gMLZ_wk,'r','linewidth',2)
patch([log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [r_pow_gMLZ_wk_lo fliplr(r_pow_gMLZ_wk_hi)], 'r','facealpha', ...
    1.0,'EdgeColor','none')

plot(log2(LZC_foi),r_pow_gMLZ_sp,'b','linewidth',2)
patch([log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [r_pow_gMLZ_sp_lo fliplr(r_pow_gMLZ_sp_hi)], 'b','facealpha', ...
    1.0,'EdgeColor','none')
legend({'Awake R^{2} x sign(r)','Awake 95% CI','Asleep R^{2} x sign(r)','Asleep 95% CI'} ,...
    'location','northeast','fontsize',18,'autoupdate','off')
legend boxoff
plot(log2(LZC_foi),zeros(1,length(LZC_foi)),'-.','color','k','linewidth',2)

ylabel('R^{2} x sign(r)')
xlabel('Center frequency (Hz)')
ylabel('R^{2} x sign(r)')
XTICK = [0:1:5];
set(h,'XTick',XTICK,'XTickLabel',2.^XTICK)
xlim([0 5])
ylim([-1 1])
makefigpretty
% SAVE FIGURE
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_correlations'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'gMLZ_correlations'))

%% Regress out effects of delta power on MSE

B_wk = nan(2,size(wk_all,1),size(wk_all,2));
R_wk = nan(size(wk_all));
B_sp = nan(2,size(sp_all,1),size(sp_all,2));
R_sp = nan(size(sp_all));

for i = 1:size(wk_all,1)
    for j = 1:size(wk_all,2)
        [b,~,r] = regress(squeeze(wk_all(i,j,:)),[log10(delta_wk(i,:))' ones(size(delta_wk,2),1)]);
        B_wk(:,i,j) = b;
        R_wk(i,j,:) = r;
    end
end

for i = 1:size(sp_all,1)
    for j = 1:size(sp_all,2)
        [b,~,r] = regress(squeeze(sp_all(i,j,:)),[log10(delta_sp(i,:))' ones(size(delta_sp,2),1)]);
        B_sp(:,i,j) = b;
        R_sp(i,j,:) = r;
    end
end

% Use residuals to examine MSE with delta power effects regressed out

regout_wk = R_wk+repmat(squeeze(B_wk(2,:,:)),1,1,size(R_wk,3)); % regress out effect of delta power (residuals + y-intercept)
regout_sp = R_sp+repmat(squeeze(B_sp(2,:,:)),1,1,size(R_sp,3)); % regress out effect of delta power (residuals + y-intercept)

% %% Also regress out Bayley score effects on MSE for subjects who have them
% 
% NDX = ~isnan(lng_score) & ~isnan(cog_score);
% 
% B_wk = nan(3,size(regout_wk,1),size(regout_wk,2));
% R_wk = nan(size(regout_wk(:,:,NDX)));
% B_sp = nan(3,size(regout_sp,1),size(regout_sp,2));
% R_sp = nan(size(regout_sp(:,:,NDX)));
% 
% for i = 1:size(regout_wk,1)
%     for j = 1:size(regout_wk,2)
%         [b,~,r] = regress(squeeze(regout_wk(i,j,NDX)),[cog_score(NDX)' lng_score(NDX)' ones(length(cog_score(NDX)),1)]);
%         B_wk(:,i,j) = b;
%         R_wk(i,j,:) = r;
%     end
% end
% 
% for i = 1:size(regout_sp,1)
%     for j = 1:size(regout_sp,2)
%         [b,~,r] = regress(squeeze(regout_sp(i,j,NDX)),[cog_score(NDX)' lng_score(NDX)' ones(length(cog_score(NDX)),1)]);
%         B_sp(:,i,j) = b;
%         R_sp(i,j,:) = r;
%     end
% end
% 
% regout_wk(:,:,NDX) = R_wk+repmat(squeeze(B_wk(3,:,:)),1,1,size(R_wk,3)); % regress out effect of delta power (residuals + y-intercept)
% regout_sp(:,:,NDX) = R_sp+repmat(squeeze(B_sp(3,:,:)),1,1,size(R_sp,3)); % regress out effect of delta power (residuals + y-intercept)
% 
% %% Also regress out STN scores from MSE for subjects who have them
% 
% 
% NDX = ~isnan(STN_TF);
% 
% B_wk = nan(2,size(regout_wk,1),size(regout_wk,2));
% R_wk = nan(size(regout_wk(:,:,NDX)));
% B_sp = nan(2,size(regout_sp,1),size(regout_sp,2));
% R_sp = nan(size(regout_sp(:,:,NDX)));
% 
% for i = 1:size(regout_wk,1)
%     for j = 1:size(regout_wk,2)
%         [b,~,r] = regress(squeeze(log10(regout_wk(i,j,NDX))),[STN_TF(NDX)' ones(length(STN_TF(NDX)),1)]);
%         B_wk(:,i,j) = b;
%         R_wk(i,j,:) = r;
%     end
% end
% 
% for i = 1:size(regout_sp,1)
%     for j = 1:size(regout_sp,2)
%         [b,~,r] = regress(squeeze(log10(regout_sp(i,j,NDX))),[STN_TF(NDX)' ones(length(STN_TF(NDX)),1)]);
%         B_sp(:,i,j) = b;
%         R_sp(i,j,:) = r;
%     end
% end
% 
% regout_wk(:,:,NDX) = 10.^(R_wk+repmat(squeeze(B_wk(2,:,:)),1,1,size(R_wk,3))); % regress out effect of delta power (residuals + y-intercept)
% regout_sp(:,:,NDX) = 10.^(R_sp+repmat(squeeze(B_sp(2,:,:)),1,1,size(R_sp,3))); % regress out effect of delta power (residuals + y-intercept)

% %% Mediation analysis (MSE)
% 
% sleep = cat(1,zeros(size(wk_scl,2),1),ones(size(wk_scl,2),1));
% options.alpha = 0.05;
% options.verbose = true;
% options.display = true;
% medP = nan(1,size(wk_scl,1));
% medCI = nan(2,size(wk_scl,1));
% 
% for iscl = 1:size(wk_scl,1)
%     out = mediationAnalysis0(sleep,[log10(mean(delta_wk)) ...
%         log10(mean(delta_sp))]',[wk_scl(iscl,:) sp_scl(iscl,:)]',options);
%     medP(iscl) = out.montecarlo.p;
%     medCI(:,iscl) = out.montecarlo.IC;
% end
% 
% medQ = mafdr(medP,'BHFDR',true);
% 
% for iscl = 1:size(wk_scl,1)
%     if medQ(iscl) < 0.05
%         fprintf('\n%s for MSE @ scale #%i, p = %1.2f FDR corrected\n','Full mediation',iscl,medQ(iscl))
%     else
%         fprintf('\n%s for MSE @ scale #%i, p = %1.2f FDR corrected\n','No mediation',iscl,medQ(iscl))
%     end
% end
% 
% myfigure
% hold on
% plot(log2(nyqfoi),-log10(medQ),'linewidth',2)
% plot(log2(nyqfoi),ones(1,length(medQ)).*-log10(0.05),'k--')
% ylim([0 4])
% xlim([log2(5) log2(100)])
% legend({'corrected p-value (nonparametric test)','alpha level'},'location','best','fontsize',14)
% legend boxoff
% ticks = 2:1:7;
% xticks(ticks)
% ticklabels = cell(1,length(ticks));
% for i = 1:length(ticklabels)
%     ticklabels{i} = num2str(2^ticks(i));
% end
% xticklabels(ticklabels)
% xlabel('Nyquist frequency (Hz)')
% ylabel('-log_{10}(p)')
% makefigpretty
% ax1 = gca; % current axes
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation',...
%     'right', 'Color','none');
% makefigpretty
% line(log2(nyqfoi),wsdelta,'Parent',ax2,'Color','none')
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
% title(sprintf('MSE mediation analysis p-values'),'fontsize',18)
% figname = sprintf('MSE mediation analysis p-values');
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
% 
% %% MSE mediation effect size (nonparametric)
% 
% myfigure
% hold on
% plot(log2(LZC_foi),mean(medCI),'k','linewidth',2)
% patch( [log2(LZC_foi) fliplr(log2(LZC_foi))] , [medCI(1,:) fliplr(medCI(2,:))], 'k','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('Effect size','Effect size 95% CI','AutoUpdate','off',...
%     'FontSize',18,'location','northeast')
% legend boxoff
% xlim([1 length(scales)])
% title('gMLZ values covary nuisance over all time scales ','fontsize',18)
% xlabel('Center frequency (Hz)')
% ylabel('LZ complexity'); ylim([0 0.8])
% set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
% xlim([0 5])
% makefigpretty
% % SAVE FIGURE
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'MSE_mediation_ES'))
% print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'MSE_mediation_ES'))

%% plot average MSE scales after regression

myfigure
hold on
[~,~,ci] = ttest(squeeze(mean(regout_wk))');
plot(log2(nyqfoi),mean(squeeze(mean(regout_wk)),2),'r','linewidth',4)
patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
    0.3,'EdgeColor','none')
[~,~,ci] = ttest(squeeze(mean(regout_sp))');
plot(log2(nyqfoi),mean(squeeze(mean(regout_sp)),2),'b','linewidth',4)
patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [ci(1,:) fliplr(ci(2,:))], 'b','facealpha', ...
    0.3,'EdgeColor','none')
legend('wake','wake 95% CI','sleep','95% sleep','AutoUpdate','off',...
    'FontSize',18,'location','northeast')
legend boxoff
%ylim([0 2])
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
ylabel('mSampEn'); ylim([0 2])
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
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'MSE_scales_with_reg'))

%% plot average MSE scales after regression compared to TD awake

myfigure
hold on
[~,~,ci] = ttest(squeeze(mean(regout_wk))');
plot(log2(nyqfoi),mean(squeeze(mean(regout_wk)),2),'r','linewidth',2)
patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
    0.3,'EdgeColor','none')

[~,~,ci] = ttest(td_scl');
plot(log2(nyqfoi),mean(td_scl,2),'color',[0 0.5 0.5],'linewidth',2)
patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [ci(1,:) fliplr(ci(2,:))],[0 0.5 0.5],'facealpha', ...
    0.3,'EdgeColor','none')
legend('AS wake (covary delta)','AS wake (covary delta) 95% CI','TD wake','TD wake 95% CI','AutoUpdate','off',...
    'FontSize',18,'location','southwest')
legend boxoff
%ylim([0 2])
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
ylabel('mSampEn'); ylim([0 2])
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
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'MSE_scales_with_reg_TD'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'MSE_scales_with_reg_TD'))

%% MSE percent change from sleep after regressing out covariates

myfigure
hold on
wsdelta = (mean(squeeze(mean(regout_wk)),2) - mean(squeeze(mean(regout_sp)),2)) ...
    ./mean(squeeze(mean(regout_sp)),2).*100;
[~,~,ci] = ttest(squeeze(mean(regout_wk)),squeeze(mean(regout_sp)),'dim',2); ci = ci';
plot(log2(nyqfoi),wsdelta ,'m','linewidth',4)
patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [ci(1,:)./mean(squeeze(mean(regout_sp)),2)'.*100 ...
    fliplr(ci(2,:)./mean(squeeze(mean(regout_sp)),2)'.*100)], 'm','facealpha', 1.0,'EdgeColor','none')
legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',18,'location','southeast')
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
ylabel('mSampEn change from sleep (%)'); ylim([0 40]);
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
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'MSE_scales_delta_with_reg'))

%% MSE percent change from sleep after regressing out covariates -- RAW VERSION

myfigure
hold on
wsdelta = squeeze(mean(mean(regout_wk,3),1) - mean(mean(regout_sp,3),1));
%wsdelta = (mean(squeeze(mean(regout_wk)),2) - mean(squeeze(mean(regout_sp))));
[~,~,ci] = ttest(squeeze(mean(regout_wk)),squeeze(mean(regout_sp)),'dim',2); ci = ci';
plot(log2(nyqfoi),wsdelta ,'m','linewidth',1)
patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [ci(1,:) fliplr(ci(2,:))], 'm','facealpha', 0.3,'EdgeColor','none')
legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',18,'location','southeast')
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
ylabel('mSampEn change from sleep (%)'); 
ylim([-0.25 0.25]);
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
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'MSE_scales_delta_with_reg_RAW'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'MSE_scales_delta_with_reg_RAW'))


%% MSE permutation clustering after regression

clst_rgt = RMpermclusttest(regout_wk,regout_sp,Nperm,p_thresh);

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
        myfigure
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
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
        
        myfigure
        plot_topo_AS(sum(clst_rgt.lbl_pos==icls)')
        title(sprintf('MSE with reg pos clust, p = %1.3f',clst_rgt.P_val_pos(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(wk_lzc_all,2)])
        mycolorbar
        figname = sprintf('MSE with reg pos clust toop #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
    end
    
    % dot plots for all clusters (even if p > 0.05)
    if any(~isnan(clst_rgt.lbl_pos(:)))            
        dot_slope_plot(regout_wk,regout_sp,clst_rgt.lbl_pos,icls,all_ages,allgenotype,p_thresh)
        title(sprintf('MSE, d_median = %1.2f',clst_rgt.d_pos(icls)))
        ylabel('SampEn')
        %ylim([0.7 1.7])
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s%i.svg',DIRFIGURE,'MSE_regout_pos_clust_dot_plot_#',icls))
        print(gcf,'-dpng','-r300','-painters',sprintf('%s%s%i.png',DIRFIGURE,'MSE_regout_pos_clust_dot_plot_#',icls))
    end
end
cnt = 0;
for icls = 1:length(sgcl_neg)
    if sgcl_neg(icls)
        cnt = cnt + 1;
        myfigure
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
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
        
        myfigure
        plot_topo_AS(sum(clst_rgt.lbl_neg==icls)')
        title(sprintf('MSE neg clust, p = %1.3f',clst_rgt.P_val_neg(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(wk_lzc_all,2)])
        mycolorbar
        figname = sprintf('MSE with reg negs clust topo #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
    end
    
    % dot plots for all clusters (even if p > 0.05)
    if any(~isnan(clst_rgt.lbl_neg(:)))            
        dot_slope_plot(regout_wk,regout_sp,clst_rgt.lbl_neg,icls,all_ages,allgenotype,p_thresh)
        title(sprintf('MSE, d_median = %1.2f',clst_rgt.d_neg(icls)))
        ylabel('SampEn')
        %ylim([0.7 1.7])
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s%i.svg',DIRFIGURE,'MSE_regout_neg_clust_dot_plot_#',icls))
        print(gcf,'-dpng','-r300','-painters',sprintf('%s%s%i.png',DIRFIGURE,'MSE_regout_neg_clust_dot_plot_#',icls))
    end
end

%% Regress out effects of delta power on gMLZ

B_wk = nan(2,size(wk_lzc_all,1),size(wk_lzc_all,2));
R_wk = nan(size(wk_lzc_all));
B_sp = nan(2,size(sp_lzc_all,1),size(sp_lzc_all,2));
R_sp = nan(size(sp_lzc_all));

% del_wk_avg = log10(mean(delta_wk))';
% del_sp_avg = log10(mean(delta_sp))';

for i = 1:size(wk_lzc_all,1)
    for j = 1:size(wk_lzc_all,2)
        [b,~,r] = regress(squeeze(wk_lzc_all(i,j,:)),[log10(delta_wk(i,:))' ones(size(delta_wk,2),1)]);
        B_wk(:,i,j) = b;
        R_wk(i,j,:) = r;
    end
end

for i = 1:size(sp_lzc_all,1)
    for j = 1:size(sp_lzc_all,2)
        [b,~,r] = regress(squeeze(sp_lzc_all(i,j,:)),[log10(delta_wk(i,:))' ones(size(delta_wk,2),1)]);
        B_sp(:,i,j) = b;
        R_sp(i,j,:) = r;
    end
end

% Use residuals to examine gMLZ with delta power effects regressed out

reg_lzc_wk = R_wk+repmat(squeeze(B_wk(2,:,:)),1,1,size(R_wk,3)); % regress out effect of delta power (residuals + y-intercept)
reg_lzc_sp = R_sp+repmat(squeeze(B_sp(2,:,:)),1,1,size(R_sp,3)); % regress out effect of delta power (residuals + y-intercept)

return
% %% Also regress out Bayley scores for subjects who have them
% 
% NDX = ~isnan(lng_score) & ~isnan(cog_score);
% 
% B_wk = nan(3,size(reg_lzc_wk,1),size(reg_lzc_wk,2));
% R_wk = nan(size(reg_lzc_wk(:,:,NDX)));
% B_sp = nan(3,size(reg_lzc_sp,1),size(reg_lzc_sp,2));
% R_sp = nan(size(reg_lzc_sp(:,:,NDX)));
% 
% for i = 1:size(reg_lzc_wk,1)
%     for j = 1:size(reg_lzc_wk,2)
%         [b,~,r] = regress(squeeze(reg_lzc_wk(i,j,NDX)),[cog_score(NDX)' lng_score(NDX)' ones(length(cog_score(NDX)),1)]);
%         B_wk(:,i,j) = b;
%         R_wk(i,j,:) = r;
%     end
% end
% 
% for i = 1:size(reg_lzc_sp,1)
%     for j = 1:size(reg_lzc_sp,2)
%         [b,~,r] = regress(squeeze(reg_lzc_sp(i,j,NDX)),[cog_score(NDX)' lng_score(NDX)' ones(length(cog_score(NDX)),1)]);
%         B_sp(:,i,j) = b;
%         R_sp(i,j,:) = r;
%     end
% end
% 
% reg_lzc_wk(:,:,NDX) = R_wk+repmat(squeeze(B_wk(3,:,:)),1,1,size(R_wk,3)); % regress out effect of delta power (residuals + y-intercept)
% reg_lzc_sp(:,:,NDX) = R_sp+repmat(squeeze(B_sp(3,:,:)),1,1,size(R_sp,3)); % regress out effect of delta power (residuals + y-intercept)
% 
% %% Also regress out STN scores from gMLZ for subjects who have them
% 
% 
% NDX = ~isnan(STN_TF);
% 
% B_wk = nan(2,size(reg_lzc_wk,1),size(reg_lzc_wk,2));
% R_wk = nan(size(reg_lzc_wk(:,:,NDX)));
% B_sp = nan(2,size(reg_lzc_sp,1),size(reg_lzc_sp,2));
% R_sp = nan(size(reg_lzc_sp(:,:,NDX)));
% 
% for i = 1:size(reg_lzc_wk,1)
%     for j = 1:size(reg_lzc_wk,2)
%         [b,~,r] = regress(squeeze(log10(reg_lzc_wk(i,j,NDX))),[STN_TF(NDX)' ones(length(STN_TF(NDX)),1)]);
%         B_wk(:,i,j) = b;
%         R_wk(i,j,:) = r;
%     end
% end
% 
% for i = 1:size(reg_lzc_sp,1)
%     for j = 1:size(reg_lzc_sp,2)
%         [b,~,r] = regress(squeeze(log10(reg_lzc_sp(i,j,NDX))),[STN_TF(NDX)' ones(length(STN_TF(NDX)),1)]);
%         B_sp(:,i,j) = b;
%         R_sp(i,j,:) = r;
%     end
% end
% 
% reg_lzc_wk(:,:,NDX) = 10.^(R_wk+repmat(squeeze(B_wk(2,:,:)),1,1,size(R_wk,3))); % regress out effect of delta power (residuals + y-intercept)
% reg_lzc_sp(:,:,NDX) = 10.^(R_sp+repmat(squeeze(B_sp(2,:,:)),1,1,size(R_sp,3))); % regress out effect of delta power (residuals + y-intercept)

%% Mediation analysis (gMLZ)

sleep = cat(1,zeros(size(wk_scl_lzc,2),1),ones(size(wk_scl_lzc,2),1));
options.alpha = 0.05;
options.verbose = false;
options.display = false;
medP = nan(1,size(wk_scl_lzc,1));
medCI = nan(2,size(wk_scl_lzc,1));

for iscl = 1:size(wk_scl_lzc,1)
    out = mediationAnalysis0(sleep,[log10(mean(delta_wk)) ...
        log10(mean(delta_sp))]',[wk_scl_lzc(iscl,:) sp_scl(iscl,:)]',options);
    medP(iscl) = out.montecarlo.p;
    medCI(:,iscl) = out.montecarlo.IC;
end

medQ = mafdr(medP,'BHFDR',true);

for iscl = 1:size(wk_scl_lzc,1)
    if medQ(iscl) < 0.05
        fprintf('\n%s for gMLZ @ scale #%i, p = %1.2f FDR corrected\n','Full mediation',iscl,medQ(iscl))
    else
        fprintf('\n%s for gMLZ @ scale #%i, p = %1.2f FDR corrected\n','No mediation',iscl,medQ(iscl))
    end
end

% %% gMLZ mediation analysis p-value
% 
% myfigure
% h=subplot(1,1,1)
% hold on
% plot(log2(LZC_foi),-log10(medQ),'linewidth',2)
% plot(log2(LZC_foi),ones(1,length(medQ)).*-log10(0.05),'k--')
% makefigpretty
% ylim([0 4])
% xlim([0 log2(LZC_foi(end))])
% legend({'corrected p-value (Sobel test)','alpha level'},'location','best','fontsize',14)
% legend boxoff
% xlabel('Center frequency (Hz)')
% ylabel('-log10(p)')
% set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
% title(sprintf('gMLZ mediation analysis p-values'),'fontsize',18)
% figname = sprintf('gMLZ mediation analysis p-values');
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
% print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
% 
% %% gMLZ mediation effect size nonparametric
% 
% myfigure
% hold on
% plot(log2(nyqfoi),mean(medCI),'k','linewidth',2)
% patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [medCI(1,:) fliplr(medCI(2,:))], 'k','facealpha', ...
%     1.0,'EdgeColor','none')
% legend('wake','wake 95% CI','sleep','95% sleep','AutoUpdate','off',...
%     'FontSize',18,'location','northeast')
% legend boxoff
% %ylim([-0.2 1])
% %title('gMLZ values over all time scales (no regression)','fontsize',18)
% xlim([log2(5) log2(100)])
% ticks = 2:1:7;
% xticks(ticks)
% ticklabels = cell(1,length(ticks));
% for i = 1:length(ticklabels)
%     ticklabels{i} = num2str(2^ticks(i));
% end
% xticklabels(ticklabels)
% xlabel('Nyquist frequency (Hz)')
% ylabel('sample entropy')
% makefigpretty
% % ADD TOP X-AXIS FOR TAU
% ax1 = gca; % current axes
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position',ax1_pos, 'XAxisLocation','top','YAxisLocation',...
%     'right', 'Color','none');
% makefigpretty
% line(log2(nyqfoi),mean(squeeze(mean(regout_wk)),2),'Parent',ax2,'Color','none')
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
% print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_mediation_ES'))
% print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'gMLZ_mediation_ES'))

%% plot average gMLZ scales after regression

myfigure
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
legend('wake','wake 95% CI','sleep','95% sleep','AutoUpdate','off','FontSize',18,'location','southeast')
legend boxoff
xlim([1 length(scales)])
title('gMLZ values covary nuisance over all time scales ','fontsize',18)
xlabel('Center frequency (Hz)')
ylabel('LZ complexity'); ylim([0 0.8])
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
xlim([0 5])
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_scales'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'gMLZ_scales'))

%% Compare awake AS gMLZ with regression to TD awake

myfigure
h=subplot(1,1,1);
hold on
[~,~,ci] = ttest(squeeze(mean(reg_lzc_wk))');
plot(log2(LZC_foi(scales)),mean(squeeze(mean(reg_lzc_wk)),2),'r','linewidth',2)
patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], 'r','facealpha', ...
    0.3,'EdgeColor','none')
[~,~,ci] = ttest(TD_lzc_scl');
plot(log2(LZC_foi(scales)),mean(TD_lzc_scl,2),'color',[0 0.5 0.5],'linewidth',2)
 patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) fliplr(ci(2,:))], [0 0.5 0.5],'facealpha', ...
     0.3,'EdgeColor','none')
legend('AS wake (covary delta)','AS wake (covary delta) 95% CI','TD wake','TD wake 95% CI','AutoUpdate','off','FontSize',18,'location','southeast')
legend boxoff
xlim([1 length(scales)])
title('gMLZ values covary nuisance over all time scales ','fontsize',18)
xlabel('Center frequency (Hz)')
ylabel('LZ complexity'); ylim([0 0.8])
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
xlim([0 5])
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_scales'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'AS_reg_TD_gMLZ'))

%% gMLZ percent change from sleep

myfigure
h=subplot(1,1,1);
hold on
wsdelta = (mean(squeeze(mean(reg_lzc_wk)),2) - mean(squeeze(mean(reg_lzc_sp)),2)) ...
    ./mean(squeeze(mean(reg_lzc_sp)),2).*100;
[~,~,ci] = ttest(squeeze(mean(reg_lzc_wk)),squeeze(mean(reg_lzc_sp)),'dim',2); ci = ci';
plot(log2(LZC_foi(scales)),wsdelta ,'m','linewidth',4)
patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:)./mean(squeeze(mean(reg_lzc_sp)),2)'.*100 ...
    fliplr(ci(2,:)./mean(squeeze(mean(reg_lzc_sp)),2)'.*100)], 'm','facealpha', 1.0,'EdgeColor','none')
legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',18,'location','southeast')
legend boxoff
plot(log2(LZC_foi(scales)),zeros(1,length(wsdelta)),'-.k')
xlim([0 5])
%ylim([-0.04 0.04])
title('gMLZ covary nuisance difference over all time scales ','fontsize',18)
xlabel('Center frequency (Hz)')
ylabel('gLZ change from sleep (%)'); ones(1,length(foi));
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_scales_delta'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'gMLZ_scales_delta'))

%% gMLZ percent change from sleep -- RAW VERSION

myfigure
h=subplot(1,1,1);
hold on

wsdelta = squeeze(mean(mean(reg_lzc_wk,1),3) -  mean(mean(reg_lzc_sp,1),3));

% (mean(squeeze(mean(reg_lzc_wk)),2) - mean(squeeze(mean(reg_lzc_sp)),2)) ...
%     ./mean(squeeze(mean(reg_lzc_sp)),2).*100;
% 
% wsdelta = (mean(squeeze(mean(reg_lzc_wk)),2) - mean(squeeze(mean(reg_lzc_sp)),2)) ...
%     ./mean(squeeze(mean(reg_lzc_sp)),2).*100;
[~,~,ci] = ttest(squeeze(mean(reg_lzc_wk)),squeeze(mean(reg_lzc_sp)),'dim',2); ci = ci';
plot(log2(LZC_foi(scales)),wsdelta ,'m','linewidth',4)
patch( [log2(LZC_foi(scales)) fliplr(log2(LZC_foi(scales)))] , [ci(1,:) ...
    fliplr(ci(2,:))], 'm','facealpha', 1.0,'EdgeColor','none')
legend('wake-sleep','wake-sleep 95% CI','AutoUpdate','off','FontSize',18,'location','southeast')
legend boxoff
plot(log2(LZC_foi(scales)),zeros(1,length(wsdelta)),'-.k')
xlim([0 5])
ylim([-0.2 0.2])
title('gMLZ covary nuisance difference over all time scales ','fontsize',18)
xlabel('Center frequency (Hz)')
ylabel('gLZ change from sleep'); ones(1,length(foi));
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_scales_delta_RAW'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'gMLZ_scales_delta_RAW'))

%% gMLZ permutation clustering after regression

clst_LZC_rgt = RMpermclusttest(reg_lzc_wk,reg_lzc_sp,Nperm,p_thresh);


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
        myfigure
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
        figname = sprintf('gMLZ with reg pos clust spect #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
        
        myfigure
        plot_topo_AS(sum(clst_LZC_rgt.lbl_pos==icls)')
        title(sprintf('gMLZ pos clust, p = %1.3f',clst_LZC_rgt.P_val_pos(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(reg_lzc_wk,2)])
        mycolorbar
        figname = sprintf('LZC with reg  pos clust topo #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
    end
    
    % dot plots for all clusters (even if p > 0.05)
    if any(~isnan(clst_LZC_rgt.lbl_pos(:)))            
        dot_slope_plot(reg_lzc_wk,reg_lzc_sp,clst_LZC_rgt.lbl_pos,icls,all_ages,allgenotype,p_thresh)
        title(sprintf('gLZC, d_median = %1.2f',clst_LZC_rgt.d_pos(icls)))
        ylabel('gMLZ')
        %ylim([0.14 0.36])
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s%i.svg',DIRFIGURE,'LZC_regout_pos_clust_dot_plot_#',icls))
        print(gcf,'-dpng','-r300','-painters',sprintf('%s%s%i.png',DIRFIGURE,'LZC_regout_post_clust_dot_plot_#',icls))
    end
end

cnt = 0;
for icls = 1:length(sgcl_neg)
    if sgcl_neg(icls)
        cnt = cnt+1;
        myfigure
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
        figname = sprintf('gMLZ with reg neg clust spect #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
        
        myfigure
        plot_topo_AS(sum(clst_LZC_rgt.lbl_neg==icls)')
        title(sprintf('gMLZ neg clust, p = %1.3f',clst_LZC_rgt.P_val_neg(icls)),'fontsize',18)
        colormap jet
        caxis([0 size(WP,2)])
        mycolorbar
        figname = sprintf('gMLZ with reg neg clust topo #%i',cnt);
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figname))
    end
    
    % dot plots for all clusters (even if p > 0.05)
    if any(~isnan(clst_LZC_rgt.lbl_neg(:)))            
        dot_slope_plot(reg_lzc_wk,reg_lzc_sp,clst_LZC_rgt.lbl_neg,icls,all_ages,allgenotype,p_thresh)
        title(sprintf('gLZC, d_median = %1.2f',clst_LZC_rgt.d_neg(icls)))
        ylabel('gMLZ')
        %ylim([0.14 0.36])
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s%i.svg',DIRFIGURE,'LZC_regout_neg_clust_dot_plot_#',icls))
        print(gcf,'-dpng','-r300','-painters',sprintf('%s%s%i.png',DIRFIGURE,'LZC_regout_neg_clust_dot_plot_#',icls))
    end
end


%%

% Report stats from all time scales after averaging across channels

[h,pval3,ci,stat] = ttest(wk_scl',sp_scl'); %,'tail','right');
Q = mafdr(pval3,'BHFDR',true); % Bejamini Hochberg FDR
ngood = sum(Q < 0.05);


myfigure
plot_topo_AS(-log10(Q)')
title(sprintf('-log10(p), FDR corrected, scales %i -%i',scales(1),scales(end)),'fontsize',18)
colormap jet
caxis([0 8])
mycolorbar
figname = sprintf('sp_wk_log10p_%i_%i',scales(1),scales(end));
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))
%%
myfigure
plot_topo_AS(tvals')
title(sprintf('t-stats, scales %i -%i',scales(1),scales(end)),'fontsize',18)
colormap jet
caxis([0 2])
mycolorbar
figname = sprintf('so_wk_tstats_%i_%i',scales(1),scales(end));
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

myfigure
plot_topo_AS(mean(wk_avg,2))
title(sprintf('Mean MSE values, wake condition, scales %i -%i',scales(1),scales(end)),'fontsize',18)
caxis([0 2.5])
colormap jet
mycolorbar
figname = sprintf('wk_mse_%i_%i',scales(1),scales(end));
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

myfigure
plot_topo_AS(mean(sp_avg,2))
caxis([0 2.5])
mycolorbar
title(sprintf('Mean MSE values, sleep condition, scales %i -%i',scales(1),scales(end)),'fontsize',18)
figname = sprintf('sp_mse_%i_%i',scales(1),scales(end));
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figname))

myfigure
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

myfigure
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

myfigure
hold on
delta = mean(wk_avg) - mean(sp_avg);
scatter(all_ages,delta,'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 1])
[r,p] = corrcoef(all_ages,delta);


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

%% Export table of all clusters
load myjet % custom colormap

% save a figure with the new colorbar
figure
imagesc(randn(10,10))
mycolorbar
colormap(myjet)
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'CustomColorbar'))
%
I = 0;
first_sf = [];
last_sf = [];
min_chan = [];
max_chan = [];
P_c = [];
size_c = [];
perc_c = []; 
pos = [];
reg = [];
d_std = [];
d = [];
whc = cell(1,1);

for icls = 1:length(clst.P_val_pos)
    if any(~isnan(clst.lbl_pos(:)))
        I = I +1;
        % This is not a typo below since MSE is "backwards"
        first_sf(I) = nyqfoi(find(sum(clst.lbl_pos==icls,2),1,'last'));
        last_sf(I) = nyqfoi(find(sum(clst.lbl_pos==icls,2),1,'first'));
        chanmemb = sum(clst.lbl_pos==icls,2);
        chanmemb(chanmemb==0) = [];
        min_chan(I) = min(chanmemb);
        max_chan(I) = max(chanmemb);
        P_c(I) = clst.P_val_pos(icls);
        size_c(I) = clst.n_pos(icls);
        perc_c(I) = size_c(I)/(size(clst.lbl_pos,1)*size(clst.lbl_pos,2));
        pos(I) = true;
        reg(I) = false;
        d_std(I) = clst.d_pos_std(icls);
        d(I) = clst.d_pos(icls);
        whc{I} = 'MSE'
        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        h=subplot(4,4,[5:7, 9:11, 13:15])
        %imagesc(~logical(clst.lbl_pos==icls))
        %imagesc(clst.lbl_pos==icls)
        imagesc(flipud(rot90(abs(squeeze(clst.d_pos_full(icls,:,:))))))
        xticks([1:19])
        xticklabels(labels)
        ylabel('mMSE Timescale')
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        yAX = get(gca,'YAxis');
        box off
        colormap(myjet); caxis([0 3])
        set(yAX,'FontSize', 50)
        set(yAX,'color','k')
        set(gca, 'TickDir', 'out')
        set(gcf,'color','w')
        set(gca,'Layer','top')
        set(xAX,'FontSize', 33)
        xlabel('EEG channel','FontSize',50)
        set(gca,'linewidth',7)
     
        % plot channel projection
        subplot(4,4,1:3); yAX = get(gca,'YAxis'); 
        plot(sum(clst.lbl_pos==icls,1),'k','linewidth',2)
        box off
        ylabel('Count')
        ylim([0 20])
        makefigpretty
        axis normal
        h = gca; h.XAxis.Visible = 'off'; set(yAX,'FontSize', 50)
        figstr = sprintf('Cluster #%i, p = %2.4f, d_{median} = %2.2f',I,clst.P_val_pos(icls),d(I));
        title(figstr,'FontSize',60)
        
        % plot timescale projection
        subplot(4,4,[8 12 16 ]); xAX = get(gca,'XAxis'); 
        plot(sum(clst.lbl_pos==icls,2),20:-1:1,'k','linewidth',4)
        box off
        xlabel('Count')
        xlim([0 20])
        makefigpretty
        axis normal
        h = gca; h.YAxis.Visible = 'off'; set(xAX,'FontSize', 50)
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figstr))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figstr))
    end
end

for icls = 1:length(clst.P_val_neg)
    if any(~isnan(clst.lbl_neg(:)))
        I = I +1;
        % This is not a typo below since MSE is "backwards"
        first_sf(I) = nyqfoi(find(sum(clst.lbl_neg==icls,2),1,'last'));
        last_sf(I) = nyqfoi(find(sum(clst.lbl_neg==icls,2),1,'first'));
        chanmemb = sum(clst.lbl_neg==icls,2);
        chanmemb(chanmemb==0) = [];
        min_chan(I) = min(chanmemb);
        max_chan(I) = max(chanmemb);
        P_c(I) = clst.P_val_neg(icls);
        size_c(I) = clst.n_neg(icls);
        perc_c(I) = size_c(I)/(size(clst.lbl_neg,1)*size(clst.lbl_neg,2));
        pos(I) = false;
        reg(I) = false;
        d_std(I) = clst.d_neg_std(icls);
        d(I) = clst.d_neg(icls);
        whc{I} = 'MSE';
        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        h=subplot(4,4,[5:7, 9:11, 13:15])
        %imagesc(clst.lbl_neg==icls)
        imagesc(flipud(rot90(abs(squeeze(clst.d_neg_full(icls,:,:))))))
        figstr = sprintf('Cluster #%i, p = %2.4f, d_{median} = %2.2f',I,clst.P_val_neg(icls),d(I));
        %colormap(gray)
        xticks([1:19])
        xticklabels(labels)
        ylabel('mMSE Timescale')
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        yAX = get(gca,'YAxis');
        box off
        colormap(myjet); caxis([0 3])
        set(yAX,'FontSize', 50)
        set(yAX,'color','k')
        set(gca, 'TickDir', 'out')
        set(gcf,'color','w')
        set(gca,'Layer','top')
        set(xAX,'FontSize', 33)
        xlabel('EEG channel','FontSize',50)
        set(gca,'linewidth',7)
     
        % plot channel projection
        subplot(4,4,1:3); yAX = get(gca,'YAxis'); 
        plot(sum(clst.lbl_pos==icls,1),'k','linewidth',2)
        box off
        ylabel('Count')
        ylim([0 20])
        set(yAX,'FontSize', 50)
        makefigpretty
        axis normal
        h = gca; h.XAxis.Visible = 'off'; set(yAX,'FontSize', 50)
        title(figstr,'FontSize',60)
        
        % plot timescale projection
        subplot(4,4,[8 12 16 ]); xAX = get(gca,'XAxis'); 
        plot(sum(clst.lbl_pos==icls,2),20:-1:1,'k','linewidth',4)
        box off
        xlabel('Count')
        xlim([0 20])
        set(xAX,'FontSize', 50)
        makefigpretty
        axis normal
        h = gca; h.YAxis.Visible = 'off'; set(xAX,'FontSize', 50)
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figstr))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figstr))
    end
end

for icls = 1:length(clst_LZC.P_val_pos)
    if any(~isnan(clst_LZC.lbl_pos(:)))
        I = I +1;
        last_sf(I) = LZC_foi(find(sum(clst_LZC.lbl_pos==icls,2),1,'last'));
        first_sf(I) = LZC_foi(find(sum(clst_LZC.lbl_pos==icls,2),1,'first'));
        chanmemb = sum(clst_LZC.lbl_pos==icls,2);
        chanmemb(chanmemb==0) = [];
        min_chan(I) = min(chanmemb);
        max_chan(I) = max(chanmemb);
        P_c(I) = clst_LZC.P_val_pos(icls);
        size_c(I) = clst_LZC.n_pos(icls);
        perc_c(I) = size_c(I)/(size(clst_LZC.lbl_pos,1)*size(clst_LZC.lbl_pos,2));
        pos(I) = true;
        reg(I) = false;
        d_std(I) = clst_LZC.d_pos_std(icls);
        d(I) = clst_LZC.d_pos(icls);
        whc{I} = 'gMLZ';
        
        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        h=subplot(4,4,[5:7, 9:11, 13:15])
        %imagesc(clst_LZC.lbl_pos==icls)
        imagesc(flipud(rot90(abs(squeeze(clst_LZC.d_pos_full(icls,:,:))))))
        figstr = sprintf('Cluster #%i, p = %2.4f, d_{median} = %2.2f',I,clst_LZC.P_val_pos(icls),d(I));
        %colormap(gray)
        xticks([1:19])
        xticklabels(labels)
        ylabel('gMLZ Timescale')
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        yAX = get(gca,'YAxis');
        box off
        colormap(myjet); caxis([0 3])
        set(yAX,'FontSize', 50)
        set(yAX,'color','k')
        set(gca, 'TickDir', 'out')
        set(gcf,'color','w')
        set(gca,'Layer','top')
        set(xAX,'FontSize', 33)
        xlabel('EEG channel','FontSize',50)
        set(gca,'linewidth',7)
     
        % plot channel projection
        subplot(4,4,1:3); yAX = get(gca,'YAxis'); 
        plot(sum(clst_LZC.lbl_pos==icls,1),'k','linewidth',2)
        box off
        ylabel('Count')
        ylim([0 20])
        makefigpretty
        axis normal
        h = gca; h.XAxis.Visible = 'off'; set(yAX,'FontSize', 50)
        title(figstr,'FontSize',60)
        
        % plot timescale projection
        subplot(4,4,[8 12 16 ]); xAX = get(gca,'XAxis'); 
        plot(sum(clst_LZC.lbl_pos==icls,2),20:-1:1,'k','linewidth',4)
        box off
        xlabel('Count')
        xlim([0 20])
        makefigpretty
        axis normal
         h = gca; h.YAxis.Visible = 'off'; set(xAX,'FontSize', 50)
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figstr))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figstr))
    end
end

for icls = 1:length(clst_LZC.P_val_neg)
    if any(~isnan(clst_LZC.lbl_neg(:)))
        I = I +1;
        last_sf(I) = LZC_foi(find(sum(clst_LZC.lbl_neg==icls,2),1,'last'));
        first_sf(I) = LZC_foi(find(sum(clst_LZC.lbl_neg==icls,2),1,'first'));
        chanmemb = sum(clst_LZC.lbl_neg==icls,2);
        chanmemb(chanmemb==0) = [];
        min_chan(I) = min(chanmemb);
        max_chan(I) = max(chanmemb);
        P_c(I) = clst_LZC.P_val_neg(icls);
        size_c(I) = clst_LZC.n_neg(icls);
        perc_c(I) = size_c(I)/(size(clst_LZC.lbl_neg,1)*size(clst_LZC.lbl_neg,2));
        pos(I) = false;
        reg(I) = false;
        d_std(I) = clst_LZC.d_neg_std(icls);
        d(I) = clst_LZC.d_neg(icls);
        whc{I} = 'gMLZ';
        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        h=subplot(4,4,[5:7, 9:11, 13:15])
        %imagesc(clst_LZC.lbl_neg==icls)
        imagesc(flipud(rot90(abs(squeeze(clst_LZC.d_neg_full(icls,:,:))))))
        figstr = sprintf('Cluster #%i, p = %2.4f, d_{median} = %2.2f',I,clst_LZC.P_val_neg(icls),d(I));
        %colormap(gray)
        xticks([1:19])
        xticklabels(labels)
        ylabel('gMLZ Timescale')
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        yAX = get(gca,'YAxis');
        box off
        colormap(myjet); caxis([0 3])
        set(yAX,'FontSize', 50)
        set(yAX,'color','k')
        set(gca, 'TickDir', 'out')
        set(gcf,'color','w')
        set(gca,'Layer','top')
        set(xAX,'FontSize', 33)
        xlabel('EEG channel','FontSize',50)
        set(gca,'linewidth',7)
     
        % plot channel projection
        subplot(4,4,1:3); yAX = get(gca,'YAxis'); 
        plot(sum(clst_LZC.lbl_neg==icls,1),'k','linewidth',2)
        box off
        ylabel('Count')
        ylim([0 20])
        makefigpretty
        axis normal
        h = gca; h.XAxis.Visible = 'off'; set(yAX,'FontSize', 50)
        title(figstr,'FontSize',60)
        
        % plot timescale projection
        subplot(4,4,[8 12 16 ]); xAX = get(gca,'XAxis'); 
        plot(sum(clst_LZC.lbl_neg==icls,2),20:-1:1,'k','linewidth',4)
        box off
        xlabel('Count')
        xlim([0 20])
        makefigpretty
        axis normal
         h = gca; h.YAxis.Visible = 'off'; set(xAX,'FontSize', 50)
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figstr))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figstr))
    end
end

for icls = 1:length(clst_rgt.P_val_pos)
    if any(~isnan(clst_rgt.lbl_pos(:)))
        I = I +1;
        % This is not a typo below since MSE is "backwards"
        first_sf(I) = nyqfoi(find(sum(clst_rgt.lbl_pos==icls,2),1,'last'));
        last_sf(I) = nyqfoi(find(sum(clst_rgt.lbl_pos==icls,2),1,'first'));
        chanmemb = sum(clst_rgt.lbl_pos==icls,2);
        chanmemb(chanmemb==0) = [];
        min_chan(I) = min(chanmemb);
        max_chan(I) = max(chanmemb);
        P_c(I) = clst_rgt.P_val_pos(icls);
        size_c(I) = clst_rgt.n_pos(icls);
        perc_c(I) = size_c(I)/(size(clst_rgt.lbl_pos,1)*size(clst_rgt.lbl_pos,2));
        pos(I) = true;
        reg(I) = true;
        d_std(I) = clst_rgt.d_pos_std(icls);
        d(I) = clst_rgt.d_pos(icls);
        whc{I} = 'MSE';
        
        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        h=subplot(4,4,[5:7, 9:11, 13:15])
        %imagesc(clst_rgt.lbl_pos==icls)
        imagesc(flipud(rot90(abs(squeeze(clst_rgt.d_pos_full(icls,:,:))))))
        figstr = sprintf('Cluster #%i, p = %2.4f, d_{median} = %2.2f',I,clst_rgt.P_val_pos(icls),d(I));
        %colormap(gray)
        xticks([1:19])
        xticklabels(labels)
        ylabel('mMSE Timescale')
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        yAX = get(gca,'YAxis');
        box off
        colormap(myjet); caxis([0 3])
        set(yAX,'FontSize', 50)
        set(yAX,'color','k')
        set(gca, 'TickDir', 'out')
        set(gcf,'color','w')
        set(gca,'Layer','top')
        set(xAX,'FontSize', 33)
        xlabel('EEG channel','FontSize',50)
        set(gca,'linewidth',7)
     
        % plot channel projection
        subplot(4,4,1:3); yAX = get(gca,'YAxis'); 
        plot(sum(clst_rgt.lbl_pos==icls,1),'k','linewidth',2)
        box off
        ylabel('Count')
        ylim([0 20])
        makefigpretty
        axis normal
        h = gca; h.XAxis.Visible = 'off'; set(yAX,'FontSize', 50)
        title(figstr,'FontSize',60)
        
        % plot timescale projection
        subplot(4,4,[8 12 16 ]); xAX = get(gca,'XAxis'); 
        plot(sum(clst_rgt.lbl_pos==icls,2),20:-1:1,'k','linewidth',4)
        box off
        xlabel('Count')
        xlim([0 20])
        makefigpretty
        axis normal
         h = gca; h.YAxis.Visible = 'off'; set(xAX,'FontSize', 50)
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figstr))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figstr))
    end
end

for icls = 1:length(clst_rgt.P_val_neg)
    if any(~isnan(clst_rgt.lbl_neg(:)))
        I = I +1;
        % This is not a typo below since MSE is "backwards"
        first_sf(I) = nyqfoi(find(sum(clst_rgt.lbl_neg==icls,2),1,'last'));
        last_sf(I) = nyqfoi(find(sum(clst_rgt.lbl_neg==icls,2),1,'first'));
        chanmemb = sum(clst_rgt.lbl_neg==icls,2);
        chanmemb(chanmemb==0) = [];
        min_chan(I) = min(chanmemb);
        max_chan(I) = max(chanmemb);
        P_c(I) = clst_rgt.P_val_neg(icls);
        size_c(I) = clst_rgt.n_neg(icls);
        perc_c(I) = size_c(I)/(size(clst_rgt.lbl_neg,1)*size(clst_rgt.lbl_neg,2));
        pos(I) = false;
        reg(I) = true;
        d_std(I) = clst_rgt.d_neg_std(icls);
        d(I) = clst_rgt.d_neg(icls);
        whc{I} = 'MSE';
        
        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        h=subplot(4,4,[5:7, 9:11, 13:15])
        %imagesc(clst_rgt.lbl_neg==icls)
        imagesc(flipud(rot90(abs(squeeze(clst_rgt.d_neg_full(icls,:,:))))))
        figstr = sprintf('Cluster #%i, p = %2.4f, d_{median} = %2.2f',I,clst_rgt.P_val_neg(icls),d(I));
        %colormap(gray)
        xticks([1:19])
        xticklabels(labels)
        ylabel('mMSE Timescale')
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        yAX = get(gca,'YAxis');
        box off
        colormap(myjet); caxis([0 3])
        set(yAX,'FontSize', 50)
        set(yAX,'color','k')
        set(gca, 'TickDir', 'out')
        set(gcf,'color','w')
        set(gca,'Layer','top')
        set(xAX,'FontSize', 33)
        xlabel('EEG channel','FontSize',50)
        set(gca,'linewidth',7)
     
        % plot channel projection
        subplot(4,4,1:3); yAX = get(gca,'YAxis'); 
        plot(sum(clst_rgt.lbl_neg==icls,1),'k','linewidth',2)
        box off
        ylabel('Count')
        ylim([0 20])
        makefigpretty
        axis normal
        h = gca; h.XAxis.Visible = 'off'; set(yAX,'FontSize', 50)
        title(figstr,'FontSize',60)
        
        % plot timescale projection
        subplot(4,4,[8 12 16 ]); xAX = get(gca,'XAxis'); 
        plot(sum(clst_rgt.lbl_neg==icls,2),20:-1:1,'k','linewidth',4)
        box off
        xlabel('Count')
        xlim([0 20])
        makefigpretty
        axis normal
         h = gca; h.YAxis.Visible = 'off'; set(xAX,'FontSize', 50)
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figstr))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figstr))
    end
end

for icls = 1:length(clst_LZC_rgt.P_val_pos)
    if any(~isnan(clst_LZC_rgt.lbl_pos(:)))
        I = I +1;
        last_sf(I) = LZC_foi(find(sum(clst_LZC_rgt.lbl_pos==icls,2),1,'last'));
        first_sf(I) = LZC_foi(find(sum(clst_LZC_rgt.lbl_pos==icls,2),1,'first'));
        chanmemb = sum(clst_LZC_rgt.lbl_pos==icls,2);
        chanmemb(chanmemb==0) = [];
        min_chan(I) = min(chanmemb);
        max_chan(I) = max(chanmemb);
        P_c(I) = clst_LZC_rgt.P_val_pos(icls);
        size_c(I) = clst_LZC_rgt.n_pos(icls);
        perc_c(I) = size_c(I)/(size(clst_LZC_rgt.lbl_pos,1)*size(clst_LZC_rgt.lbl_pos,2));
        pos(I) = true;
        reg(I) = true;
        d_std(I) = clst_LZC_rgt.d_pos_std(icls);
        d(I) = clst_LZC_rgt.d_pos(icls);
        whc{I} = 'gMLZ';
        
        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        h=subplot(4,4,[5:7, 9:11, 13:15])
        %imagesc(clst_LZC_rgt.lbl_pos==icls)
        imagesc(flipud(rot90(abs(squeeze(clst_LZC_rgt.d_pos_full(icls,:,:))))))
        figstr = sprintf('Cluster #%i, p = %2.4f, d_{median} = %2.2f',I,clst_LZC_rgt.P_val_pos(icls),d(I));
        %colormap(gray)
        xticks([1:19])
        xticklabels(labels)
        ylabel('gMLZ Timescale')
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        yAX = get(gca,'YAxis');
        box off
        colormap(myjet); caxis([0 3])
        set(yAX,'FontSize', 50)
        set(yAX,'color','k')
        set(gca, 'TickDir', 'out')
        set(gcf,'color','w')
        set(gca,'Layer','top')
        set(xAX,'FontSize', 33)
        xlabel('EEG channel','FontSize',50)
        set(gca,'linewidth',7)
     
        % plot channel projection
        subplot(4,4,1:3); yAX = get(gca,'YAxis'); 
        plot(sum(clst_LZC_rgt.lbl_pos==icls,1),'k','linewidth',2)
        box off
        ylabel('Count')
        ylim([0 20])
        makefigpretty
        axis normal
        h = gca; h.XAxis.Visible = 'off'; set(yAX,'FontSize', 50)
        title(figstr,'FontSize',60)
        
        % plot timescale projection
        subplot(4,4,[8 12 16 ]); xAX = get(gca,'XAxis'); 
        plot(sum(clst_LZC_rgt.lbl_pos==icls,2),20:-1:1,'k','linewidth',4)
        box off
        xlabel('Count')
        xlim([0 20])
        makefigpretty
        axis normal
         h = gca; h.YAxis.Visible = 'off'; set(xAX,'FontSize', 50)
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figstr))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figstr))
    end
end

for icls = 1:length(clst_LZC_rgt.P_val_neg)
    if any(~isnan(clst_LZC_rgt.lbl_neg(:)))
        I = I +1;
        last_sf(I) = LZC_foi(find(sum(clst_LZC_rgt.lbl_neg==icls,2),1,'last'));
        first_sf(I) = LZC_foi(find(sum(clst_LZC_rgt.lbl_neg==icls,2),1,'first'));
        chanmemb = sum(clst_LZC_rgt.lbl_neg==icls,2);
        chanmemb(chanmemb==0) = [];
        min_chan(I) = min(chanmemb);
        max_chan(I) = max(chanmemb);
        P_c(I) = clst_LZC_rgt.P_val_neg(icls);
        size_c(I) = clst_LZC_rgt.n_neg(icls);
        perc_c(I) = size_c(I)/(size(clst_LZC_rgt.lbl_neg,1)*size(clst_LZC_rgt.lbl_neg,2));
        pos(I) = false;
        reg(I) = true;
        d_std(I) = clst_LZC_rgt.d_neg_std(icls);
        d(I) = clst_LZC_rgt.d_neg(icls);
        whc{I} = 'gMLZ';
        
        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        h=subplot(4,4,[5:7, 9:11, 13:15])
        %imagesc(clst_LZC_rgt.lbl_neg==icls)
        imagesc(flipud(rot90(abs(squeeze(clst_LZC_rgt.d_neg_full(icls,:,:))))))
        figstr = sprintf('Cluster #%i, p = %2.4f, d_{median} = %2.2f',I,clst_LZC_rgt.P_val_neg(icls),d(I));
        %colormap(gray)
        xticks([1:19])
        xticklabels(labels)
        ylabel('gMLZ Timescale')
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        yAX = get(gca,'YAxis');
        box off
        colormap(myjet); caxis([0 3])
        set(yAX,'FontSize', 50)
        set(yAX,'color','k')
        set(gca, 'TickDir', 'out')
        set(gcf,'color','w')
        set(gca,'Layer','top')
        set(xAX,'FontSize', 33)
        xlabel('EEG channel','FontSize',50)
        set(gca,'linewidth',7)
     
        % plot channel projection
        subplot(4,4,1:3); yAX = get(gca,'YAxis'); 
        plot(sum(clst_LZC_rgt.lbl_neg==icls,1),'k','linewidth',2)
        box off
        ylabel('Count')
        ylim([0 20])
        makefigpretty
        axis normal
        h = gca; h.XAxis.Visible = 'off'; set(yAX,'FontSize', 50)
        title(figstr,'FontSize',60)
        
        % plot timescale projection
        subplot(4,4,[8 12 16 ]); xAX = get(gca,'XAxis'); 
        plot(sum(clst_LZC_rgt.lbl_neg==icls,2),20:-1:1,'k','linewidth',4)
        box off
        xlabel('Count')
        xlim([0 20])
        makefigpretty
        axis normal
         h = gca; h.YAxis.Visible = 'off'; set(xAX,'FontSize', 50)
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figstr))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figstr))
    end
end

for icls = 1:length(clst_pow.P_val_pos)
    if any(~isnan(clst_pow.lbl_pos(:)))
        I = I +1;
        last_sf(I) = foi(find(sum(clst_pow.lbl_pos==icls,2),1,'last'));
        first_sf(I) = foi(find(sum(clst_pow.lbl_pos==icls,2),1,'first'));
        chanmemb = sum(clst_pow.lbl_pos==icls,2);
        chanmemb(chanmemb==0) = [];
        min_chan(I) = min(chanmemb);
        max_chan(I) = max(chanmemb);
        P_c(I) = clst_pow.P_val_pos(icls);
        size_c(I) = clst_pow.n_pos(icls);
        perc_c(I) = size_c(I)/(size(clst_pow.lbl_pos,1)*size(clst_pow.lbl_pos,2));
        pos(I) = true;
        reg(I) = false;
        d_std(I) = clst_pow.d_pos_std(icls);
        d(I) = clst_pow.d_pos(icls);
        whc{I} = 'power';
        
        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        h=subplot(4,4,[5:7, 9:11, 13:15])
        %imagesc(clst_pow.lbl_pos==icls)
        imagesc(flipud(rot90(abs(squeeze(clst_pow.d_pos_full(icls,:,:))))))
        figstr = sprintf('Cluster #%i, p = %2.4f, d_{median} = %2.2f',I,clst_pow.P_val_pos(icls),d(I));
        %colormap(gray)
        xticks([1:19])
        xticklabels(labels)
        ylabel('Frequency (Hz)')
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        yAX = get(gca,'YAxis');
        yticks([1:8:length(foi)])
        yticklabels(foi([1:8:length(foi)]))
        box off
        colormap(myjet); caxis([0 3])
        set(yAX,'FontSize', 50)
        set(yAX,'color','k')
        set(gca, 'TickDir', 'out')
        set(gcf,'color','w')
        set(gca,'Layer','top')
        set(xAX,'FontSize', 33)
        xlabel('EEG channel','FontSize',50)
        set(gca,'linewidth',7)
     
        % plot channel projection
        subplot(4,4,1:3); yAX = get(gca,'YAxis'); 
        plot(sum(clst_pow.lbl_pos==icls,1),'k','linewidth',2)
        box off
        ylabel('Count')
        ylim([0 40])
        makefigpretty
        axis normal
        h = gca; h.XAxis.Visible = 'off'; set(yAX,'FontSize', 50)
        title(figstr,'FontSize',60)
        
        % plot timescale projection
        subplot(4,4,[8 12 16 ]); xAX = get(gca,'XAxis'); 
        plot(sum(clst_pow.lbl_pos==icls,2),length(foi):-1:1,'k','linewidth',4)
        box off
        xlabel('Count')
        xlim([0 20])
        makefigpretty
        axis normal
        h = gca; h.YAxis.Visible = 'off'; set(xAX,'FontSize', 50)
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figstr))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figstr))
    end
end

for icls = 1:length(clst_pow.P_val_neg)
    if any(~isnan(clst_pow.lbl_neg(:)))
        I = I +1;
        last_sf(I) = foi(find(sum(clst_pow.lbl_neg==icls,2),1,'last'));
        first_sf(I) = foi(find(sum(clst_pow.lbl_neg==icls,2),1,'first'));
        chanmemb = sum(clst_pow.lbl_neg==icls,2);
        chanmemb(chanmemb==0) = [];
        min_chan(I) = min(chanmemb);
        max_chan(I) = max(chanmemb);
        P_c(I) = clst_pow.P_val_neg(icls);
        size_c(I) = clst_pow.n_neg(icls);
        perc_c(I) = size_c(I)/(size(clst_pow.lbl_neg,1)*size(clst_pow.lbl_neg,2));
        pos(I) = false;
        reg(I) = false;
        d_std(I) = clst_pow.d_neg_std(icls);
        d(I) = clst_pow.d_neg(icls);
        whc{I} = 'power';
        
        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        h=subplot(4,4,[5:7, 9:11, 13:15])
        %imagesc(clst_pow.lbl_neg==icls)
        imagesc(flipud(rot90(abs(squeeze(clst_pow.d_neg_full(icls,:,:))))))
        figstr = sprintf('Cluster #%i, p = %2.4f, d_{median} = %2.2f',I,clst_pow.P_val_neg(icls),d(I));
        %colormap(gray)
        xticks([1:19])
        xticklabels(labels)
        ylabel('Frequency (Hz)')
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        yAX = get(gca,'YAxis');
        yticks([1:8:length(foi)])
        yticklabels(foi([1:8:length(foi)]))
        box off
        colormap(myjet); caxis([0 3])
        set(yAX,'FontSize', 50)
        set(yAX,'color','k')
        set(gca, 'TickDir', 'out')
        set(gcf,'color','w')
        set(gca,'Layer','top')
        set(xAX,'FontSize', 33)
        xlabel('EEG channel','FontSize',50)
        set(gca,'linewidth',7)
     
        % plot channel projection
        subplot(4,4,1:3); yAX = get(gca,'YAxis'); 
        plot(sum(clst_pow.lbl_neg==icls,1),'k','linewidth',2)
        box off
        ylabel('Count')
        ylim([0 40])
        makefigpretty
        axis normal
        h = gca; h.XAxis.Visible = 'off'; set(yAX,'FontSize', 50)
        title(figstr,'FontSize',60)
        
        % plot timescale projection
        subplot(4,4,[8 12 16 ]); xAX = get(gca,'XAxis'); 
        plot(sum(clst_pow.lbl_neg==icls,2),length(foi):-1:1,'k','linewidth',4)
        box off
        xlabel('Count')
        xlim([0 20])
        makefigpretty
        axis normal
         h = gca; h.YAxis.Visible = 'off'; set(xAX,'FontSize', 50)
        print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,figstr))
        print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,figstr))
    end
end

% 
% for icls = 1:length(clst_pow_reg.P_val_pos)
%     if any(~isnan(clst_pow_reg.lbl_pos(:)))
%         I = I +1;
%         last_sf(I) = foi(find(sum(clst_pow_reg.lbl_pos==icls,2),1,'last'));
%         first_sf(I) = foi(find(sum(clst_pow_reg.lbl_pos==icls,2),1,'first'));
%         chanmemb = sum(clst_pow_reg.lbl_pos==icls,2);
%         chanmemb(chanmemb==0) = [];
%         min_chan(I) = min(chanmemb);
%         max_chan(I) = max(chanmemb);
%         P_c(I) = clst_pow_reg.P_val_pos(icls);
%         size_c(I) = clst_pow_reg.n_pos(icls);
%         pos(I) = true;
%         reg(I) = true;
%         whc{I} = 'power';
%     end
% end
% 
% for icls = 1:length(clst_pow_reg.P_val_neg)
%     if any(~isnan(clst_pow_reg.lbl_neg(:)))
%         I = I +1;
%         last_sf(I) = foi(find(sum(clst_pow_reg.lbl_neg==icls,2),1,'last'));
%         first_sf(I) = foi(find(sum(clst_pow_reg.lbl_neg==icls,2),1,'first'));
%         chanmemb = sum(clst_pow_reg.lbl_neg==icls,2);
%         chanmemb(chanmemb==0) = [];
%         min_chan(I) = min(chanmemb);
%         max_chan(I) = max(chanmemb);
%         P_c(I) = clst_pow_reg.P_val_neg(icls);
%         size_c(I) = clst_pow_reg.n_neg(icls);
%         pos(I) = false;
%         reg(I) = true;
%         whc{I} = 'power';
%     end
% end

perc_c = perc_c.*100; % express as percentage

switch match
    case true
        strcll = repmat({'Targeted'},1,length(pos));
        cltb = table(whc',strcll',pos',reg',P_c',d',d_std',size_c',perc_c',round(first_sf,1)',...
    round(last_sf,1)',min_chan',max_chan');
        writetable(cltb,sprintf('%s%s',DIRFIGURE,'Stat_cluster_table_targeted.csv'))
    case false
        strcll = repmat({'Broad'},1,length(pos));
        cltb = table(whc',strcll',pos',reg',P_c',d',d_std',size_c',perc_c',round(first_sf,1)',...
    round(last_sf,1)',min_chan',max_chan');
        writetable(cltb,sprintf('%s%s',DIRFIGURE,'Stat_cluster_table_broad.csv'))
end


%% Write data to CSV

SID = [1:size(delta_wk,2)]';

d_wk = mean(delta_wk)';
d_sp = mean(delta_sp)';

for i = 1:size(wk_scl,1)
    v = genvarname(sprintf('%s%i','M_wk_',i));
    eval([v ' = wk_scl(i,:)''']);
end

for i = 1:size(sp_scl,1)
    v = genvarname(sprintf('%s%i','M_sp_',i));
    eval([v ' = sp_scl(i,:)''']);
end

for i = 1:size(wk_scl_lzc,1)
    v = genvarname(sprintf('%s%i','L_wk_',i));
    eval([v ' = wk_scl_lzc(i,:)''']);
end

for i = 1:size(sp_scl_lzc,1)
    v = genvarname(sprintf('%s%i','L_sp_',i));
    eval([v ' = sp_scl_lzc(i,:)''']);
end

delta_TB = table(SID,d_wk,d_sp);

MSE_TB_WK = table(M_wk_1,M_wk_2,M_wk_3,M_wk_4,M_wk_5,M_wk_6,M_wk_7,M_wk_8,...
    M_wk_9,M_wk_10,M_wk_11,M_wk_12,M_wk_13,M_wk_14,M_wk_15,M_wk_16,M_wk_17,M_wk_18,...
    M_wk_19,M_wk_20);

MSE_TB_SP = table(M_sp_1,M_sp_2,M_sp_3,M_sp_4,M_sp_5,M_sp_6,M_sp_7,M_sp_8,...
    M_sp_9,M_sp_10,M_sp_11,M_sp_12,M_sp_13,M_sp_14,M_sp_15,M_sp_16,M_sp_17,M_sp_18,...
    M_sp_19,M_sp_20);

LZC_TB_WK = table(L_wk_1,L_wk_2,L_wk_3,L_wk_4,L_wk_5,L_wk_6,L_wk_7,L_wk_8,...
    L_wk_9,L_wk_10,L_wk_11,L_wk_12,L_wk_13,L_wk_14,L_wk_15,L_wk_16,L_wk_17,L_wk_18,...
    L_wk_19,L_wk_20);

LZC_TB_SP = table(L_sp_1,L_sp_2,L_sp_3,L_sp_4,L_sp_5,L_sp_6,L_sp_7,L_sp_8,...
    L_sp_9,L_sp_10,L_sp_11,L_sp_12,L_sp_13,L_sp_14,L_sp_15,L_sp_16,L_sp_17,L_sp_18,...
    L_sp_19,L_sp_20);

TABLE = delta_TB;

for i = 1:20
    TABLE = [TABLE MSE_TB_WK(:,i) MSE_TB_SP(:,i)];
end

for i = 1:20
    TABLE = [TABLE LZC_TB_WK(:,i) LZC_TB_SP(:,i)];
end

writetable(TABLE,'AS.csv','WriteVariableNames',true)

%% Supplemental figure (mediation analysis)

M = readtable('./MediationAnalysisCIsSimple.xlsx');

myfigure
hold on
patch( [log2(nyqfoi) fliplr(log2(nyqfoi))] , [M.mMSE_lo' fliplr(M.mMSE_hi')], 'b','facealpha', ...
    0.3,'EdgeColor','none')
plot(log2(nyqfoi),zeros(1,length(nyqfoi)),'--k')
legend('95% CI mMSE','fontsize',18)
legend boxoff
%ylim([0 2])
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
ylabel('mediation'); ylim([-0.02 0.04]);
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

print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'MSE_delta_mediation'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'MSE_delta_mediation'))

%%
myfigure
h=subplot(1,1,1);
hold on
patch( [log2(LZC_foi) fliplr(log2(LZC_foi))] , [M.gMLZ_lo' fliplr(M.gMLZ_hi')], 'r','facealpha', ...
    0.3,'EdgeColor','none')
plot(log2(LZC_foi),zeros(1,length(LZC_foi)),'--k')
legend('95% CI gMLZ','fontsize',18)
legend boxoff
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
xlim([0 5])
xlabel('Center frequency (Hz)')
ylabel('mediation'); ylim([-0.02 0.04]); 
makefigpretty
print(gcf,'-dpng',sprintf('%s%s.png',DIRFIGURE,'gMLZ_delta_mediation'))
print(gcf,'-dsvg','-r300','-painters',sprintf('%s%s.svg',DIRFIGURE,'gMLZ_delta_mediation'))

%% Save output

mysave(sprintf('stats_output_%s',date))

%% Generate meds report for reviewers

ReportMeds

%end

%% helper function

function[] = dot_slope_plot(WK,SP,lbl,icls,all_ages,allgenotype,p_thresh)

% Creates plots for AS manuscript showing how individual subjects change
% from wake to sleep

nsubj = 35; % number of subjects to anticipate

% auto detect whether these are power values
if size(WK,2) > 20
    pow = true;
else
    pow = false;
end

myfigure
clsndx = lbl==icls;
clsndx = repmat(clsndx',1,1,35);

avgme = WK(clsndx); % extract only those elements belonging to cluster
ravgme = squeeze(reshape(avgme,1,size(avgme,1)/size(WK,3),size(WK,3))); % rehape so last dim in subjects

if isempty(ravgme)
    fprintf('No points belong to this cluster ... \n')
    return
end

if length(ravgme(:)) == nsubj % if the cluster is only one channel/timescale
    doavg = false; % don't do anymore averaging
else
    doavg = true; % do averaging across channels/timescales
end

switch pow
    case true
        switch doavg
            case true
                W_MSE = log10(nanmean(ravgme))'; % mean across channels/timescales
            case false
                W_MSE = log10(ravgme);
        end
    case false
        switch doavg
            case true
                W_MSE = nanmean(ravgme)'; % mean across channels/timescales
            case false
                W_MSE = ravgme;
        end
end

avgme = SP(clsndx); % extract only those elements belonging to cluster
ravgme = squeeze(reshape(avgme,1,size(avgme,1)/size(SP,3),size(SP,3))); % rehape so last dim in subjects
switch pow
    case true
        switch doavg
            case true
                S_MSE = log10(nanmean(ravgme))'; % mean across channels/timescales
            case false
                S_MSE = log10(ravgme);
        end
    case false
        switch doavg
            case true
                S_MSE = nanmean(ravgme)'; % mean across channels/timescales
            case false
                S_MSE = ravgme;
        end
end

palette = hsv; % the color map we use
step = round(size(palette,1)/5); % spacing between colors
scaling = 1.33; % make color map darker

% select colors for age groupings
cl1 = palette(1,:)./scaling;
cl2 = palette(step,:)./scaling;
cl3 = palette(step*2,:)./scaling;
cl4 = palette(step*3,:)./scaling;
cl5 = palette(step*4,:)./scaling;

% Create age bins
bins = round(logspace(log10(min(all_ages)),log10(max(all_ages)),6));

plotme = [W_MSE S_MSE];
scatlbl = [ones(length(W_MSE),1); ones(length(W_MSE),1).*2];

scatter(scatlbl(repmat(allgenotype==0,1,2)),[plotme(allgenotype==0,1); plotme(allgenotype==0,2)],...
    30,'markeredgecolor','k','markerfacecolor','m','marker','o','linewidth',0.5)
scatter(scatlbl(repmat(allgenotype==1,1,2)),[plotme(allgenotype==1,1); plotme(allgenotype==1,2)],...
    30,'markeredgecolor','k','markerfacecolor','c','marker','o','linewidth',0.5)

% Group by age
age1 = plotme(all_ages >= bins(1) & all_ages < bins(2),:);
age2 = plotme(all_ages >= bins(2) & all_ages < bins(3),:);
age3 = plotme(all_ages >= bins(3) & all_ages < bins(4),:);
age4 = plotme(all_ages >= bins(4) & all_ages < bins(5),:);
age5 = plotme(all_ages >= bins(5) & all_ages < bins(6),:);
plot(age1(1,:)','color',cl1,'linewidth',1.5)
plot(age2(1,:)','color',cl2,'linewidth',1.5)
plot(age3(1,:)','color',cl3,'linewidth',1.5)
plot(age4(1,:)','color',cl4,'linewidth',1.5)
plot(age5(1,:)','color',cl5,'linewidth',1.5)
legend({'non-deletion','deletion',sprintf('%i - %i months',bins(1),bins(2)-1),...
    sprintf('%i - %i months',bins(2),bins(3)-1),sprintf('%i - %i months',bins(3),bins(4)-1),...
    sprintf('%i - %i months',bins(4),bins(5)-1),sprintf('%i - %i months',bins(5),bins(6)-1)},...
    'autoupdate','off','fontsize',18,'location','northeastoutside')
legend boxoff
plot(age1(2:end,:)','color',cl1,'linewidth',1.5)
plot(age2(2:end,:)','color',cl2,'linewidth',1.5)
plot(age3(2:end,:)','color',cl3,'linewidth',1.5)
plot(age4(2:end,:)','color',cl4,'linewidth',1.5)
plot(age5(2:end,:)','color',cl5,'linewidth',1.5)
set(gca,'children',flipud(get(gca,'children'))) % reverse object order so that scatter points are on top
violin([W_MSE S_MSE])

makefigpretty
% if p_thresh >= 0.01
%     axis normal % the figures for p = 0.01 should have long tall dot plot
% end

xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Awake','Asleep'})
ax1 = gca;                   % gca = get current axis
ax1.XAxis.Visible = 'off';   % remove x-axis

end

%% Generate topoplots for reviewer

function[] = ClstAvgTopo(WP,SP,lbl,icls)

NDX = (lbl == icls)'; % which points belong to this cluster

WPavg = squeeze(mean(WP,3));
clstavg_wk = nan(size(WPavg,1),1);
clstmem_wk = cell(size(WPavg,1),1);
for ich = 1:size(clstavg_wk,1)
    if sum(NDX(ich,:)) > 0 % if any of these channels belong to the cluster
        clstavg_wk(ich) = mean(WPavg(ich,logical(NDX(ich,:))));
        clstmem_wk{ich} = WPavg(ich,logical(NDX(ich,:)));
    else
        clstavg_wk(ich) = 0;
        clstmem_wk{ich} = [];
    end
end

myfigure
subplot(2,2,1)
plot_topo_AS_classic(clstavg_wk)
title('awake')
%caxis[-1.5 1.5])
colormap jet
mycolorbar

SPavg = squeeze(mean(SP,3));
clstavg_sp = nan(size(SPavg,1),1);
clstmem_sp = cell(size(SPavg,1),1);
for ich = 1:size(clstavg_sp,1)
    if sum(NDX(ich,:)) > 0 % if any of these channels belong to the cluster
        clstavg_sp(ich) = mean(SPavg(ich,logical(NDX(ich,:))));
        clstmem_sp{ich} = SPavg(ich,logical(NDX(ich,:)));
    else
        clstavg_sp(ich) = 0;
        clstmem_sp{ich} = [];
    end
end

%myfigure
subplot(2,2,2)
plot_topo_AS_classic(clstavg_sp)
title('asleep')
%caxis[-1.5 1.5])
colormap jet
mycolorbar

%myfigure
subplot(2,2,3)
plot_topo_AS_classic(clstavg_wk - clstavg_sp)
title('wake - sleep')
%caxis[-1.5 1.5])
colormap jet
mycolorbar


t = nan(size(WPavg,1),1);
for ich = 1:size(t,1)
    [~,~,~,stat] = ttest(clstmem_wk{ich},clstmem_sp{ich});
    t(ich) = stat.tstat;
end

subplot(2,2,4)
plot_topo_AS_classic(t)
%caxis[-1.5 1.5])
colormap jet
mycolorbar

end