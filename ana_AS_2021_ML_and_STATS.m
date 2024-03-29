%%% STATS 2021 analysis %%%
% Joel Frohlich
% University of California, Los Angeles (UCLA)
% Monti Lab, Psychology Department
%
% University of Tuebingen, Germany
% Institue for Neuromodulation and Neurotechnology 
%
% Last update: 04 Dec, 2022 (cleaned up code and comments)

%%% What does this code do? %%%

% The code first loads in data from files, which contain power and other 
% metrics extracted from the EEG data. The code then initializes arrays to 
% store the data and extracts the data from the files. Next, the code 
% conducts machine learning using the extracted data and produces plots to 
% visualize the results. The code is designed to be used with data from
% patients with Angelman Syndrome (AS), Dup15q Syndrome (DS), and typically 
% developing (TD) children. 

%%% Overview %%%

%     Initializes arrays to store the data and extracts the data from the files.
%     Conducts machine learning using the extracted data.
%     Produces plots to visualize the results.
%     Allows users to choose whether to use linear mixed models (LMMs) or principal component analysis (PCA) for feature selection.
%     Provides a parameter for setting the minimum effect size for feature selection with LMMs.
%     Provides a parameter for setting the minimum amount of variance to be explained by PCs in PCA.
%     Allows users to specify whether to plot the results or not.
%     Uses a random number generator with a fixed seed for reproducibility.
%     Loads data for power, weighted spectral mutual information (wSMI), permutation entropy (PermEn), and directed weighted phase lag index (dwPLI)
%     Uses the above data to train and test a regularized logistic regression classifer 
%     Produces diagnostic plots of the classifier performance.
%     Evaluates the classifer and reports the resulting performance metrics.

% Optional, to clear all variables and close all figures
clearvars
close all

rng(4876) % For reproducibility

% This lay file below is used to produce topoplots 
load angelman_lay.mat lay % lay file 

plotting = true; % do plotting?

% Change lines below for your own directory
% Each of the three variables below contain the locations of the EEG power
% spectra .mat files ... the stem of these file name can then be adapted 
% using sprintf to load other related files 
dup15q_EEG_power_files = dir('./freq_out/Dup15q/nan=0.2/wake/*.mat');
TD_EEG_power_files = dir('./freq_out/TD/nan=0.2/wake/*.mat');
AS_EEG_power_files = dir('./freq_out/AS/nan=0.2/wake/*.mat');

% The variable job tells us whether we are running the script to look
% for features or to do the machine learning
%job = 'FindFeatures';
job = 'MachineLearning';

% The method variables specifies whether to choose features based on LMMs or PCA
%method = 'LMM'
method = 'PCA'

% minimum effect size for feature selection with LMMs
betathresh = 0.5;
% use enough PCs to explain at least this much variance
varthresh = 90; 

% load EEG frequencies of interest
load foi foi
foi = foi(1:find(foi==32)); % only used frequencies up to 32 Hz
nchan = 19; % channels
nfreq = 49; % frequency bins
ntau = 5; % timescales for PermEn and wSMI

%% Extract variables
% This section loads data for power, weighted spectral mutual information 
% (wSMI), permutation entropy (PermEn), and directed weighted phase lag 
% index (dwPLI)


% initialize power

ASWakePow = nan(nchan,nfreq,length(AS_EEG_power_files));
ASSleepPow = nan(nchan,nfreq,length(AS_EEG_power_files));

DSWakePow = nan(nchan,nfreq,length(dup15q_EEG_power_files));
DSSleepPow = nan(nchan,nfreq,length(dup15q_EEG_power_files));

TDWakePow = nan(nchan,nfreq,length(TD_EEG_power_files));
TDSleepPow = nan(nchan,nfreq,length(TD_EEG_power_files));

% initialize wSMI and PermEn

ASWakewSMI = nan(nchan,nchan,ntau,length(AS_EEG_power_files));
ASSleepwSMI = nan(nchan,nchan,ntau,length(AS_EEG_power_files));
ASWakePermEn = nan(nchan,ntau,length(AS_EEG_power_files));
ASSleepPermEn = nan(nchan,ntau,length(AS_EEG_power_files));

DSWakewSMI = nan(nchan,nchan,ntau,length(dup15q_EEG_power_files));
DSSleepwSMI= nan(nchan,nchan,ntau,length(dup15q_EEG_power_files));
DSWakePermEn = nan(nchan,ntau,length(dup15q_EEG_power_files));
DSSleepPermEn = nan(nchan,ntau,length(dup15q_EEG_power_files));

TDWakewSMI = nan(nchan,nchan,ntau,length(TD_EEG_power_files));
TDSleepwSMI = nan(nchan,nchan,ntau,length(TD_EEG_power_files));
TDWakePermEn = nan(nchan,ntau,length(TD_EEG_power_files));
TDSleepPermEn = nan(nchan,ntau,length(TD_EEG_power_files));

% initialize dwPLI

ASWakedwPLI = nan(nchan,nchan,6,length(AS_EEG_power_files));
ASSleepdwPLI = nan(nchan,nchan,6,length(AS_EEG_power_files));

TDWakedwPLI = nan(nchan,nchan,6,length(TD_EEG_power_files));
TDSleepdwPLI = nan(nchan,nchan,6,length(TD_EEG_power_files));

DSWakedwPLI = nan(nchan,nchan,6,length(dup15q_EEG_power_files));
DSSleepdwPLI = nan(nchan,nchan,6,length(dup15q_EEG_power_files));

DSWakedwPLI_SR = nan(length(dup15q_EEG_power_files),length(foi));
DSSleepdwPLI_SR = nan(length(dup15q_EEG_power_files),length(foi));
DSWakedwPLI_LR = nan(length(dup15q_EEG_power_files),length(foi));
DSSleepdwPLI_LR = nan(length(dup15q_EEG_power_files),length(foi));

TDWakedwPLI_SR = nan(length(TD_EEG_power_files),length(foi));
TDSleepdwPLI_SR = nan(length(TD_EEG_power_files),length(foi));
TDWakedwPLI_LR = nan(length(TD_EEG_power_files),length(foi));
TDSleepdwPLI_LR = nan(length(TD_EEG_power_files),length(foi));

swit = nan(length(AS_EEG_power_files),length(foi));
ASSleepdwPLI_SR = nan(length(AS_EEG_power_files),length(foi));
ASWakedwPLI_LR = nan(length(AS_EEG_power_files),length(foi));
ASSleepdwPLI_LR = nan(length(AS_EEG_power_files),length(foi));


% initialize LZc and CTW

TDWakeLZc = nan(nchan,length(TD_EEG_power_files));
TDSleepLZc = nan(nchan,length(TD_EEG_power_files));
DSWakeLZc = nan(nchan,length(dup15q_EEG_power_files));
DSSleepLZc = nan(nchan,length(dup15q_EEG_power_files));
ASWakeLZc = nan(nchan,length(AS_EEG_power_files));
ASSleepLZc = nan(nchan,length(AS_EEG_power_files));

TDWakeCTW = nan(nchan,length(TD_EEG_power_files));
TDSleepCTW = nan(nchan,length(TD_EEG_power_files));
DSWakeCTW = nan(nchan,length(dup15q_EEG_power_files));
DSSleepCTW = nan(nchan,length(dup15q_EEG_power_files));
ASWakeCTW = nan(nchan,length(AS_EEG_power_files));
ASSleepCTW = nan(nchan,length(AS_EEG_power_files));

% initialize mMSE

TDWakeMSE = nan(nchan,length(TD_EEG_power_files));
TDSleepMSE = nan(nchan,length(TD_EEG_power_files));
DSWakeMSE = nan(nchan,length(dup15q_EEG_power_files));
DSSleepMSE = nan(nchan,length(dup15q_EEG_power_files));
ASWakeMSE = nan(nchan,length(AS_EEG_power_files));
ASSleepMSE = nan(nchan,length(AS_EEG_power_files));

% everything else

ASsids = cell(1,length(AS_EEG_power_files));
ASages = nan(1,length(AS_EEG_power_files));
ASsex = nan(1,length(AS_EEG_power_files));
AS15qcn = nan(1,length(AS_EEG_power_files)); % 15q copy number

TDsids = cell(1,length(TD_EEG_power_files));
TDages = nan(1,length(TD_EEG_power_files));
TDsex = nan(1,length(TD_EEG_power_files));
TD15qcn = ones(1,length(TD_EEG_power_files)).*2; % 15q copy number

DSsids = cell(1,length(dup15q_EEG_power_files));
DSages = nan(1,length(dup15q_EEG_power_files));
DSsex = nan(1,length(dup15q_EEG_power_files));
DS15qcn = nan(1,length(dup15q_EEG_power_files)); % 15q copy number

ngood = 15; % min number of valid f Hz wavelet windows
fidx = find(foi == 0.5); % lowest frequency analyzed

minwin = 5; % minimum window size for wSMI
overlap = 0.5; % window overlap for wSMI

% script below identifies short-range and long-range channel pairings
connectivity_distances_plotting

for ifile = 1:length(dup15q_EEG_power_files)
    
    % Load datasets to get other info (age, sex, dup type, etc.)
    tmp = load(sprintf('./postprocessed/Dup15q/%s_pstprc_12-Aug-2021',dup15q_EEG_power_files(ifile).name(1:end-14)));
    DSages(ifile) = tmp.data.age;
    DSsids{ifile} = tmp.data.fstr(8:11);
        switch tmp.data.genotype
            case 'Idic'
                DS15qcn(ifile) = 4;
            case 'Int'
                DS15qcn(ifile) = 3;
        end
    clear tmp
    
    % Load Spectral Power
    waketmp = load(sprintf('./freq_out/Dup15q/nan=0.2/wake/%s',dup15q_EEG_power_files(ifile).name));
    sleeptmp = load(sprintf('./freq_out/Dup15q/nan=0.2/sleep/%s_sleep_freq',dup15q_EEG_power_files(ifile).name(1:end-14)));
    if ~exist('foi','var'), foi = waketmp.foi; end
    if waketmp.n(fidx) >= ngood && sleeptmp.n(fidx) >= ngood  % if there are at least n valid windows for the 0.5 Hz wavelet
        DSWakePow(:,:,ifile) = waketmp.pow;
        DSSleepPow(:,:,ifile) = sleeptmp.pow;
    else
        fprintf('Skipping %s, not enough valid wavelets windows ...\n',dup15q_EEG_power_files(ifile).name)
        continue % go on to next subject or file
    end
    clear waketmp sleeptmp
    
    
    % Load LZc
    waketmp = load(sprintf('./Dup15q_output/LZCv/wake/%s_wake_LZCv_wlen=60s',dup15q_EEG_power_files(ifile).name(1:end-14)));
    sleeptmp = load(sprintf('./Dup15q_output/LZCv/sleep/%s_sleep_LZCv_wlen=60s',dup15q_EEG_power_files(ifile).name(1:end-14)));
    
    assert(~waketmp.LZCout.surrogate,'Wake LZc data is from surrogate')
    assert(~sleeptmp.LZCout.surrogate,'Sleep LZc data is from surrogate')
    
    wakebad = waketmp.LZCout.n_valid_vLZC<2000;
    waketmp.LZCout.vLZC(wakebad) = nan;
    
    sleepbad = sleeptmp.LZCout.n_valid_vLZC<2000;
    sleeptmp.LZCout.vLZC(sleepbad) = nan;
    
    wakebad = waketmp.LZCout.n_valid_CTW<2000;
    waketmp.LZCout.CTW(wakebad) = nan;
    
    sleepbad = sleeptmp.LZCout.n_valid_CTW<2000;
    sleeptmp.LZCout.CTW(sleepbad) = nan;
    
    DSWakeLZc(:,ifile) = nanmean(waketmp.LZCout.vLZC,2);
    DSSleepLZc(:,ifile) =  nanmean(sleeptmp.LZCout.vLZC,2);
    
    DSWakeCTW(:,ifile) = nanmean(squeeze(waketmp.LZCout.CTW),2);
    DSSleepCTW(:,ifile) = nanmean(squeeze(sleeptmp.LZCout.CTW),2);
    
    DSWakeDatLen(ifile) = waketmp.LZCout.dur_used/60; % in minutes
    DSSleepDatLen(ifile) = sleeptmp.LZCout.dur_used/60; % in minutes
    
    DSWakePrcUsed(ifile) = waketmp.LZCout.prc_used*100;
    DSSleepPrcUsed(ifile) = sleeptmp.LZCout.prc_used*100; 
    
%     wakebad = waketmp.LZCout.n_valid_gLZC<2000;
%     waketmp.LZCout.gLZC(wakebad) = nan;
%     
%     sleepbad = sleeptmp.LZCout.n_valid_gLZC<2000;
%     sleeptmp.LZCout.gLZC(sleepbad) = nan;
    
%     % take just delta frequencies from gMLZC (following our last paper)
%     pickme = waketmp.LZCout.cfg.foi >0.99 & waketmp.LZCout.cfg.foi <4;
%     DSWakeLZcclst(:,ifile) = nanmean(waketmp.LZCout.gLZC(:,pickme,:),[2 3]);
%     pickme = sleeptmp.LZCout.cfg.foi >0.99 & sleeptmp.LZCout.cfg.foi <4;
%     DSSleepLZcclst(:,ifile) =  nanmean(sleeptmp.LZCout.gLZC(:,pickme,:),[2 3]);
    
    % Load mMSE
    
    mse_sidx = 1:10; % Use first 10 time-scales, based on findings in Frohlich et al. 2020
    
    waketmp = load(sprintf('./Dup15q_output/Xie/r=0.15/dynr/wake/%s_wake_MSE',dup15q_EEG_power_files(ifile).name(1:end-14)));
    sleeptmp = load(sprintf('./Dup15q_output/Xie/r=0.15/dynr/sleep/%s_sleep_MSE',dup15q_EEG_power_files(ifile).name(1:end-14)));
    
    assert(~waketmp.MSEout.surrogate,'Wake mSME data is from surrogate')
    assert(~sleeptmp.MSEout.surrogate,'Sleep mMSE data is from surrogate')
    assert(waketmp.MSEout.win_len == 30,'Wrong mMSE window length')
    assert(sleeptmp.MSEout.win_len == 30,'Wrong mMSE window length')
    assert(waketmp.MSEout.window_shift == 1,'There should be no mMSE window overlap')
    assert(sleeptmp.MSEout.window_shift == 1,'There should be no mMSE window overlap')
    
    wakebad = waketmp.MSEout.n_valid<100;
    waketmp.MSEout.mse(wakebad) = nan;
    
    sleepbad = sleeptmp.MSEout.n_valid<100;
    sleeptmp.MSEout.mse(sleepbad) = nan;
    
    DSWakeMSE(:,ifile) = nanmean(waketmp.MSEout.mse(:,mse_sidx,:),[2 3]);
    DSSleepMSE(:,ifile) =  nanmean(sleeptmp.MSEout.mse(:,mse_sidx,:),[2 3]);
    
    % Use the timescales/channel identified in the cluster from our last
    % paper (posterior channels, fast timescales)
    
%     DSWakeMSEclst(:,ifile) = nanmean(waketmp.MSEout.mse(postchan,1:10,:),[2 3]);
%     DSSleepMSEclst(:,ifile) =  nanmean(sleeptmp.MSEout.mse(postchan,1:10,:),[2 3]);
    clear waketmp sleeptmp
    
%     % Load criticality
%     waketmp = load(sprintf('./Dup15q_output/CRIT/wake/%s_wake_CRIT',dup15q_EEG_power_files(ifile).name(1:end-14)));
%     sleeptmp = load(sprintf('./Dup15q_output/CRIT/sleep/%s_sleep_CRIT',dup15q_EEG_power_files(ifile).name(1:end-14)));
%     assert(~waketmp.Crit.surrogate,'Wake criticality data is from surrogate')
%     assert(~sleeptmp.Crit.surrogate,'Sleep criticality data is from surrogate')
%     
%     % don't need to look at valid field, because criticality only computed
%     % for windows with all good data
%     DSWakeCrit(:,ifile) = nanmean(waketmp.Crit.Cr,2);
%     DSSleepCrit(:,ifile) = nanmean(sleeptmp.Crit.Cr,2);
%     
%     % Chaoticity (k)
%     DSWakeChaos(:,ifile) = nanmean(waketmp.Crit.K,2);
%     DSSleepChaos(:,ifile) = nanmean(sleeptmp.Crit.K,2);
    
    % Load dwPLI
    
    waketmp = load(sprintf('./Dup15q_output/wPLI/wake/%swake_wPLI_fieldtrip',dup15q_EEG_power_files(ifile).name(1:end-14)));
    sleeptmp = load(sprintf('./Dup15q_output/wPLI/sleep/%ssleep_wPLI_fieldtrip',dup15q_EEG_power_files(ifile).name(1:end-14)));
    
    assert(waketmp.wPLIout.overlap,'No overlap for wPLI windows')
    assert(waketmp.wPLIout.win_len == minwin,'wrong window size')
    assert(sleeptmp.wPLIout.overlap,'No overlap for wPLI windows')
    assert(sleeptmp.wPLIout.win_len == minwin,'wrong window size')   
    
    % Take the integral of each frequency band to get the dwPLI
    flims = 2.^[-1:5];
    for ifrq = 2:length(flims)
        idf = foi >= flims(ifrq-1) & foi < flims(ifrq);
        DSWakedwPLI(:,:,ifrq-1,ifile)  = squeeze(trapz(foi(idf),waketmp.wPLIout.dwPLI(:,:,idf),3));
        DSSleepdwPLI(:,:,ifrq-1,ifile) = squeeze(trapz(foi(idf),sleeptmp.wPLIout.dwPLI(:,:,idf),3));
    end
    
    for ifrq = 1:length(foi)
        % average wake shortrange dwPLI, whole spectrum
        tmp = waketmp.wPLIout.dwPLI(:,:,ifrq);
        DSWakedwPLI_SR(ifile,ifrq) = squeeze(nanmean(tmp(shortrange),[1 2]));
        % average sleep shortrange dwPLI, whole spectrum
        tmp = sleeptmp.wPLIout.dwPLI(:,:,ifrq);
        DSSleepdwPLI_SR(ifile,ifrq) = squeeze(nanmean(tmp(shortrange),[1 2]));
        % average wake longrange dwPLI, whole spectrum
        tmp = waketmp.wPLIout.dwPLI(:,:,ifrq);
        DSWakedwPLI_LR(ifile,ifrq) = squeeze(nanmean(tmp(longrange),[1 2]));
        % average sleep longrange dwPLI, whole spectrum
        tmp = sleeptmp.wPLIout.dwPLI(:,:,ifrq);
        DSSleepdwPLI_LR(ifile,ifrq) = squeeze(nanmean(tmp(longrange),[1 2]));
    end
    
    % Uncomment below to use the average dwPLI rather than the integral
%     DSWakedwPLI(:,:,1,ifile) = waketmp.wPLIout.dwPLI_slow;
%     DSWakedwPLI(:,:,2,ifile) = waketmp.wPLIout.dwPLI_delta1;
%     DSWakedwPLI(:,:,3,ifile) = waketmp.wPLIout.dwPLI_delta2;
%     DSWakedwPLI(:,:,4,ifile) = waketmp.wPLIout.dwPLI_theta;
%     DSWakedwPLI(:,:,5,ifile) = waketmp.wPLIout.dwPLI_alpha;
%     DSWakedwPLI(:,:,6,ifile) = waketmp.wPLIout.dwPLI_beta;
%     
%     DSSleepdwPLI(:,:,1,ifile) = sleeptmp.wPLIout.dwPLI_slow;
%     DSSleepdwPLI(:,:,2,ifile) = sleeptmp.wPLIout.dwPLI_delta1;
%     DSSleepdwPLI(:,:,3,ifile) = sleeptmp.wPLIout.dwPLI_delta2;
%     DSSleepdwPLI(:,:,4,ifile) = sleeptmp.wPLIout.dwPLI_theta;
%     DSSleepdwPLI(:,:,5,ifile) = sleeptmp.wPLIout.dwPLI_alpha;
%     DSSleepdwPLI(:,:,6,ifile) = sleeptmp.wPLIout.dwPLI_beta;
    
    
    % Load wSMI
    waketmp = load(sprintf('./Dup15q_output/wSMI/wake/%s_wake_wSMI',dup15q_EEG_power_files(ifile).name(1:end-14)));
    sleeptmp = load(sprintf('./Dup15q_output/wSMI/sleep/%s_sleep_wSMI',dup15q_EEG_power_files(ifile).name(1:end-14)));
    
    assert(waketmp.wSMIout.win_len == minwin,'Wrong wSMI window size')
    assert(waketmp.wSMIout.window_shift == overlap,'Wrong wSMI window overlap')
    assert(sleeptmp.wSMIout.win_len == minwin,'Wrong wSMI window size')
    assert(sleeptmp.wSMIout.window_shift == overlap,'Wrong wSMI window overlap')
    
    assert(all(waketmp.wSMIout.cfg.cfg.taus_ms==sleeptmp.wSMIout.cfg.cfg.taus_ms),'Kernels lengths don''t match between wake and sleep ')
    assert(all(waketmp.wSMIout.cfg.cfg.taus==sleeptmp.wSMIout.cfg.cfg.taus),'Kernels lengths don''t match between wake and sleep ')
    assert(~waketmp.wSMIout.surrogate,'Wake wSMI data is from surrogate')
    assert(~sleeptmp.wSMIout.surrogate,'Sleep wSMI data is from surrogate')
    
    idx8 = waketmp.wSMIout.cfg.cfg.taus(waketmp.wSMIout.cfg.cfg.taus_ms==8);
    idx16 = waketmp.wSMIout.cfg.cfg.taus(waketmp.wSMIout.cfg.cfg.taus_ms==16);
    idx32 = waketmp.wSMIout.cfg.cfg.taus(waketmp.wSMIout.cfg.cfg.taus_ms==32);
    idx64 = waketmp.wSMIout.cfg.cfg.taus(waketmp.wSMIout.cfg.cfg.taus_ms==64);
    idx128 = waketmp.wSMIout.cfg.cfg.taus(waketmp.wSMIout.cfg.cfg.taus_ms==128);
    
    %     % if the wSMI values are all 0s, skip it
    %     if nanmean(waketmp.wSMIout.wSMI_32) == 0
    %         fprintf('%s does not have wSMI ...\n',dup15q_EEG_power_files(ifile).name)
    %         continue
    %     end
    
    eval(sprintf('DSWakewSMI(:,:,1,ifile) = nanmean(waketmp.wSMIout.wSMI_%i,3);',idx8));
    eval(sprintf('DSWakewSMI(:,:,2,ifile) = nanmean(waketmp.wSMIout.wSMI_%i,3);',idx16));   
    eval(sprintf('DSWakewSMI(:,:,3,ifile) = nanmean(waketmp.wSMIout.wSMI_%i,3);',idx32));
    eval(sprintf('DSWakewSMI(:,:,4,ifile) = nanmean(waketmp.wSMIout.wSMI_%i,3);',idx64));
    eval(sprintf('DSWakewSMI(:,:,5,ifile) = nanmean(waketmp.wSMIout.wSMI_%i,3);',idx128));
    
    eval(sprintf('DSWakePermEn(:,1,ifile) = nanmean(waketmp.wSMIout.PE_%i,2);',idx8));
    eval(sprintf('DSWakePermEn(:,2,ifile) = nanmean(waketmp.wSMIout.PE_%i,2);',idx16));   
    eval(sprintf('DSWakePermEn(:,3,ifile) = nanmean(waketmp.wSMIout.PE_%i,2);',idx32));
    eval(sprintf('DSWakePermEn(:,4,ifile) = nanmean(waketmp.wSMIout.PE_%i,2);',idx64));
    eval(sprintf('DSWakePermEn(:,5,ifile) = nanmean(waketmp.wSMIout.PE_%i,2);',idx128));

    eval(sprintf('DSSleepwSMI(:,:,1,ifile) = nanmean(sleeptmp.wSMIout.wSMI_%i,3);',idx8));
    eval(sprintf('DSSleepwSMI(:,:,2,ifile) = nanmean(sleeptmp.wSMIout.wSMI_%i,3);',idx16));    
    eval(sprintf('DSSleepwSMI(:,:,3,ifile) = nanmean(sleeptmp.wSMIout.wSMI_%i,3);',idx32));
    eval(sprintf('DSSleepwSMI(:,:,4,ifile) = nanmean(sleeptmp.wSMIout.wSMI_%i,3);',idx64));
    eval(sprintf('DSSleepwSMI(:,:,5,ifile) = nanmean(sleeptmp.wSMIout.wSMI_%i,3);',idx128));
    
    eval(sprintf('DSSleepPermEn(:,1,ifile) = nanmean(sleeptmp.wSMIout.PE_%i,2);',idx8));
    eval(sprintf('DSSleepPermEn(:,2,ifile) = nanmean(sleeptmp.wSMIout.PE_%i,2);',idx16));
    eval(sprintf('DSSleepPermEn(:,3,ifile) = nanmean(sleeptmp.wSMIout.PE_%i,2);',idx32));
    eval(sprintf('DSSleepPermEn(:,4,ifile) = nanmean(sleeptmp.wSMIout.PE_%i,2);',idx64));
    eval(sprintf('DSSleepPermEn(:,5,ifile) = nanmean(sleeptmp.wSMIout.PE_%i,2);',idx128));
    
    
    clear waketmp sleeptmp
    
end


for ifile = 1:length(TD_EEG_power_files)
    
    % Load datasets to get other info (age, sex, dup type, etc.)
    tmp = load(sprintf('./postprocessed/TD/%s_pstprc_12-Aug-2021',TD_EEG_power_files(ifile).name(1:end-14)));
    TDages(ifile) = tmp.data.age;
    if ~isnan(str2double(tmp.data.fstr(1:2)))
        % if it's a 2 digit SID
        TDsids{ifile} = tmp.data.fstr(1:2);
    else
        % if it's a 1 digit SID
        TDsids{ifile} = tmp.data.fstr(1);
    end
    assert(strcmp(tmp.data.sex,'F') || strcmp(tmp.data.sex,'M'),'Error reading sex')
    TDsex(ifile) = strcmp(tmp.data.sex,'F'); % isfemale
    
    waketmp = load(sprintf('./freq_out/TD/nan=0.2/wake/%s',TD_EEG_power_files(ifile).name));
    sleeptmp = load(sprintf('./freq_out/TD/nan=0.2/sleep/%s_sleep_freq',TD_EEG_power_files(ifile).name(1:end-14)));
    assert(all(waketmp.foi==foi),'frequency bins don''t match between conditions')
    if waketmp.n(fidx) >= ngood  && sleeptmp.n(fidx) >= ngood  % if there are at least n valid windows for the 0.5 Hz wavelet
        TDWakePow(:,:,ifile) = waketmp.pow;
        TDSleepPow(:,:,ifile) = sleeptmp.pow;
    else
        fprintf('Skipping %s, not enough valid wavelets windows ...\n',TD_EEG_power_files(ifile).name)
        continue % go on to next subject or file
    end
    clear waketmp sleeptmp
    
    % Load LZc
    waketmp = load(sprintf('./TD_output/LZCv/wake/%s_wake_LZCv_wlen=60s',TD_EEG_power_files(ifile).name(1:end-14)));
    sleeptmp = load(sprintf('./TD_output/LZCv/sleep/%s_sleep_LZCv_wlen=60s',TD_EEG_power_files(ifile).name(1:end-14)));
    
    assert(~waketmp.LZCout.surrogate,'Wake LZc data is from surrogate')
    assert(~sleeptmp.LZCout.surrogate,'Sleep LZc data is from surrogate')
    
    wakebad = waketmp.LZCout.n_valid_vLZC<2000;
    waketmp.LZCout.vLZC(wakebad) = nan;
    
    sleepbad = sleeptmp.LZCout.n_valid_vLZC<2000;
    sleeptmp.LZCout.vLZC(sleepbad) = nan;
    
    wakebad = waketmp.LZCout.n_valid_CTW<2000;
    waketmp.LZCout.CTW(wakebad) = nan;
    
    sleepbad = sleeptmp.LZCout.n_valid_CTW<2000;
    sleeptmp.LZCout.CTW(sleepbad) = nan;
    
    TDWakeLZc(:,ifile) = nanmean(waketmp.LZCout.vLZC,2);
    TDSleepLZc(:,ifile) =  nanmean(sleeptmp.LZCout.vLZC,2);
    
    TDWakeCTW(:,ifile) = nanmean(waketmp.LZCout.CTW,2);
    TDSleepCTW(:,ifile) = nanmean(sleeptmp.LZCout.CTW,2);
    
    TDWakeDatLen(ifile) = waketmp.LZCout.dur_used/60; % in minutes
    TDSleepDatLen(ifile) = sleeptmp.LZCout.dur_used/60; % in minutes
    
    TDWakePrcUsed(ifile) = waketmp.LZCout.prc_used*100;
    TDSleepPrcUsed(ifile) = sleeptmp.LZCout.prc_used*100; 
    
    
    clear waketmp sleeptmp
    
    % Load mMSE
    
    waketmp = load(sprintf('./TD_output/Xie/r=0.15/dynr/wake/%s_wake_MSE',TD_EEG_power_files(ifile).name(1:end-14)));
    sleeptmp = load(sprintf('./TD_output/Xie/r=0.15/dynr/sleep/%s_sleep_MSE',TD_EEG_power_files(ifile).name(1:end-14)));
    
    assert(~waketmp.MSEout.surrogate,'Wake mSME data is from surrogate')
    assert(~sleeptmp.MSEout.surrogate,'Sleep mMSE data is from surrogate')
    assert(waketmp.MSEout.win_len == 30,'Wrong mMSE window length')
    assert(sleeptmp.MSEout.win_len == 30,'Wrong mMSE window length')
    assert(waketmp.MSEout.window_shift == 1,'There should be no mMSE window overlap')
    assert(sleeptmp.MSEout.window_shift == 1,'There should be no mMSE window overlap')
    
    wakebad = waketmp.MSEout.n_valid<100;
    waketmp.MSEout.mse(wakebad) = nan;
    
    sleepbad = sleeptmp.MSEout.n_valid<100;
    sleeptmp.MSEout.mse(sleepbad) = nan;
    
    TDWakeMSE(:,ifile) = nanmean(waketmp.MSEout.mse(:,mse_sidx,:),[2 3]);
    TDSleepMSE(:,ifile) =  nanmean(sleeptmp.MSEout.mse(:,mse_sidx,:),[2 3]);
    
    clear waketmp sleeptmp   
    
    % Load dwPLI
    
    waketmp = load(sprintf('./TD_output/wPLI/wake/%swake_wPLI_fieldtrip',TD_EEG_power_files(ifile).name(1:end-14)));
    sleeptmp = load(sprintf('./TD_output/wPLI/sleep/%ssleep_wPLI_fieldtrip',TD_EEG_power_files(ifile).name(1:end-14)));
    
    assert(waketmp.wPLIout.overlap,'No overlap for wPLI windows')
    assert(waketmp.wPLIout.win_len == minwin,'wrong window size')
    assert(sleeptmp.wPLIout.overlap,'No overlap for wPLI windows')
    assert(sleeptmp.wPLIout.win_len == minwin,'wrong window size')
    
    % Take the integral of each frequency band to get the dwPLI
    flims = 2.^[-1:5];
    for ifrq = 2:length(flims)
        idf = foi >= flims(ifrq-1) & foi < flims(ifrq);
        TDWakedwPLI(:,:,ifrq-1,ifile)  = squeeze(trapz(foi(idf),waketmp.wPLIout.dwPLI(:,:,idf),3));
        TDSleepdwPLI(:,:,ifrq-1,ifile) = squeeze(trapz(foi(idf),sleeptmp.wPLIout.dwPLI(:,:,idf),3));
    end
    
    for ifrq = 1:length(foi)
        % average wake shortrange dwPLI, whole spectrum
        tmp = waketmp.wPLIout.dwPLI(:,:,ifrq);
        TDWakedwPLI_SR(ifile,ifrq) = squeeze(nanmean(tmp(shortrange),[1 2]));
        % average sleep shortrange dwPLI, whole spectrum
        tmp = sleeptmp.wPLIout.dwPLI(:,:,ifrq);
        TDSleepdwPLI_SR(ifile,ifrq) = squeeze(nanmean(tmp(shortrange),[1 2]));
        % average wake longrange dwPLI, whole spectrum
        tmp = waketmp.wPLIout.dwPLI(:,:,ifrq);
        TDWakedwPLI_LR(ifile,ifrq) = squeeze(nanmean(tmp(longrange),[1 2]));
        % average sleep longrange dwPLI, whole spectrum
        tmp = sleeptmp.wPLIout.dwPLI(:,:,ifrq);
        TDSleepdwPLI_LR(ifile,ifrq) = squeeze(nanmean(tmp(longrange),[1 2]));
    end
    
    % Uncomment below to use average rather than integral
%     TDWakedwPLI(:,:,1,ifile) = waketmp.wPLIout.dwPLI_slow;
%     TDWakedwPLI(:,:,2,ifile) = waketmp.wPLIout.dwPLI_delta1;
%     TDWakedwPLI(:,:,3,ifile) = waketmp.wPLIout.dwPLI_delta2;
%     TDWakedwPLI(:,:,4,ifile) = waketmp.wPLIout.dwPLI_theta;
%     TDWakedwPLI(:,:,5,ifile) = waketmp.wPLIout.dwPLI_alpha;
%     TDWakedwPLI(:,:,6,ifile) = waketmp.wPLIout.dwPLI_beta;
%     
%     TDSleepdwPLI(:,:,1,ifile) = sleeptmp.wPLIout.dwPLI_slow;
%     TDSleepdwPLI(:,:,2,ifile) = sleeptmp.wPLIout.dwPLI_delta1;
%     TDSleepdwPLI(:,:,3,ifile) = sleeptmp.wPLIout.dwPLI_delta2;
%     TDSleepdwPLI(:,:,4,ifile) = sleeptmp.wPLIout.dwPLI_theta;
%     TDSleepdwPLI(:,:,5,ifile) = sleeptmp.wPLIout.dwPLI_alpha;
%     TDSleepdwPLI(:,:,6,ifile) = sleeptmp.wPLIout.dwPLI_beta;
    
    
    
    % Load wSMI
    waketmp = load(sprintf('./TD_output/wSMI/wake/%s_wake_wSMI',TD_EEG_power_files(ifile).name(1:end-14)));
    sleeptmp = load(sprintf('./TD_output/wSMI/sleep/%s_sleep_wSMI',TD_EEG_power_files(ifile).name(1:end-14)));
    
    assert(waketmp.wSMIout.win_len == minwin,'Wrong wSMI window size')
    assert(waketmp.wSMIout.window_shift == overlap,'Wrong wSMI window overlap')
    assert(sleeptmp.wSMIout.win_len == minwin,'Wrong wSMI window size')
    assert(sleeptmp.wSMIout.window_shift == overlap,'Wrong wSMI window overlap')
    
    assert(all(waketmp.wSMIout.cfg.cfg.taus_ms==sleeptmp.wSMIout.cfg.cfg.taus_ms),'Kernels lengths don''t match between wake and sleep ')
    assert(all(waketmp.wSMIout.cfg.cfg.taus==sleeptmp.wSMIout.cfg.cfg.taus),'Kernels lengths don''t match between wake and sleep ')
    assert(~waketmp.wSMIout.surrogate,'Wake wSMI data is from surrogate')
    assert(~sleeptmp.wSMIout.surrogate,'Sleep wSMI data is from surrogate')
    
    idx8 = waketmp.wSMIout.cfg.cfg.taus(waketmp.wSMIout.cfg.cfg.taus_ms==8);
    idx16 = waketmp.wSMIout.cfg.cfg.taus(waketmp.wSMIout.cfg.cfg.taus_ms==16);
    idx32 = waketmp.wSMIout.cfg.cfg.taus(waketmp.wSMIout.cfg.cfg.taus_ms==32);
    idx64 = waketmp.wSMIout.cfg.cfg.taus(waketmp.wSMIout.cfg.cfg.taus_ms==64);
    idx128 = waketmp.wSMIout.cfg.cfg.taus(waketmp.wSMIout.cfg.cfg.taus_ms==128);

    
    eval(sprintf('TDWakewSMI(:,:,1,ifile) = nanmean(waketmp.wSMIout.wSMI_%i,3);',idx8));
    eval(sprintf('TDWakewSMI(:,:,2,ifile) = nanmean(waketmp.wSMIout.wSMI_%i,3);',idx16));
    eval(sprintf('TDWakewSMI(:,:,3,ifile) = nanmean(waketmp.wSMIout.wSMI_%i,3);',idx32));
    eval(sprintf('TDWakewSMI(:,:,4,ifile) = nanmean(waketmp.wSMIout.wSMI_%i,3);',idx64));
    eval(sprintf('TDWakewSMI(:,:,5,ifile) = nanmean(waketmp.wSMIout.wSMI_%i,3);',idx128));
    
    eval(sprintf('TDWakePermEn(:,1,ifile) = nanmean(waketmp.wSMIout.PE_%i,2);',idx8));
    eval(sprintf('TDWakePermEn(:,2,ifile) = nanmean(waketmp.wSMIout.PE_%i,2);',idx16));   
    eval(sprintf('TDWakePermEn(:,3,ifile) = nanmean(waketmp.wSMIout.PE_%i,2);',idx32));
    eval(sprintf('TDWakePermEn(:,4,ifile) = nanmean(waketmp.wSMIout.PE_%i,2);',idx64));
    eval(sprintf('TDWakePermEn(:,5,ifile) = nanmean(waketmp.wSMIout.PE_%i,2);',idx128));
    
    eval(sprintf('TDSleepwSMI(:,:,1,ifile) = nanmean(sleeptmp.wSMIout.wSMI_%i,3);',idx8));    
    eval(sprintf('TDSleepwSMI(:,:,2,ifile) = nanmean(sleeptmp.wSMIout.wSMI_%i,3);',idx16));    
    eval(sprintf('TDSleepwSMI(:,:,3,ifile) = nanmean(sleeptmp.wSMIout.wSMI_%i,3);',idx32));
    eval(sprintf('TDSleepwSMI(:,:,4,ifile) = nanmean(sleeptmp.wSMIout.wSMI_%i,3);',idx64));
    eval(sprintf('TDSleepwSMI(:,:,5,ifile) = nanmean(sleeptmp.wSMIout.wSMI_%i,3);',idx128));

    eval(sprintf('TDSleepPermEn(:,1,ifile) = nanmean(sleeptmp.wSMIout.PE_%i,2);',idx8));
    eval(sprintf('TDSleepPermEn(:,2,ifile) = nanmean(sleeptmp.wSMIout.PE_%i,2);',idx16));
    eval(sprintf('TDSleepPermEn(:,3,ifile) = nanmean(sleeptmp.wSMIout.PE_%i,2);',idx32));
    eval(sprintf('TDSleepPermEn(:,4,ifile) = nanmean(sleeptmp.wSMIout.PE_%i,2);',idx64));
    eval(sprintf('TDSleepPermEn(:,5,ifile) = nanmean(sleeptmp.wSMIout.PE_%i,2);',idx128));
    
    clear waketmp sleeptmp
end
%%
Tages = readtable('./AS_table_2019_JF.csv');
for ifile = 1:length(AS_EEG_power_files)
    nsur = 20; % number of surrogates for the LZC decomp
    seclen = 256*60; % length EEG sections to use for LZC decomp
    
    % Load datasets to get other info (age, sex, dup type, etc.)
    tmp = load(sprintf('./postprocessed/AS/%spstprc_12-Aug-2021',AS_EEG_power_files(ifile).name(1:end-14)));

    us = find(tmp.data.cfg.fileout=='_'); % find underscores
    slsh = find(tmp.data.cfg.fileout=='/'); % find slashes
    fd = str2double(tmp.data.cfg.fileout(us(1)+1:us(2)-1)); % file date
    if ~isfield(tmp.data,'ID')
        if isempty(slsh)
            tmp.data.ID = tmp.data.cfg.fileout(1:us(1)-1);
        else
            tmp.data.ID = tmp.data.cfg.fileout(slsh(end)+1:us(1)-1);
        end
    end
    try
        assert(isa(tmp.data.ID,'double'))
        rw = find(tmp.data.ID == Tages.ID & fd == Tages.TimeRecording);
    catch
        rw = find(str2double(tmp.data.ID) == Tages.ID & fd == Tages.TimeRecording);
    end
    ASages(ifile) = Tages.AgeRecording(rw);
    ASsex(ifile) = strcmp(Tages.Gender(rw),'F'); % Is female?
    if strcmpi(Tages.Deletion(rw),'Yes')
        AS15qcn(ifile) = 1;
    elseif strcmpi(Tages.Deletion(rw),'No')
        AS15qcn(ifile) = 2;
    end


    ASsids{ifile} = tmp.data.fstr(1:6);
    
    
    try
        waketmp = load(sprintf('./freq_out/AS/nan=0.2/wake/%s',AS_EEG_power_files(ifile).name));
        sleeptmp = load(sprintf('./freq_out/AS/nan=0.2/sleep/%s_sleep_freq',AS_EEG_power_files(ifile).name(1:end-14)));
    catch
        % DO NOTHING AND SKIP
        AS_EEG_power_files(ifile).name
        continue
    end
    assert(all(waketmp.foi==foi),'frequency bins don''t match between conditions')
    
    if waketmp.n(fidx) >= ngood  && sleeptmp.n(fidx) >= ngood  % if there are at least n valid windows for the 0.5 Hz wavelet
        ASWakePow(:,:,ifile) = waketmp.pow;
        ASSleepPow(:,:,ifile) = sleeptmp.pow;
    else
        fprintf('Skipping %s, not enough valid wavelets windows ...\n',AS_EEG_power_files(ifile).name)
        continue % go on to next subject or file
    end
    
    clear waketmp sleeptmp
    
    % Load LZc
    try
        waketmp = load(sprintf('./AS_output/LZCv/wake/%s_wake_LZCv_wlen=60s',AS_EEG_power_files(ifile).name(1:end-14)));
        sleeptmp = load(sprintf('./AS_output/LZCv/sleep/%s_sleep_LZCv_wlen=60s',AS_EEG_power_files(ifile).name(1:end-14)));
    catch
        AS_EEG_power_files(ifile).name
        fprintf('Missing LZc, skipping complexity measures for this one ... \n')
        continue
    end
    
    assert(~waketmp.LZCout.surrogate,'Wake LZc data is from surrogate')
    assert(~sleeptmp.LZCout.surrogate,'Sleep LZc data is from surrogate')
    
    wakebad = waketmp.LZCout.n_valid_vLZC<2000;
    waketmp.LZCout.vLZC(wakebad) = nan;
    
    sleepbad = sleeptmp.LZCout.n_valid_vLZC<2000;
    sleeptmp.LZCout.vLZC(sleepbad) = nan;
    
    wakebad = waketmp.LZCout.n_valid_CTW<2000;
    waketmp.LZCout.CTW(wakebad) = nan;
    
    sleepbad = sleeptmp.LZCout.n_valid_CTW<2000;
    sleeptmp.LZCout.CTW(sleepbad) = nan;
    
    ASWakeLZc(:,ifile) = nanmean(waketmp.LZCout.vLZC,2);
    ASSleepLZc(:,ifile) =  nanmean(sleeptmp.LZCout.vLZC,2);
    
    ASWakeCTW(:,ifile) = nanmean(waketmp.LZCout.CTW,2);
    ASSleepCTW(:,ifile) = nanmean(sleeptmp.LZCout.CTW,2);
    
    ASWakeDatLen(ifile) = waketmp.LZCout.dur_used/60; % in minutes
    ASSleepDatLen(ifile) = sleeptmp.LZCout.dur_used/60; % in minutes
    
    ASWakePrcUsed(ifile) = waketmp.LZCout.prc_used*100;
    ASSleepPrcUsed(ifile) = sleeptmp.LZCout.prc_used*100; 

    clear waketmp sleeptmp
    
    % Load mMSE
    
    waketmp = load(sprintf('./AS_output/Xie/r=0.15/dynr/wake/%s_wake_MSE',AS_EEG_power_files(ifile).name(1:end-14)));
    sleeptmp = load(sprintf('./AS_output/Xie/r=0.15/dynr/sleep/%s_sleep_MSE',AS_EEG_power_files(ifile).name(1:end-14)));
    
    assert(~waketmp.MSEout.surrogate,'Wake mSME data is from surrogate')
    assert(~sleeptmp.MSEout.surrogate,'Sleep mMSE data is from surrogate')
    assert(waketmp.MSEout.win_len == 30,'Wrong mMSE window length')
    assert(sleeptmp.MSEout.win_len == 30,'Wrong mMSE window length')
    assert(waketmp.MSEout.window_shift == 1,'There should be no mMSE window overlap')
    assert(sleeptmp.MSEout.window_shift == 1,'There should be no mMSE window overlap')
    
    wakebad = waketmp.MSEout.n_valid<100;
    waketmp.MSEout.mse(wakebad) = nan;
    
    sleepbad = sleeptmp.MSEout.n_valid<100;
    sleeptmp.MSEout.mse(sleepbad) = nan;
       
    ASWakeMSE(:,ifile) = nanmean(waketmp.MSEout.mse(:,mse_sidx,:),[2 3]);
    ASSleepMSE(:,ifile) =  nanmean(sleeptmp.MSEout.mse(:,mse_sidx,:),[2 3]);
    
    clear waketmp sleeptmp
     
    %Load dwPLI
    
    waketmp = load(sprintf('./AS_output/wPLI/wake/%swake_wPLI_fieldtrip',AS_EEG_power_files(ifile).name(1:end-14)));
    sleeptmp = load(sprintf('./AS_output/wPLI/sleep/%ssleep_wPLI_fieldtrip',AS_EEG_power_files(ifile).name(1:end-14)));
    
    assert(waketmp.wPLIout.overlap,'No overlap for wPLI windows')
    assert(waketmp.wPLIout.win_len == minwin,'wrong window size')
    assert(sleeptmp.wPLIout.overlap,'No overlap for wPLI windows')
    assert(sleeptmp.wPLIout.win_len == minwin,'wrong window size')
    
    % Take the integral of each frequency band to get the dwPLI
    flims = 2.^[-1:5];
    for ifrq = 2:length(flims)
        idf = foi >= flims(ifrq-1) & foi < flims(ifrq);
        ASWakedwPLI(:,:,ifrq-1,ifile)  = squeeze(trapz(foi(idf),waketmp.wPLIout.dwPLI(:,:,idf),3));
        ASSleepdwPLI(:,:,ifrq-1,ifile) = squeeze(trapz(foi(idf),sleeptmp.wPLIout.dwPLI(:,:,idf),3));
    end
    
    for ifrq = 1:length(foi)
        % average wake shortrange dwPLI, whole spectrum
        tmp = waketmp.wPLIout.dwPLI(:,:,ifrq);
        ASWakedwPLI_SR(ifile,ifrq) = squeeze(nanmean(tmp(shortrange),[1 2]));
        % average sleep shortrange dwPLI, whole spectrum
        tmp = sleeptmp.wPLIout.dwPLI(:,:,ifrq);
        ASSleepdwPLI_SR(ifile,ifrq) = squeeze(nanmean(tmp(shortrange),[1 2]));
        % average wake longrange dwPLI, whole spectrum
        tmp = waketmp.wPLIout.dwPLI(:,:,ifrq);
        ASWakedwPLI_LR(ifile,ifrq) = squeeze(nanmean(tmp(longrange),[1 2]));
        % average sleep longrange dwPLI, whole spectrum
        tmp = sleeptmp.wPLIout.dwPLI(:,:,ifrq);
        ASSleepdwPLI_LR(ifile,ifrq) = squeeze(nanmean(tmp(longrange),[1 2]));
    end
    
    % Uncomment below to use the average rather than the intgral 
%     ASWakedwPLI(:,:,1,ifile) = waketmp.wPLIout.dwPLI_slow;
%     ASWakedwPLI(:,:,2,ifile) = waketmp.wPLIout.dwPLI_delta1;
%     ASWakedwPLI(:,:,3,ifile) = waketmp.wPLIout.dwPLI_delta2;
%     ASWakedwPLI(:,:,4,ifile) = waketmp.wPLIout.dwPLI_theta;
%     ASWakedwPLI(:,:,5,ifile) = waketmp.wPLIout.dwPLI_alpha;
%     ASWakedwPLI(:,:,6,ifile) = waketmp.wPLIout.dwPLI_beta;
%     
%     ASSleepdwPLI(:,:,1,ifile) = sleeptmp.wPLIout.dwPLI_slow;
%     ASSleepdwPLI(:,:,2,ifile) = sleeptmp.wPLIout.dwPLI_delta1;
%     ASSleepdwPLI(:,:,3,ifile) = sleeptmp.wPLIout.dwPLI_delta2;
%     ASSleepdwPLI(:,:,4,ifile) = sleeptmp.wPLIout.dwPLI_theta;
%     ASSleepdwPLI(:,:,5,ifile) = sleeptmp.wPLIout.dwPLI_alpha;
%     ASSleepdwPLI(:,:,6,ifile) = sleeptmp.wPLIout.dwPLI_beta;
    
    
    % Load wSMI
        waketmp = load(sprintf('./AS_output/wSMI/wake/%s_wake_wSMI',AS_EEG_power_files(ifile).name(1:end-14)));
        sleeptmp = load(sprintf('./AS_output/wSMI/sleep/%s_sleep_wSMI',AS_EEG_power_files(ifile).name(1:end-14)));

    
    assert(waketmp.wSMIout.win_len == minwin,'Wrong wSMI window size')
    assert(waketmp.wSMIout.window_shift == overlap,'Wrong wSMI window overlap')
    assert(sleeptmp.wSMIout.win_len == minwin,'Wrong wSMI window size')
    assert(sleeptmp.wSMIout.window_shift == overlap,'Wrong wSMI window overlap')
    
    assert(all(waketmp.wSMIout.cfg.cfg.taus_ms==sleeptmp.wSMIout.cfg.cfg.taus_ms),'Kernels lengths don''t match between wake and sleep ')
    assert(all(waketmp.wSMIout.cfg.cfg.taus==sleeptmp.wSMIout.cfg.cfg.taus),'Kernels lengths don''t match between wake and sleep ')
    assert(~waketmp.wSMIout.surrogate,'Wake wSMI data is from surrogate')
    assert(~sleeptmp.wSMIout.surrogate,'Sleep wSMI data is from surrogate')

    idx8 = waketmp.wSMIout.cfg.cfg.taus(waketmp.wSMIout.cfg.cfg.taus_ms==8);
    idx16 = waketmp.wSMIout.cfg.cfg.taus(waketmp.wSMIout.cfg.cfg.taus_ms==16);    
    idx32 = waketmp.wSMIout.cfg.cfg.taus(waketmp.wSMIout.cfg.cfg.taus_ms==32);
    idx64 = waketmp.wSMIout.cfg.cfg.taus(waketmp.wSMIout.cfg.cfg.taus_ms==64);
    idx128 = waketmp.wSMIout.cfg.cfg.taus(waketmp.wSMIout.cfg.cfg.taus_ms==128);

    
    eval(sprintf('ASWakewSMI(:,:,1,ifile) = nanmean(waketmp.wSMIout.wSMI_%i,3);',idx8));
    eval(sprintf('ASWakewSMI(:,:,2,ifile) = nanmean(waketmp.wSMIout.wSMI_%i,3);',idx16));
    eval(sprintf('ASWakewSMI(:,:,3,ifile) = nanmean(waketmp.wSMIout.wSMI_%i,3);',idx32));
    eval(sprintf('ASWakewSMI(:,:,4,ifile) = nanmean(waketmp.wSMIout.wSMI_%i,3);',idx64));
    eval(sprintf('ASWakewSMI(:,:,5,ifile) = nanmean(waketmp.wSMIout.wSMI_%i,3);',idx128));
    
    eval(sprintf('ASWakePermEn(:,1,ifile) = nanmean(waketmp.wSMIout.PE_%i,2);',idx8));
    eval(sprintf('ASWakePermEn(:,2,ifile) = nanmean(waketmp.wSMIout.PE_%i,2);',idx16));   
    eval(sprintf('ASWakePermEn(:,3,ifile) = nanmean(waketmp.wSMIout.PE_%i,2);',idx32));
    eval(sprintf('ASWakePermEn(:,4,ifile) = nanmean(waketmp.wSMIout.PE_%i,2);',idx64));
    eval(sprintf('ASWakePermEn(:,5,ifile) = nanmean(waketmp.wSMIout.PE_%i,2);',idx128));
    
    eval(sprintf('ASSleepwSMI(:,:,1,ifile) = nanmean(sleeptmp.wSMIout.wSMI_%i,3);',idx8));    
    eval(sprintf('ASSleepwSMI(:,:,2,ifile) = nanmean(sleeptmp.wSMIout.wSMI_%i,3);',idx16));    
    eval(sprintf('ASSleepwSMI(:,:,3,ifile) = nanmean(sleeptmp.wSMIout.wSMI_%i,3);',idx32));
    eval(sprintf('ASSleepwSMI(:,:,4,ifile) = nanmean(sleeptmp.wSMIout.wSMI_%i,3);',idx64));
    eval(sprintf('ASSleepwSMI(:,:,5,ifile) = nanmean(sleeptmp.wSMIout.wSMI_%i,3);',idx128));

    eval(sprintf('ASSleepPermEn(:,1,ifile) = nanmean(sleeptmp.wSMIout.PE_%i,2);',idx8));
    eval(sprintf('ASSleepPermEn(:,2,ifile) = nanmean(sleeptmp.wSMIout.PE_%i,2);',idx16));
    eval(sprintf('ASSleepPermEn(:,3,ifile) = nanmean(sleeptmp.wSMIout.PE_%i,2);',idx32));
    eval(sprintf('ASSleepPermEn(:,4,ifile) = nanmean(sleeptmp.wSMIout.PE_%i,2);',idx64));
    eval(sprintf('ASSleepPermEn(:,5,ifile) = nanmean(sleeptmp.wSMIout.PE_%i,2);',idx128));
    
    clear waketmp sleeptmp
end


%% integrate power within bands

%%% GET RELATIVE POWER %%%

ASWakeTot = squeeze(trapz(foi,ASWakePow,2));
ASWakeTotRS= repmat(reshape(ASWakeTot,19,1,size(ASWakePow,3)),1,size(ASWakePow,2),1); % reshape to repeat columns for frequency bins
ASWakePowR = ASWakePow./ASWakeTotRS; % relative power normalization

ASSleepTot = squeeze(trapz(foi,ASSleepPow,2));
ASSleepTotRS= repmat(reshape(ASSleepTot,19,1,size(ASSleepPow,3)),1,size(ASSleepPow,2),1); % reshape to repeat columns for frequency bins
ASSleepPowR = ASSleepPow./ASSleepTotRS; % relative power normalization

TDWakeTot = squeeze(trapz(foi,TDWakePow,2));
TDWakeTotRS= repmat(reshape(TDWakeTot,19,1,size(TDWakePow,3)),1,size(TDWakePow,2),1); % reshape to repeat columns for frequency bins
TDWakePowR = TDWakePow./TDWakeTotRS; % relative power normalization

TDSleepTot = squeeze(trapz(foi,TDSleepPow,2));
TDSleepTotRS= repmat(reshape(TDSleepTot,19,1,size(TDSleepPow,3)),1,size(TDSleepPow,2),1); % reshape to repeat columns for frequency bins
TDSleepPowR = TDSleepPow./TDSleepTotRS; % relative power normalization


DSWakeTot = squeeze(trapz(foi,DSWakePow,2));
DSWakeTotRS= repmat(reshape(DSWakeTot,19,1,size(DSWakePow,3)),1,size(DSWakePow,2),1); % reshape to repeat columns for frequency bins
DSWakePowR = DSWakePow./DSWakeTotRS; % relative power normalization

DSSleepTot = squeeze(trapz(foi,DSSleepPow,2));
DSSleepTotRS= repmat(reshape(DSSleepTot,19,1,size(DSSleepPow,3)),1,size(DSSleepPow,2),1); % reshape to repeat columns for frequency bins
DSSleepPowR = DSSleepPow./DSSleepTotRS; % relative power normalization



%% Do averaging

ASWakePowCA = squeeze(log10(nanmean(ASWakePow))); % average across channels BEFORE log-scaling
ASSleepPowCA = squeeze(log10(nanmean(ASSleepPow))); % average across channels BEFORE log-scaling
DSWakePowCA = squeeze(log10(nanmean(DSWakePow))); % average across channels BEFORE log-scaling
DSSleepPowCA = squeeze(log10(nanmean(DSSleepPow))); % average across channels BEFORE log-scaling
TDWakePowCA = squeeze(log10(nanmean(TDWakePow))); % average across channels BEFORE log-scaling
TDSleepPowCA = squeeze(log10(nanmean(TDSleepPow))); % average across channels BEFORE log-scaling

ASWakePowRCA = squeeze(log10(nanmean(ASWakePowR))); % average across channels BEFORE log-scaling
ASSleepPowRCA = squeeze(log10(nanmean(ASSleepPowR))); % average across channels BEFORE log-scaling
DSWakePowRCA = squeeze(log10(nanmean(DSWakePowR))); % average across channels BEFORE log-scaling
DSSleepPowRCA = squeeze(log10(nanmean(DSSleepPowR))); % average across channels BEFORE log-scaling
TDWakePowRCA = squeeze(log10(nanmean(TDWakePowR))); % average across channels BEFORE log-scaling
TDSleepPowRCA = squeeze(log10(nanmean(TDSleepPowR))); % average across channels BEFORE log-scaling


ASPowDiff = ASSleepPowCA-ASWakePowCA;
DSPowDiff = DSSleepPowCA-DSWakePowCA;
TDPowDiff = TDSleepPowCA-TDWakePowCA;

ASPowRDiff = ASSleepPowRCA-ASWakePowRCA;
DSPowRDiff = DSSleepPowRCA-DSWakePowRCA;
TDPowRDiff = TDSleepPowRCA-TDWakePowRCA;

foi_hd = 2.^[linspace(-1,5,200)];

% Absolute power: For AS, first average within subjects
[~,ASWakePowCAnm] = nestedmean(ASWakePowCA,ASsids,2);
[~,ASSleepPowCAnm] = nestedmean(ASSleepPowCA,ASsids,2);
[~,ASPowDiffnm] = nestedmean(ASPowDiff,ASsids,2);

[~,~,ASWakeci] = ttest(ASWakePowCAnm');
[~,~,ASSleepci] = ttest(ASSleepPowCAnm');
[~,~,ASDiffci] = ttest(ASPowDiffnm');

[~,~,DSWakeci] = ttest(DSWakePowCA');
[~,~,DSSleepci] = ttest(DSSleepPowCA');
[~,~,DSDiffci] = ttest(DSPowDiff');

[~,~,TDWakeci] = ttest(TDWakePowCA');
[~,~,TDSleepci] = ttest(TDSleepPowCA');
[~,~,TDDiffci] = ttest(TDPowDiff');

% Relative power: For AS, first average within subjects
[~,ASWakePowRCAnm] = nestedmean(ASWakePowRCA,ASsids,2);
[~,ASSleepPowRCAnm] = nestedmean(ASSleepPowRCA,ASsids,2);
[~,ASPowRDiffnm] = nestedmean(ASPowRDiff,ASsids,2);

[~,~,ASWakeRci] = ttest(ASWakePowRCAnm');
[~,~,ASSleepRci] = ttest(ASSleepPowRCAnm');
[~,~,ASDiffRci] = ttest(ASPowRDiffnm');

[~,~,DSWakeRci] = ttest(DSWakePowRCA');
[~,~,DSSleepRci] = ttest(DSSleepPowRCA');
[~,~,DSDiffRci] = ttest(DSPowRDiff');

[~,~,TDWakeRci] = ttest(TDWakePowRCA');
[~,~,TDSleepRci] = ttest(TDSleepPowRCA');
[~,~,TDDiffRci] = ttest(TDPowRDiff');


ASWakeLZcCA = mean(ASWakeLZc);
ASsleepLZcCA = mean(ASSleepLZc);
TDWakeLZcCA = mean(TDWakeLZc);
TDsleepLZcCA = mean(TDSleepLZc);
DSWakeLZcCA = mean(DSWakeLZc);
DSSleepLZcCA = mean(DSSleepLZc);

ASWakeCTWCA = mean(ASWakeCTW);
ASsleepCTWCA = mean(ASSleepCTW);
TDWakeCTWCA = mean(TDWakeCTW);
TDsleepCTWCA = mean(TDSleepCTW);
DSWakeCTWCA = mean(DSWakeCTW);
DSSleepCTWCA = mean(DSSleepCTW);


ASWakeMSECA = mean(ASWakeMSE);
ASsleepMSECA = mean(ASSleepMSE);
TDWakeMSECA = mean(TDWakeMSE);
TDsleepMSECA = mean(TDSleepMSE);
DSWakeMSECA = mean(DSWakeMSE);
DSSleepMSECA = mean(DSSleepMSE);


%% Plot power

switch plotting
    case true
        
        wkclr = [1 0 0];
        spclr = [0 0 1];
        dsclr = [1 0 0.5];
        tdclr = [0 0.75 0.5];
        asclr = [0 0.5 1];
        dsclr2 = [1 0 0.5]./1.5;
        tdclr2 = [0 0.75 0.5]./1.5;
        asclr2 = [0 0.5 1]./1.5;

        myfigure
        plot(log2(foi_hd),interp1(foi,nanmean(DSWakePowCA,2)',foi_hd,'spline'),'Color',wkclr,'LineWidth',2)
        hold on
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,DSWakeci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,DSWakeci(2,:),foi_hd,'spline'))],wkclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(DSSleepPowCA,2)',foi_hd,'spline'),'Color',spclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,DSSleepci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,DSSleepci(2,:),foi_hd,'spline'))],spclr,'facealpha',0.3,'EdgeColor','none')
        
        legend({'Dup15q Wake','Dup15q Wake 95% CI','Dup15q Sleep','Dup15q Sleep 95% CI'},'fontsize',16','location',...
            'northeastoutside')
        legend boxoff
        ylim([-1 4])
        xlim([-1 5])
        xlabel('Frequency (Hz)','FontSize',12)
        ylabel('Power log_1_0(\muV^{2}/Hz)','FontSize',18)
        title(sprintf('Dup15q absolute Power, N = %i',sum(~isnan(DSWakePowCA(1,:)))),'FontSize',20)
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        box off
        legend boxoff
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', 22)
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', 22)
        set(gca, 'TickDir', 'out');
        makefighandsome
        print('-dpng','./figures/Dup15qAbsPower.png')
        print('-dsvg','./figures/Dup15qAbsPower.svg')
        
        
        myfigure
        plot(log2(foi_hd),interp1(foi,nanmean(TDWakePowCA,2)',foi_hd,'spline'),'Color',wkclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,TDWakeci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,TDWakeci(2,:),foi_hd,'spline'))],wkclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(TDSleepPowCA,2)',foi_hd,'spline'),'Color',spclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,TDSleepci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,TDSleepci(2,:),foi_hd,'spline'))],spclr,'facealpha',0.3,'EdgeColor','none')
        
        legend({'TD Wake','TD Wake 95% CI','TD Sleep','TD Sleep 95% CI'},'fontsize',16','location',...
            'northeastoutside')
        legend boxoff
        ylim([-1 4])
        xlim([-1 5])
        xlabel('Frequency (Hz)','FontSize',12)
        ylabel('Power log_1_0(\muV^{2}/Hz)','FontSize',18)
        title(sprintf('TD absolute Power, N = %i',sum(~isnan(TDWakePowCA(1,:)))),'FontSize',20)
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        box off
        legend boxoff
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', 22)
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', 22)
        set(gca, 'TickDir', 'out');
        makefighandsome
        print('-dpng','./figures/TDAbsPower.png')
        print('-dsvg','./figures/TDAbsPower.svg')
        
        
        myfigure
        plot(log2(foi_hd),interp1(foi,nestedmean(ASWakePowCA,ASsids,2)',foi_hd,'spline'),'Color',wkclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ASWakeci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,ASWakeci(2,:),foi_hd,'spline'))],wkclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nestedmean(ASSleepPowCA,ASsids,2)',foi_hd,'spline'),'Color',spclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ASSleepci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,ASSleepci(2,:),foi_hd,'spline'))],spclr,'facealpha',0.3,'EdgeColor','none')
        
        legend({'AS Wake','AS Wake 95% CI','AS Sleep','AS Sleep 95% CI'},'fontsize',16','location',...
            'northeastoutside')
        legend boxoff
        ylim([-1 4])
        xlim([-1 5])
        xlabel('Frequency (Hz)','FontSize',12)
        ylabel('Power log_1_0(\muV^{2}/Hz)','FontSize',18)
        title(sprintf('AS absolute Power, N = %i',sum(~isnan(ASWakePowCAnm(1,:)))),'FontSize',20)
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        box off
        legend boxoff
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', 22)
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', 22)
        set(gca, 'TickDir', 'out');
        makefighandsome
        print('-dpng','./figures/ASAbsPower.png')
        print('-dsvg','./figures/ASAbsPower.svg')
        
        %% Export info for figure 3 
        
        Freq = foi';
        Tfig3 = table(Freq,ASWakePowCA,ASSleepPowCA,TDWakePowCA,TDSleepPowCA,DSWakePowCA,DSSleepPowCA);
        writetable(Tfig3,'./Table_fig3.csv')
        
        %% Invidual PSD traces

        myfigure
        plot(log2(foi_hd),interp1(foi,DSWakePowCA,foi_hd),'k')
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        ylim([-1 4])
        xlabel('Frequency (Hz)','FontSize',12)
        ylabel('Power log_{10}(\muV^{2}/Hz)')
        title('Dup15q wake')
        makefigpretty
        print('-dsvg','./Figures/Dup15q_wake_all_traces.svg')

        myfigure
        plot(log2(foi_hd),interp1(foi,DSSleepPowCA,foi_hd),'k')
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        ylim([-1 4])
        xlabel('Frequency (Hz)','FontSize',12)
        ylabel('Power log_{10}(\muV^{2}/Hz)')
        title('Dup15q sleep')
        makefigpretty
        print('-dsvg','./Figures/Dup15q_sleep_all_traces.svg')

                myfigure
        plot(log2(foi_hd),interp1(foi,ASWakePowCA,foi_hd),'k')
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        ylim([-1 4])
        xlabel('Frequency (Hz)','FontSize',12)
        ylabel('Power log_{10}(\muV^{2}/Hz)')
        title('AS wake')
        makefigpretty
        print('-dsvg','./Figures/AS_wake_all_traces.svg')

                myfigure
        plot(log2(foi_hd),interp1(foi,ASSleepPowCA,foi_hd),'k')
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        ylim([-1 5])
        xlabel('Frequency (Hz)','FontSize',12)
        ylabel('Power log_{10}(\muV^{2}/Hz)')
        title('AS sleep')
        makefigpretty
        print('-dsvg','./Figures/AS_Sleep_all_traces.svg')

                myfigure
        plot(log2(foi_hd),interp1(foi,TDWakePowCA,foi_hd),'k')
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        ylim([-1 4])
        xlabel('Frequency (Hz)','FontSize',12)
        ylabel('Power log_{10}(\muV^{2}/Hz)')
        title('TD wake')
        makefigpretty
        print('-dsvg','./Figures/TD_wake_all_traces.svg')

                myfigure
        plot(log2(foi_hd),interp1(foi,TDSleepPowCA,foi_hd),'k')
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        ylim([-1 4])
        xlabel('Frequency (Hz)','FontSize',12)
        ylabel('Power log_{10}(\muV^{2}/Hz)')
        title('TD sleep')
        makefigpretty
        print('-dsvg','./Figures/TD_sleep_all_traces.svg')

        save dataPeaks 
        %keyboard
        
        %% Power differences, absolute power
        
        myfigure
        
        plot(log2(foi_hd),interp1(foi,nanmean(ASPowDiff,2)',foi_hd,'spline'),'Color',asclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ASDiffci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,ASDiffci(2,:),foi_hd,'spline'))],asclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(TDPowDiff,2)',foi_hd,'spline'),'Color',tdclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,TDDiffci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,TDDiffci(2,:),foi_hd,'spline'))],tdclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(DSPowDiff,2)',foi_hd,'spline'),'Color',dsclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,DSDiffci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,DSDiffci(2,:),foi_hd,'spline'))],dsclr,'facealpha',0.3,'EdgeColor','none')
        
        
        xlim([-1 5])
        ylim([-1.5 1.5])
        
        xlabel('Frequency (Hz)','FontSize',12)
        ylabel('Sleep - wake power [log_1_0(\muV^{2}/Hz)]','FontSize',18)
        title('Sleep-wake absolute Power','FontSize',20)
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        box off
        legend boxoff
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', 22)
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', 22)
        set(gca, 'TickDir', 'out');
        makefighandsome
        
        % add highlights to show frequency bands
        fill([log2(0.5) log2(1) log2(1) log2(0.5)],[-3 -3 3 3],'k','facealpha',0.1,'edgealpha',0)
        fill([log2(1) log2(2) log2(2) log2(1)],[-3 -3 3 3],'b','facealpha',0.1,'edgealpha',0)
        fill([log2(2) log2(4) log2(4) log2(2)],[-3 -3 3 3],'r','facealpha',0.1,'edgealpha',0)
        fill([log2(4) log2(8) log2(8) log2(4)],[-3 -3 3 3],'y','facealpha',0.1,'edgealpha',0)
        fill([log2(8) log2(16) log2(16) log2(8)],[-3 -3 3 3],'c','facealpha',0.1,'edgealpha',0)
        fill([log2(16) log2(32) log2(32) log2(16)],[-3 -3 3 3],'g','facealpha',0.1,'edgealpha',0)
        
        legend({'AS Sleep-Wake','AS Sleep-Wake 95% CI','TD Sleep-Wake','TD Sleep-Wake 95% CI',...
            'Dup15q Sleep-Wake','Dup15q Sleep-Wake 95% CI','s','\delta1',...
            '\delta2','\theta','\alpha-\sigma','\beta'},'fontsize',16','location',...
            'northeastoutside','autoupdate','off')
        legend boxoff
        axis square
        
        plot(log2(foi_hd),zeros(1,length(foi_hd)),'k:')
        print('-dpng','./figures/AbsPowerWakeMinusSleep.png')
        print('-dsvg','./figures/AbsPowerWakeMinusSleep.svg')
        
        %% Power differences, relative power
        
        myfigure
        
        plot(log2(foi_hd),interp1(foi,nanmean(ASPowRDiff,2)',foi_hd,'spline'),'Color',asclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ASDiffRci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,ASDiffRci(2,:),foi_hd,'spline'))],asclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(TDPowRDiff,2)',foi_hd,'spline'),'Color',tdclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,TDDiffRci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,TDDiffRci(2,:),foi_hd,'spline'))],tdclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(DSPowRDiff,2)',foi_hd,'spline'),'Color',dsclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,DSDiffRci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,DSDiffRci(2,:),foi_hd,'spline'))],dsclr,'facealpha',0.3,'EdgeColor','none')
        
        
        xlim([-1 5])
        ylim([-1.5 1.5])
        
        xlabel('Frequency (Hz)','FontSize',12)
        ylabel('Sleep - wake power [log_1_0(1/Hz)]','FontSize',18)
        title('Sleep-wake relative Power','FontSize',20)
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        box off
        legend boxoff
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', 22)
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', 22)
        set(gca, 'TickDir', 'out');
        makefighandsome
        
        % add highlights to show frequency bands
        fill([log2(0.5) log2(1) log2(1) log2(0.5)],[-3 -3 3 3],'k','facealpha',0.1,'edgealpha',0)
        fill([log2(1) log2(2) log2(2) log2(1)],[-3 -3 3 3],'b','facealpha',0.1,'edgealpha',0)
        fill([log2(2) log2(4) log2(4) log2(2)],[-3 -3 3 3],'r','facealpha',0.1,'edgealpha',0)
        fill([log2(4) log2(8) log2(8) log2(4)],[-3 -3 3 3],'y','facealpha',0.1,'edgealpha',0)
        fill([log2(8) log2(16) log2(16) log2(8)],[-3 -3 3 3],'c','facealpha',0.1,'edgealpha',0)
        fill([log2(16) log2(32) log2(32) log2(16)],[-3 -3 3 3],'g','facealpha',0.1,'edgealpha',0)
        
        legend({'AS Sleep-Wake','AS Sleep-Wake 95% CI','TD Sleep-Wake','TD Sleep-Wake 95% CI',...
            'Dup15q Sleep-Wake','Dup15q Sleep-Wake 95% CI','s','\delta1',...
            '\delta2','\theta','\alpha-\sigma','\beta'},'fontsize',16','location',...
            'northeastoutside','autoupdate','off')
        legend boxoff
        axis square
        
        plot(log2(foi_hd),zeros(1,length(foi_hd)),'k:')
        print('-dpng','./figures/RelPowerWakeMinusSleep.png')
        print('-dsvg','./figures/RelPowerWakeMinusSleep.svg')
        %% Absolute power, all cohorts
        myfigure
        
        plot(log2(foi_hd),interp1(foi,nanmean(ASWakePowCA,2)',foi_hd,'spline'),'Color',asclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ASWakeci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,ASWakeci(2,:),foi_hd,'spline'))],asclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(ASSleepPowCA,2)',foi_hd,'spline'),'Color',asclr2,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ASSleepci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,ASSleepci(2,:),foi_hd,'spline'))],asclr2,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(DSWakePowCA,2)',foi_hd,'spline'),'Color',dsclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,DSWakeci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,DSWakeci(2,:),foi_hd,'spline'))],dsclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(DSSleepPowCA,2)',foi_hd,'spline'),'Color',dsclr2,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,DSSleepci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,DSSleepci(2,:),foi_hd,'spline'))],dsclr2,'facealpha',0.3,'EdgeColor','none')
        
        
        plot(log2(foi_hd),interp1(foi,nanmean(TDWakePowCA,2)',foi_hd,'spline'),'Color',tdclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,TDWakeci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,TDWakeci(2,:),foi_hd,'spline'))],tdclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(TDSleepPowCA,2)',foi_hd,'spline'),'Color',tdclr2,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,TDSleepci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,TDSleepci(2,:),foi_hd,'spline'))],tdclr2,'facealpha',0.3,'EdgeColor','none')
        
        
        xlim([-1 5])
        ylim([-1 4])
        
        xlabel('Frequency (Hz)','FontSize',12)
        ylabel('Power [log_1_0(\muV^{2}/Hz)]','FontSize',18)
        title('Absolute Power','FontSize',20)
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        box off
        legend boxoff
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', 22)
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', 22)
        set(gca, 'TickDir', 'out');
        makefighandsome
        
        % add highlights to show frequency bands
        fill([log2(0.5) log2(1) log2(1) log2(0.5)],[-3 -3 4 4],'m','facealpha',0.1,'edgealpha',0)
        fill([log2(1) log2(2) log2(2) log2(1)],[-3 -3 4 4],'b','facealpha',0.1,'edgealpha',0)
        fill([log2(2) log2(4) log2(4) log2(2)],[-3 -3 4 4],'r','facealpha',0.1,'edgealpha',0)
        fill([log2(4) log2(8) log2(8) log2(4)],[-3 -3 4 4],'y','facealpha',0.1,'edgealpha',0)
        fill([log2(8) log2(16) log2(16) log2(8)],[-3 -3 4 4],'c','facealpha',0.1,'edgealpha',0)
        fill([log2(16) log2(32) log2(32) log2(16)],[-3 -3 4 4],'g','facealpha',0.1,'edgealpha',0)
        
        legend({'AS Wake','AS Wake 95% CI',...
            'AS Sleep','AS Sleep 95% CI',...
            'DS Wake','DS Wake 95% CI',...
            'DS Sleep','DS Sleep 95% CI',...
            'TD Wake','TD Wake 95% CI',...
            'TD Sleep','TD Sleep 95% CI',...
            's','\delta1','\delta2','\theta',...
            '\alpha-\sigma','\beta'},'fontsize',16','location',...
            'northeastoutside','autoupdate','off')
        legend boxoff
        
        %plot(log2(foi_hd),zeros(1,length(foi_hd)),'k:')
        print('-dpng','./figures/AbsPowerWakeSleep.png')
        print('-dsvg','./figures/AbsPowerWakeSleep.svg')
        
        %% Relative power, all cohorts
        myfigure
        
        plot(log2(foi_hd),interp1(foi,nanmean(ASWakePowRCA,2)',foi_hd,'spline'),'Color',asclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ASWakeRci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,ASWakeRci(2,:),foi_hd,'spline'))],asclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(ASSleepPowRCA,2)',foi_hd,'spline'),'Color',asclr2,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ASSleepRci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,ASSleepRci(2,:),foi_hd,'spline'))],asclr2,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(DSWakePowRCA,2)',foi_hd,'spline'),'Color',dsclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,DSWakeRci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,DSWakeRci(2,:),foi_hd,'spline'))],dsclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(DSSleepPowRCA,2)',foi_hd,'spline'),'Color',dsclr2,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,DSSleepRci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,DSSleepRci(2,:),foi_hd,'spline'))],dsclr2,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(TDWakePowRCA,2)',foi_hd,'spline'),'Color',tdclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,TDWakeRci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,TDWakeRci(2,:),foi_hd,'spline'))],tdclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(TDSleepPowRCA,2)',foi_hd,'spline'),'Color',tdclr2,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,TDSleepRci(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,TDSleepRci(2,:),foi_hd,'spline'))],tdclr2,'facealpha',0.3,'EdgeColor','none')
        
        
        xlim([-1 5])
        ylim([-4.5 0.5])
        
        xlabel('Frequency (Hz)','FontSize',12)
        ylabel('Power [log_1_0(1/Hz)]','FontSize',18)
        title('Relative Power','FontSize',20)
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        box off
        legend boxoff
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', 22)
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', 22)
        set(gca, 'TickDir', 'out');
        makefighandsome
        
        % add highlights to show frequency bands
        fill([log2(0.5) log2(1) log2(1) log2(0.5)],[-5 -5 1 1 ],'m','facealpha',0.1,'edgealpha',0)
        fill([log2(1) log2(2) log2(2) log2(1)],[-5 -5 1 1 ],'b','facealpha',0.1,'edgealpha',0)
        fill([log2(2) log2(4) log2(4) log2(2)],[-5 -5 1 1 ],'r','facealpha',0.1,'edgealpha',0)
        fill([log2(4) log2(8) log2(8) log2(4)],[-5 -5 1 1 ],'y','facealpha',0.1,'edgealpha',0)
        fill([log2(8) log2(16) log2(16) log2(8)],[-5 -5 1 1 ],'c','facealpha',0.1,'edgealpha',0)
        fill([log2(16) log2(32) log2(32) log2(16)],[-5 -5 1 1 ],'g','facealpha',0.1,'edgealpha',0)
        
        legend({'AS Wake','AS Wake 95% CI',...
            'AS Sleep','AS Sleep 95% CI',...
            'DS Wake','DS Wake 95% CI',...
            'DS Sleep','DS Sleep 95% CI',...
            'TD Wake','TD Wake 95% CI',...
            'TD Sleep','TD Sleep 95% CI',...
            's','\delta1','\delta2','\theta',...
            '\alpha-\sigma','\beta'},'fontsize',16','location',...
            'northeastoutside','autoupdate','off')
        legend boxoff
        
        %plot(log2(foi_hd),zeros(1,length(foi_hd)),'k:')
        print('-dpng','./figures/RelPowerWakeSleep.png')
        print('-dsvg','./figures/RelPowerWakeSleep.svg')
        
        %% Export info for figure 3 (relative power
        
        Freq = foi';
        Tfig3_rel = table(Freq,ASWakePowRCA,ASSleepPowRCA,TDWakePowRCA,TDSleepPowRCA,DSWakePowRCA,DSSleepPowRCA);
        writetable(Tfig3_rel,'./Table_fig3_rel.csv')
             
        
        %% dwPLI CIs
        [~,~,ASWakeciSR] = ttest(ASWakedwPLI_SR);
        [~,~,ASSleepciSR] = ttest(ASSleepdwPLI_SR);

        [~,~,DSWakeciSR] = ttest(DSWakedwPLI_SR);
        [~,~,DSSleepciSR] = ttest(DSSleepdwPLI_SR);

        [~,~,TDWakeciSR] = ttest(TDWakedwPLI_SR);
        [~,~,TDSleepciSR] = ttest(TDSleepdwPLI_SR);
        
        [~,~,ASWakeciLR] = ttest(ASWakedwPLI_LR);
        [~,~,ASSleepciLR] = ttest(ASSleepdwPLI_LR);

        [~,~,DSWakeciLR] = ttest(DSWakedwPLI_LR);
        [~,~,DSSleepciLR] = ttest(DSSleepdwPLI_LR);

        [~,~,TDWakeciLR] = ttest(TDWakedwPLI_LR);
        [~,~,TDSleepciLR] = ttest(TDSleepdwPLI_LR);
        
        %% Plot shortrange dwPLI spectrum
        myfigure
        
        plot(log2(foi_hd),interp1(foi,nanmean(ASWakedwPLI_SR)',foi_hd,'spline'),'Color',asclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ASWakeciSR(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,ASWakeciSR(2,:),foi_hd,'spline'))],asclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(ASSleepdwPLI_SR)',foi_hd,'spline'),'Color',asclr2,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ASSleepciSR(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,ASSleepciSR(2,:),foi_hd,'spline'))],asclr2,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(DSWakedwPLI_SR)',foi_hd,'spline'),'Color',dsclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,DSWakeciSR(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,DSWakeciSR(2,:),foi_hd,'spline'))],dsclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(DSSleepdwPLI_SR)',foi_hd,'spline'),'Color',dsclr2,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,DSSleepciSR(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,DSSleepciSR(2,:),foi_hd,'spline'))],dsclr2,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(TDWakedwPLI_SR)',foi_hd,'spline'),'Color',tdclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,TDWakeciSR(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,TDWakeciSR(2,:),foi_hd,'spline'))],tdclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(TDSleepdwPLI_SR)',foi_hd,'spline'),'Color',tdclr2,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,TDSleepciSR(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,TDSleepciSR(2,:),foi_hd,'spline'))],tdclr2,'facealpha',0.3,'EdgeColor','none')
        
        
        xlim([-1 5])
        ylim([0 0.4])
        
        xlabel('Frequency (Hz)','FontSize',12)
        ylabel('dwPLI','FontSize',18)
        title('Shortrange spectral connectivity','FontSize',20)
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        box off
        legend boxoff
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', 22)
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', 22)
        set(gca, 'TickDir', 'out');
        makefighandsome
        
        % add highlights to show frequency bands
        fill([log2(0.5) log2(1) log2(1) log2(0.5)],[-5 -5 1 1 ],'m','facealpha',0.1,'edgealpha',0)
        fill([log2(1) log2(2) log2(2) log2(1)],[-5 -5 1 1 ],'b','facealpha',0.1,'edgealpha',0)
        fill([log2(2) log2(4) log2(4) log2(2)],[-5 -5 1 1 ],'r','facealpha',0.1,'edgealpha',0)
        fill([log2(4) log2(8) log2(8) log2(4)],[-5 -5 1 1 ],'y','facealpha',0.1,'edgealpha',0)
        fill([log2(8) log2(16) log2(16) log2(8)],[-5 -5 1 1 ],'c','facealpha',0.1,'edgealpha',0)
        fill([log2(16) log2(32) log2(32) log2(16)],[-5 -5 1 1 ],'g','facealpha',0.1,'edgealpha',0)
        
        legend({'AS Wake','AS Wake 95% CI',...
            'AS Sleep','AS Sleep 95% CI',...
            'DS Wake','DS Wake 95% CI',...
            'DS Sleep','DS Sleep 95% CI',...
            'TD Wake','TD Wake 95% CI',...
            'TD Sleep','TD Sleep 95% CI',...
            's','\delta1','\delta2','\theta',...
            '\alpha-\sigma','\beta'},'fontsize',16','location',...
            'northeastoutside','autoupdate','off')
        legend boxoff
        
        print('-dpng','./figures/SRdwPLIWakeSleep.png')
        print('-dsvg','./figures/SRdwPLIWakeSleep.svg')
        
        %% Export info for figure 4 (short range)
        
        Freq = foi';
        Tfig4_dwPLI_SR = table(Freq,ASWakedwPLI_SR',ASSleepdwPLI_SR',TDWakedwPLI_SR',TDSleepdwPLI_SR',DSWakedwPLI_SR',DSSleepdwPLI_SR');
        writetable(Tfig4_dwPLI_SR,'./Table_fig4_wPLI_SR.csv')
        
        
        
        %% Plot longrange dwPLI spectrum
        myfigure
        
        plot(log2(foi_hd),interp1(foi,nanmean(ASWakedwPLI_LR)',foi_hd,'spline'),'Color',asclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ASWakeciLR(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,ASWakeciLR(2,:),foi_hd,'spline'))],asclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(ASSleepdwPLI_LR)',foi_hd,'spline'),'Color',asclr2,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,ASSleepciLR(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,ASSleepciLR(2,:),foi_hd,'spline'))],asclr2,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(DSWakedwPLI_LR)',foi_hd,'spline'),'Color',dsclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,DSWakeciLR(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,DSWakeciLR(2,:),foi_hd,'spline'))],dsclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(DSSleepdwPLI_LR)',foi_hd,'spline'),'Color',dsclr2,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,DSSleepciLR(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,DSSleepciLR(2,:),foi_hd,'spline'))],dsclr2,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(TDWakedwPLI_LR)',foi_hd,'spline'),'Color',tdclr,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,TDWakeciLR(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,TDWakeciLR(2,:),foi_hd,'spline'))],tdclr,'facealpha',0.3,'EdgeColor','none')
        
        plot(log2(foi_hd),interp1(foi,nanmean(TDSleepdwPLI_LR)',foi_hd,'spline'),'Color',tdclr2,'LineWidth',2)
        patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,TDSleepciLR(1,:),foi_hd,'spline') ...
            fliplr(interp1(foi,TDSleepciLR(2,:),foi_hd,'spline'))],tdclr2,'facealpha',0.3,'EdgeColor','none')
        
        
        xlim([-1 5])
        ylim([0 0.4])
        
        xlabel('Frequency (Hz)','FontSize',12)
        ylabel('dwPLI','FontSize',18)
        title('Longrange spectral connectivity','FontSize',20)
        set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
        box off
        legend boxoff
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', 22)
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', 22)
        set(gca, 'TickDir', 'out');
        makefighandsome
        
        % add highlights to show frequency bands
        fill([log2(0.5) log2(1) log2(1) log2(0.5)],[-5 -5 1 1 ],'m','facealpha',0.1,'edgealpha',0)
        fill([log2(1) log2(2) log2(2) log2(1)],[-5 -5 1 1 ],'b','facealpha',0.1,'edgealpha',0)
        fill([log2(2) log2(4) log2(4) log2(2)],[-5 -5 1 1 ],'r','facealpha',0.1,'edgealpha',0)
        fill([log2(4) log2(8) log2(8) log2(4)],[-5 -5 1 1 ],'y','facealpha',0.1,'edgealpha',0)
        fill([log2(8) log2(16) log2(16) log2(8)],[-5 -5 1 1 ],'c','facealpha',0.1,'edgealpha',0)
        fill([log2(16) log2(32) log2(32) log2(16)],[-5 -5 1 1 ],'g','facealpha',0.1,'edgealpha',0)
        
        legend({'AS Wake','AS Wake 95% CI',...
            'AS Sleep','AS Sleep 95% CI',...
            'DS Wake','DS Wake 95% CI',...
            'DS Sleep','DS Sleep 95% CI',...
            'TD Wake','TD Wake 95% CI',...
            'TD Sleep','TD Sleep 95% CI',...
            's','\delta1','\delta2','\theta',...
            '\alpha-\sigma','\beta'},'fontsize',16','location',...
            'northeastoutside','autoupdate','off')
        legend boxoff
        
        print('-dpng','./figures/LRdwPLIWakeSleep.png')
        print('-dsvg','./figures/LRdwPLIWakeSleep.svg')
        
        %% Export info for figure (long range)
        
        Freq = foi';
        Tfig4_dwPLI_LR = table(Freq,ASWakedwPLI_LR',ASSleepdwPLI_LR',TDWakedwPLI_LR',TDSleepdwPLI_LR',DSWakedwPLI_LR',DSSleepdwPLI_LR');
        writetable(Tfig4_dwPLI_LR,'./Table_fig4_wPLI_LR.csv')
        
        
    case false
        % do nothing
end

%% Extract slow power (0.5 - 1.0 Hz)

f1=0.5; f2 = 1;
idf1 = foi>= f1 & foi<f2;

slow_wk_AS = squeeze(trapz(foi(idf1),ASWakePow(:,idf1,:),2));
slow_sp_AS = squeeze(trapz(foi(idf1),ASSleepPow(:,idf1,:),2));
slow_rel_wk_AS = squeeze(trapz(foi(idf1),ASWakePowR(:,idf1,:),2));
slow_rel_sp_AS = squeeze(trapz(foi(idf1),ASSleepPowR(:,idf1,:),2));

slow_wk_TD = squeeze(trapz(foi(idf1),TDWakePow(:,idf1,:),2));
slow_sp_TD = squeeze(trapz(foi(idf1),TDSleepPow(:,idf1,:),2));
slow_rel_wk_TD = squeeze(trapz(foi(idf1),TDWakePowR(:,idf1,:),2));
slow_rel_sp_TD = squeeze(trapz(foi(idf1),TDSleepPowR(:,idf1,:),2));

slow_wk_DS = squeeze(trapz(foi(idf1),DSWakePow(:,idf1,:),2));
slow_sp_DS = squeeze(trapz(foi(idf1),DSSleepPow(:,idf1,:),2));
slow_rel_wk_DS = squeeze(trapz(foi(idf1),DSWakePowR(:,idf1,:),2));
slow_rel_sp_DS = squeeze(trapz(foi(idf1),DSSleepPowR(:,idf1,:),2));

%% Extract delta1 power (1 - 2 Hz)

f1=1; f2 = 2;
idf2 = foi>= f1 & foi<f2;

delta1_wk_AS = squeeze(trapz(foi(idf2),ASWakePow(:,idf2,:),2));
delta1_sp_AS = squeeze(trapz(foi(idf2),ASSleepPow(:,idf2,:),2));
delta1_rel_wk_AS = squeeze(trapz(foi(idf2),ASWakePowR(:,idf2,:),2));
delta1_rel_sp_AS = squeeze(trapz(foi(idf2),ASSleepPowR(:,idf2,:),2));

delta1_wk_TD = squeeze(trapz(foi(idf2),TDWakePow(:,idf2,:),2));
delta1_sp_TD = squeeze(trapz(foi(idf2),TDSleepPow(:,idf2,:),2));
delta1_rel_wk_TD = squeeze(trapz(foi(idf2),TDWakePowR(:,idf2,:),2));
delta1_rel_sp_TD = squeeze(trapz(foi(idf2),TDSleepPowR(:,idf2,:),2));

delta1_wk_DS = squeeze(trapz(foi(idf2),DSWakePow(:,idf2,:),2));
delta1_sp_DS = squeeze(trapz(foi(idf2),DSSleepPow(:,idf2,:),2));
delta1_rel_wk_DS = squeeze(trapz(foi(idf2),DSWakePowR(:,idf2,:),2));
delta1_rel_sp_DS = squeeze(trapz(foi(idf2),DSSleepPowR(:,idf2,:),2));

%% Extract delta2 power (2 - 4 Hz)

f1=2; f2 = 4;
idf3 = foi>= f1 & foi<f2;

delta2_wk_AS = squeeze(trapz(foi(idf3),ASWakePow(:,idf3,:),2));
delta2_sp_AS = squeeze(trapz(foi(idf3),ASSleepPow(:,idf3,:),2));
delta2_rel_wk_AS = squeeze(trapz(foi(idf3),ASWakePowR(:,idf3,:),2));
delta2_rel_sp_AS = squeeze(trapz(foi(idf3),ASSleepPowR(:,idf3,:),2));

delta2_wk_TD = squeeze(trapz(foi(idf3),TDWakePow(:,idf3,:),2));
delta2_sp_TD = squeeze(trapz(foi(idf3),TDSleepPow(:,idf3,:),2));
delta2_rel_wk_TD = squeeze(trapz(foi(idf3),TDWakePowR(:,idf3,:),2));
delta2_rel_sp_TD = squeeze(trapz(foi(idf3),TDSleepPowR(:,idf3,:),2));

delta2_wk_DS = squeeze(trapz(foi(idf3),DSWakePow(:,idf3,:),2));
delta2_sp_DS = squeeze(trapz(foi(idf3),DSSleepPow(:,idf3,:),2));
delta2_rel_wk_DS = squeeze(trapz(foi(idf3),DSWakePowR(:,idf3,:),2));
delta2_rel_sp_DS = squeeze(trapz(foi(idf3),DSSleepPowR(:,idf3,:),2));

%% Extract theta power (4 - 8 Hz)

f1=4; f2 = 8;
idf4 = foi>= f1 & foi<f2;

theta_wk_AS = squeeze(trapz(foi(idf4),ASWakePow(:,idf4,:),2));
theta_sp_AS = squeeze(trapz(foi(idf4),ASSleepPow(:,idf4,:),2));
theta_rel_wk_AS = squeeze(trapz(foi(idf4),ASWakePowR(:,idf4,:),2));
theta_rel_sp_AS = squeeze(trapz(foi(idf4),ASSleepPowR(:,idf4,:),2));

theta_wk_TD = squeeze(trapz(foi(idf4),TDWakePow(:,idf4,:),2));
theta_sp_TD = squeeze(trapz(foi(idf4),TDSleepPow(:,idf4,:),2));
theta_rel_wk_TD = squeeze(trapz(foi(idf4),TDWakePowR(:,idf4,:),2));
theta_rel_sp_TD = squeeze(trapz(foi(idf4),TDSleepPowR(:,idf4,:),2));

theta_wk_DS = squeeze(trapz(foi(idf4),DSWakePow(:,idf4,:),2));
theta_sp_DS = squeeze(trapz(foi(idf4),DSSleepPow(:,idf4,:),2));
theta_rel_wk_DS = squeeze(trapz(foi(idf4),DSWakePowR(:,idf4,:),2));
theta_rel_sp_DS = squeeze(trapz(foi(idf4),DSSleepPowR(:,idf4,:),2));

%% Extract alpha-sigma power (8 - 16 Hz)

f1=8; f2 = 16;
idf5 = foi>= f1 & foi<f2;

alpha_wk_AS = squeeze(trapz(foi(idf5),ASWakePow(:,idf5,:),2));
alpha_sp_AS = squeeze(trapz(foi(idf5),ASSleepPow(:,idf5,:),2));
alpha_rel_wk_AS = squeeze(trapz(foi(idf5),ASWakePowR(:,idf5,:),2));
alpha_rel_sp_AS = squeeze(trapz(foi(idf5),ASSleepPowR(:,idf5,:),2));

alpha_wk_TD = squeeze(trapz(foi(idf5),TDWakePow(:,idf5,:),2));
alpha_sp_TD = squeeze(trapz(foi(idf5),TDSleepPow(:,idf5,:),2));
alpha_rel_wk_TD = squeeze(trapz(foi(idf5),TDWakePowR(:,idf5,:),2));
alpha_rel_sp_TD = squeeze(trapz(foi(idf5),TDSleepPowR(:,idf5,:),2));

alpha_wk_DS = squeeze(trapz(foi(idf5),DSWakePow(:,idf5,:),2));
alpha_sp_DS = squeeze(trapz(foi(idf5),DSSleepPow(:,idf5,:),2));
alpha_rel_wk_DS = squeeze(trapz(foi(idf5),DSWakePowR(:,idf5,:),2));
alpha_rel_sp_DS = squeeze(trapz(foi(idf5),DSSleepPowR(:,idf5,:),2));


%% Extract beta power (16 - 32 Hz)

f1=16; f2 = 32;
idf6 = foi>= f1 & foi<f2;

beta_wk_AS = squeeze(trapz(foi(idf6),ASWakePow(:,idf6,:),2));
beta_sp_AS = squeeze(trapz(foi(idf6),ASSleepPow(:,idf6,:),2));
beta_rel_wk_AS = squeeze(trapz(foi(idf6),ASWakePowR(:,idf6,:),2));
beta_rel_sp_AS = squeeze(trapz(foi(idf6),ASSleepPowR(:,idf6,:),2));

beta_wk_TD = squeeze(trapz(foi(idf6),TDWakePow(:,idf6,:),2));
beta_sp_TD = squeeze(trapz(foi(idf6),TDSleepPow(:,idf6,:),2));
beta_rel_wk_TD = squeeze(trapz(foi(idf6),TDWakePowR(:,idf6,:),2));
beta_rel_sp_TD = squeeze(trapz(foi(idf6),TDSleepPowR(:,idf6,:),2));

beta_wk_DS = squeeze(trapz(foi(idf6),DSWakePow(:,idf6,:),2));
beta_sp_DS = squeeze(trapz(foi(idf6),DSSleepPow(:,idf6,:),2));
beta_rel_wk_DS = squeeze(trapz(foi(idf6),DSWakePowR(:,idf6,:),2));
beta_rel_sp_DS = squeeze(trapz(foi(idf6),DSSleepPowR(:,idf6,:),2));


%% wSMI

%%% normalize wSMI by permutation entropy as 2*[wSMI(X,Y)/(PE(X)+PE(Y))]

for itau = 1:ntau % for each timescale
    
    for ifile = 1:size(ASWakewSMI,4)
        ASWakewSMI(:,:,itau,ifile) = ASWakewSMI(:,:,itau,ifile).*2 ./ (ASWakePermEn(:,itau,ifile) + ASWakePermEn(:,itau,ifile)');
        ASSleepwSMI(:,:,itau,ifile) = ASSleepwSMI(:,:,itau,ifile).*2 ./ (ASSleepPermEn(:,itau,ifile) + ASSleepPermEn(:,itau,ifile)');
    end
    
    for ifile = 1:size(DSWakewSMI,4)
        DSWakewSMI(:,:,itau,ifile) = DSWakewSMI(:,:,itau,ifile).*2 ./ (DSWakePermEn(:,itau,ifile) + DSWakePermEn(:,itau,ifile)');
        DSSleepwSMI(:,:,itau,ifile) = DSSleepwSMI(:,:,itau,ifile).*2 ./ (DSSleepPermEn(:,itau,ifile) + DSSleepPermEn(:,itau,ifile)');
    end
    
    for ifile = 1:size(TDWakewSMI,4)
        TDWakewSMI(:,:,itau,ifile) = TDWakewSMI(:,:,itau,ifile).*2 ./ (TDWakePermEn(:,itau,ifile) + TDWakePermEn(:,itau,ifile)');
        TDSleepwSMI(:,:,itau,ifile) = TDSleepwSMI(:,:,itau,ifile).*2 ./ (TDSleepPermEn(:,itau,ifile) + TDSleepPermEn(:,itau,ifile)');
    end
end

SR = reshape(shortrange,nchan^2,1);
LR = reshape(longrange,nchan^2,1);

% replace 0s with nans

ASWakewSMI(ASWakewSMI==0) = nan;
ASSleepwSMI(ASSleepwSMI==0) = nan;
TDWakewSMI(TDWakewSMI==0) = nan;
TDSleepwSMI(TDSleepwSMI==0) = nan;
DSWakewSMI(DSWakewSMI==0) = nan;
DSSleepwSMI(DSSleepwSMI==0) = nan;

ASWakewSMIrs = reshape(ASWakewSMI,nchan^2,ntau,size(ASWakewSMI,4));
ASSleepwSMIrs = reshape(ASSleepwSMI,nchan^2,ntau,size(ASSleepwSMI,4));
TDWakewSMIrs = reshape(TDWakewSMI,nchan^2,ntau,size(TDWakewSMI,4));
TDSleepwSMIrs = reshape(TDSleepwSMI,nchan^2,ntau,size(TDSleepwSMI,4));
DSWakewSMIrs = reshape(DSWakewSMI,nchan^2,ntau,size(DSWakewSMI,4));
DSSleepwSMIrs = reshape(DSSleepwSMI,nchan^2,ntau,size(DSSleepwSMI,4));

kernels = 2.^[3:7];

for itau = 1:ntau
    % Short range wSMI
    eval(sprintf('ASspTau%iSR = nanmean(squeeze(ASSleepwSMIrs(SR,%i,:)));',kernels(itau),itau));
    eval(sprintf('ASwkTau%iSR = nanmean(squeeze(ASWakewSMIrs(SR,%i,:)));',kernels(itau),itau));
    eval(sprintf('TDspTau%iSR = nanmean(squeeze(TDSleepwSMIrs(SR,%i,:)));',kernels(itau),itau));
    eval(sprintf('TDwkTau%iSR = nanmean(squeeze(TDWakewSMIrs(SR,%i,:)));',kernels(itau),itau));
    eval(sprintf('DSspTau%iSR = nanmean(squeeze(DSSleepwSMIrs(SR,%i,:)));',kernels(itau),itau));
    eval(sprintf('DSwkTau%iSR = nanmean(squeeze(DSWakewSMIrs(SR,%i,:)));',kernels(itau),itau)); 
    % Long range wSMI
    eval(sprintf('ASspTau%iLR = nanmean(squeeze(ASSleepwSMIrs(LR,%i,:)));',kernels(itau),itau));
    eval(sprintf('ASwkTau%iLR = nanmean(squeeze(ASWakewSMIrs(LR,%i,:)));',kernels(itau),itau));
    eval(sprintf('TDspTau%iLR = nanmean(squeeze(TDSleepwSMIrs(LR,%i,:)));',kernels(itau),itau));
    eval(sprintf('TDwkTau%iLR = nanmean(squeeze(TDWakewSMIrs(LR,%i,:)));',kernels(itau),itau));
    eval(sprintf('DSspTau%iLR = nanmean(squeeze(DSSleepwSMIrs(LR,%i,:)));',kernels(itau),itau));
    eval(sprintf('DSwkTau%iLR = nanmean(squeeze(DSWakewSMIrs(LR,%i,:)));',kernels(itau),itau));   
    % PermEn
    eval(sprintf('ASPermEnwkTau%i = squeeze(nanmean(ASWakePermEn(:,%i,:),1));',kernels(itau),itau));
    eval(sprintf('ASPermEnwkTau%i = squeeze(nanmean(ASWakePermEn(:,%i,:),1));',kernels(itau),itau));
    eval(sprintf('TDPermEnwkTau%i = squeeze(nanmean(TDWakePermEn(:,%i,:),1));',kernels(itau),itau));
    eval(sprintf('TDPermEnwkTau%i = squeeze(nanmean(TDWakePermEn(:,%i,:),1));',kernels(itau),itau));
    eval(sprintf('DSPermEnwkTau%i = squeeze(nanmean(DSWakePermEn(:,%i,:),1));',kernels(itau),itau));
    eval(sprintf('DSPermEnwkTau%i = squeeze(nanmean(DSWakePermEn(:,%i,:),1));',kernels(itau),itau));
    
    eval(sprintf('ASPermEnspTau%i = squeeze(nanmean(ASSleepPermEn(:,%i,:),1));',kernels(itau),itau));
    eval(sprintf('ASPermEnspTau%i = squeeze(nanmean(ASSleepPermEn(:,%i,:),1));',kernels(itau),itau));
    eval(sprintf('TDPermEnspTau%i = squeeze(nanmean(TDSleepPermEn(:,%i,:),1));',kernels(itau),itau));
    eval(sprintf('TDPermEnspTau%i = squeeze(nanmean(TDSleepPermEn(:,%i,:),1));',kernels(itau),itau));
    eval(sprintf('DSPermEnspTau%i = squeeze(nanmean(DSSleepPermEn(:,%i,:),1));',kernels(itau),itau));
    eval(sprintf('DSPermEnspTau%i = squeeze(nanmean(DSSleepPermEn(:,%i,:),1));',kernels(itau),itau));
end


%% dwPLI

ASWakedwPLIrs = reshape(ASWakedwPLI,nchan^2,6,size(ASWakedwPLI,4));
ASSleepdwPLIrs = reshape(ASSleepdwPLI,nchan^2,6,size(ASSleepdwPLI,4));
TDWakedwPLIrs = reshape(TDWakedwPLI,nchan^2,6,size(TDWakedwPLI,4));
TDSleepdwPLIrs = reshape(TDSleepdwPLI,nchan^2,6,size(TDSleepdwPLI,4));
DSWakedwPLIrs = reshape(DSWakedwPLI,nchan^2,6,size(DSWakedwPLI,4));
DSSleepdwPLIrs = reshape(DSSleepdwPLI,nchan^2,6,size(DSSleepdwPLI,4));

bands = {'Slow','Delta1','Delta2','Theta','Alpha','Beta'};

for ibnd = 1:length(bands)
    
    eval(sprintf('ASsp%sSR = nanmean(squeeze(ASSleepdwPLIrs(SR,ibnd,:)));',bands{ibnd}));
    eval(sprintf('ASwk%sSR = nanmean(squeeze(ASWakedwPLIrs(SR,ibnd,:)));',bands{ibnd}));
    eval(sprintf('DSsp%sSR = nanmean(squeeze(DSSleepdwPLIrs(SR,ibnd,:)));',bands{ibnd}));
    eval(sprintf('DSwk%sSR = nanmean(squeeze(DSWakedwPLIrs(SR,ibnd,:)));',bands{ibnd}));
    eval(sprintf('TDsp%sSR = nanmean(squeeze(TDSleepdwPLIrs(SR,ibnd,:)));',bands{ibnd}));
    eval(sprintf('TDwk%sSR = nanmean(squeeze(TDWakedwPLIrs(SR,ibnd,:)));',bands{ibnd}));
    
    eval(sprintf('ASsp%sLR = nanmean(squeeze(ASSleepdwPLIrs(LR,ibnd,:)));',bands{ibnd}));
    eval(sprintf('ASwk%sLR = nanmean(squeeze(ASWakedwPLIrs(LR,ibnd,:)));',bands{ibnd}));
    eval(sprintf('DSsp%sLR = nanmean(squeeze(DSSleepdwPLIrs(LR,ibnd,:)));',bands{ibnd}));
    eval(sprintf('DSwk%sLR = nanmean(squeeze(DSWakedwPLIrs(LR,ibnd,:)));',bands{ibnd}));
    eval(sprintf('TDsp%sLR = nanmean(squeeze(TDSleepdwPLIrs(LR,ibnd,:)));',bands{ibnd}));
    eval(sprintf('TDwk%sLR = nanmean(squeeze(TDWakedwPLIrs(LR,ibnd,:)));',bands{ibnd}));
end

% Uncomment below to use the wPLI without debiasing 

% %% wPLI
% 
% ASwakewPLIrs = reshape(ASwakewPLI,nchan^2,6,size(ASwakewPLI,4));
% ASsleepwPLIrs = reshape(ASsleepwPLI,nchan^2,6,size(ASsleepwPLI,4));
% TDwakewPLIrs = reshape(TDwakewPLI,nchan^2,6,size(TDwakewPLI,4));
% TDsleepwPLIrs = reshape(TDsleepwPLI,nchan^2,6,size(TDsleepwPLI,4));
% DSwakewPLIrs = reshape(DSwakewPLI,nchan^2,6,size(DSwakewPLI,4));
% DSsleepwPLIrs = reshape(DSsleepwPLI,nchan^2,6,size(DSsleepwPLI,4));
% 
% bands = {'Slow','Delta1','Delta2','Theta','Alpha','Beta'};
% 
% for ibnd = 1:length(bands)
%     
%     eval(sprintf('ASsp%sSR = nanmean(squeeze(ASsleepwPLIrs(SR,ibnd,:)));',bands{ibnd}));
%     eval(sprintf('ASwk%sSR = nanmean(squeeze(ASwakewPLIrs(SR,ibnd,:)));',bands{ibnd}));
%     eval(sprintf('DSsp%sSR = nanmean(squeeze(DSsleepwPLIrs(SR,ibnd,:)));',bands{ibnd}));
%     eval(sprintf('DSwk%sSR = nanmean(squeeze(DSwakewPLIrs(SR,ibnd,:)));',bands{ibnd}));
%     eval(sprintf('TDsp%sSR = nanmean(squeeze(TDsleepwPLIrs(SR,ibnd,:)));',bands{ibnd}));
%     eval(sprintf('TDwk%sSR = nanmean(squeeze(TDwakewPLIrs(SR,ibnd,:)));',bands{ibnd}));
%     
%     eval(sprintf('ASsp%sLR = nanmean(squeeze(ASsleepwPLIrs(LR,ibnd,:)));',bands{ibnd}));
%     eval(sprintf('ASwk%sLR = nanmean(squeeze(ASwakewPLIrs(LR,ibnd,:)));',bands{ibnd}));
%     eval(sprintf('DSsp%sLR = nanmean(squeeze(DSsleepwPLIrs(LR,ibnd,:)));',bands{ibnd}));
%     eval(sprintf('DSwk%sLR = nanmean(squeeze(DSwakewPLIrs(LR,ibnd,:)));',bands{ibnd}));
%     eval(sprintf('TDsp%sLR = nanmean(squeeze(TDsleepwPLIrs(LR,ibnd,:)));',bands{ibnd}));
%     eval(sprintf('TDwk%sLR = nanmean(squeeze(TDwakewPLIrs(LR,ibnd,:)));',bands{ibnd}));
% end



%%
nTD = length(TDwkTau32SR);
nDS = length(DSwkTau32SR);
nAS = length(ASwkTau32SR);

keepDup15q = true; % retain the dup15q group

varnames = {'Sleep','Group','Conscious','Subject','Sleep2','Group2','Conscious2','Age','mMSE','LZc','CTW',...
    'PermEn8','PermEn16','PermEn32','PermEn64',...
    'PermEn128','SRwSMI8','SRwSMI16','SRwSMI32','SRwSMI64','SRwSMI128',...
    'LRwSMI8','LRwSMI16','LRwSMI32','LRwSMI64','LRwSMI128','slow','delta1','delta2',...
    'theta','alpha','beta','slowR','delta1R','delta2R',...
    'thetaR','alphaR','betaR','SRdwPLIslow','SRdwPLIdelta1','SRdwPLIdelta2',...
    'SRdwPLItheta','SRdwPLIalpha','SRdwPLIbeta','LRdwPLIslow','LRdwPLIdelta1','LRdwPLIdelta2',...
    'LRdwPLItheta','LRdwPLIalpha','LRdwPLIbeta'};

vartype = {'IV','IV','IV','IV','IV','IV','IV','IV','scEntropy','scEntropy','scEntropy',...
    'scEntropy','scEntropy','scEntropy','scEntropy','scEntropy',...
    'fcEntropy','fcEntropy','fcEntropy','fcEntropy',...
    'fcEntropy','fcEntropy','fcEntropy','fcEntropy',...
    'fcEntropy','fcEntropy','scSpectralA','scSpectralA','scSpectralA',...
    'scSpectralA','scSpectralA','scSpectralA','scSpectralR','scSpectralR',...
    'scSpectralR','scSpectralR','scSpectralR','scSpectralR','fcSpectral',...
    'fcSpectral','fcSpectral','fcSpectral','fcSpectral','fcSpectral',...
    'fcSpectral','fcSpectral','fcSpectral','fcSpectral','fcSpectral','fcSpectral'};

% Identify categorical variables
iscat = contains(varnames,'Sleep') | contains(varnames,'Group') | contains(varnames,'Subject') | contains(varnames,'Conscious');

% number of EEG features, disregard first 5 variables which aren't EEG
ntest = length(varnames)-5;

%Create categorical variables
issleep = categorical([repmat({'No'},nAS+nTD+nDS,1); repmat({'Yes'},nAS+nTD+nDS,1)]);
isconscious = categorical([repmat({'Yes'},nAS+nTD+nDS,1); repmat({'No'},nAS+nTD+nDS,1)]);
group = categorical([repmat({'AS'},nAS,1); repmat({'TD'},nTD,1); repmat({'DS'},nDS,1);...
    repmat({'AS'},nAS,1); repmat({'TD'},nTD,1); repmat({'DS'},nDS,1)]);
subjects = categorical([ASsids'; TDsids'; DSsids'; ASsids'; TDsids'; DSsids']);

% also add non-categorical versions of variables (compatible with mes
% toolbox)

sleepBin = [zeros(nAS+nTD+nDS,1); ones(nAS+nTD+nDS,1)];
groupBin = [zeros(nAS,1); ones(nTD,1); repmat(2,nDS,1);...
    zeros(nAS,1); ones(nTD,1); repmat(2,nDS,1)];
consciousBin = [ones(nAS+nTD+nDS,1); zeros(nAS+nTD+nDS,1)];


assert(~isordinal(issleep),'Categorical var is treated as ordinal')
assert(~isordinal(group),'Categorical var is treated as ordinal')
assert(~isordinal(subjects),'Categorical var is treated as ordinal')

% Now make the HUGE table :)
T = table(issleep,group,isconscious,subjects,sleepBin,groupBin,consciousBin,...
    [ASages'; TDages'; DSages'; ASages'; TDages'; DSages'],...
    [ASWakeMSECA'; TDWakeMSECA'; DSWakeMSECA'; ASsleepMSECA'; TDsleepMSECA'; DSSleepMSECA'],...
    [ASWakeLZcCA'; TDWakeLZcCA'; DSWakeLZcCA'; ASsleepLZcCA'; TDsleepLZcCA'; DSSleepLZcCA'],...
    [ASWakeCTWCA'; TDWakeCTWCA'; DSWakeCTWCA'; ASsleepCTWCA'; TDsleepCTWCA'; DSSleepCTWCA'],...
    [ASPermEnwkTau8; TDPermEnwkTau8; DSPermEnwkTau8; ASPermEnspTau8; TDPermEnspTau8; DSPermEnspTau8],...
    [ASPermEnwkTau16; TDPermEnwkTau16; DSPermEnwkTau16; ASPermEnspTau16; TDPermEnspTau16; DSPermEnspTau16],...
    [ASPermEnwkTau32; TDPermEnwkTau32; DSPermEnwkTau32; ASPermEnspTau32; TDPermEnspTau32; DSPermEnspTau32],...
    [ASPermEnwkTau64; TDPermEnwkTau64; DSPermEnwkTau64; ASPermEnspTau64; TDPermEnspTau64; DSPermEnspTau64],...
    [ASPermEnwkTau128; TDPermEnwkTau128; DSPermEnwkTau128; ASPermEnspTau128; TDPermEnspTau128; DSPermEnspTau128],...
    [ASwkTau8SR'; TDwkTau8SR'; DSwkTau8SR'; ASspTau8SR'; TDspTau8SR'; DSspTau8SR'],...
    [ASwkTau16SR'; TDwkTau16SR'; DSwkTau16SR'; ASspTau16SR'; TDspTau16SR'; DSspTau16SR'],...
    [ASwkTau32SR'; TDwkTau32SR'; DSwkTau32SR'; ASspTau32SR'; TDspTau32SR'; DSspTau32SR'],...
    [ASwkTau64SR'; TDwkTau64SR'; DSwkTau64SR'; ASspTau64SR'; TDspTau64SR'; DSspTau64SR'],...
    [ASwkTau128SR'; TDwkTau128SR'; DSwkTau128SR'; ASspTau128SR'; TDspTau128SR'; DSspTau128SR'],...
    [ASwkTau8LR'; TDwkTau8LR'; DSwkTau8LR'; ASspTau8LR'; TDspTau8LR'; DSspTau8LR'],...
    [ASwkTau16LR'; TDwkTau16LR'; DSwkTau16LR'; ASspTau16LR'; TDspTau16LR'; DSspTau16LR'],...    
    [ASwkTau32LR'; TDwkTau32LR'; DSwkTau32LR'; ASspTau32LR'; TDspTau32LR'; DSspTau32LR'],...
    [ASwkTau64LR'; TDwkTau64LR'; DSwkTau64LR'; ASspTau64LR'; TDspTau64LR'; DSspTau64LR'],...
    [ASwkTau128LR'; TDwkTau128LR'; DSwkTau128LR'; ASspTau128LR'; TDspTau128LR'; DSspTau128LR'],...
    [mean(slow_wk_AS)'; mean(slow_wk_TD)'; mean(slow_wk_DS)'; mean(slow_wk_AS)'; mean(slow_sp_TD)'; mean(slow_sp_DS)'],...
    [mean(delta1_wk_AS)'; mean(delta1_wk_TD)'; mean(delta1_wk_DS)'; mean(delta1_sp_AS)'; mean(delta1_sp_TD)'; mean(delta1_sp_DS)'],...
    [mean(delta2_wk_AS)'; mean(delta2_wk_TD)'; mean(delta2_wk_DS)'; mean(delta2_sp_AS)'; mean(delta2_sp_TD)'; mean(delta2_sp_DS)'],...
    [mean(theta_wk_AS)'; mean(theta_wk_TD)'; mean(theta_wk_DS)'; mean(theta_sp_AS)'; mean(theta_sp_TD)'; mean(theta_sp_DS)'],...
    [mean(alpha_wk_AS)'; mean(alpha_wk_TD)'; mean(alpha_wk_DS)'; mean(alpha_sp_AS)'; mean(alpha_sp_TD)'; mean(alpha_sp_DS)'],...
    [mean(beta_wk_AS)'; mean(beta_wk_TD)'; mean(beta_wk_DS)'; mean(beta_sp_AS)'; mean(beta_sp_TD)'; mean(beta_sp_DS)'],...
    [mean(slow_rel_wk_AS)'; mean(slow_rel_wk_TD)'; mean(slow_rel_wk_DS)'; mean(slow_rel_sp_AS)';mean(slow_rel_sp_TD)'; mean(slow_rel_sp_DS)'],...
    [mean(delta1_rel_wk_AS)'; mean(delta1_rel_wk_TD)'; mean(delta1_rel_wk_DS)'; mean(delta1_rel_sp_AS)';mean(delta1_rel_sp_TD)'; mean(delta1_rel_sp_DS)'],...
    [mean(delta2_rel_wk_AS)'; mean(delta2_rel_wk_TD)'; mean(delta2_rel_wk_DS)'; mean(delta2_rel_sp_AS)';mean(delta2_rel_sp_TD)'; mean(delta2_rel_sp_DS)'],...
    [mean(theta_rel_wk_AS)'; mean(theta_rel_wk_TD)'; mean(theta_rel_wk_DS)'; mean(theta_rel_sp_AS)';mean(theta_rel_sp_TD)'; mean(theta_rel_sp_DS)'],...
    [mean(alpha_rel_wk_AS)'; mean(alpha_rel_wk_TD)'; mean(alpha_rel_wk_DS)'; mean(alpha_rel_sp_AS)';mean(alpha_rel_sp_TD)'; mean(alpha_rel_sp_DS)'],...
    [mean(beta_rel_wk_AS)'; mean(beta_rel_wk_TD)'; mean(beta_rel_wk_DS)'; mean(beta_rel_sp_AS)';mean(beta_rel_sp_TD)'; mean(beta_rel_sp_DS)'],...
    [ASwkSlowSR'; TDwkSlowSR'; DSwkSlowSR'; ASspSlowSR'; TDspSlowSR'; DSspSlowSR'],...
    [ASwkDelta1SR'; TDwkDelta1SR'; DSwkDelta1SR'; ASspDelta1SR'; TDspDelta1SR'; DSspDelta1SR'],...
    [ASwkDelta2SR'; TDwkDelta2SR'; DSwkDelta2SR'; ASspDelta2SR'; TDspDelta2SR'; DSspDelta2SR'],...
    [ASwkThetaSR'; TDwkThetaSR'; DSwkThetaSR'; ASspThetaSR'; TDspThetaSR'; DSspThetaSR'],...
    [ASwkAlphaSR'; TDwkAlphaSR'; DSwkAlphaSR'; ASspAlphaSR'; TDspAlphaSR'; DSspAlphaSR'],...
    [ASwkBetaSR'; TDwkBetaSR'; DSwkBetaSR'; ASspBetaSR'; TDspBetaSR'; DSspBetaSR'],...
    [ASwkSlowLR'; TDwkSlowLR'; DSwkSlowLR'; ASspSlowLR'; TDspSlowLR'; DSspSlowLR'],...
    [ASwkDelta1LR'; TDwkDelta1LR'; DSwkDelta1LR'; ASspDelta1LR'; TDspDelta1LR'; DSspDelta1LR'],...
    [ASwkDelta2LR'; TDwkDelta2LR'; DSwkDelta2LR'; ASspDelta2LR'; TDspDelta2LR'; DSspDelta2LR'],...
    [ASwkThetaLR'; TDwkThetaLR'; DSwkThetaLR'; ASspThetaLR'; TDspThetaLR'; DSspThetaLR'],...
    [ASwkAlphaLR'; TDwkAlphaLR'; DSwkAlphaLR'; ASspAlphaLR'; TDspAlphaLR'; DSspAlphaLR'],...
    [ASwkBetaLR'; TDwkBetaLR'; DSwkBetaLR'; ASspBetaLR'; TDspBetaLR'; DSspBetaLR'],...
    'VariableNames',varnames);

writetable(T,'./EEGfeatures.csv')


switch keepDup15q % do we want to keep Dup15q in analysis?
    case true
        fprintf('Keeping all three groups (AS, TD, Dup15q)\n')
    case false
        
        dsidx = find(all(char(group)=='DS',2)); % identify Dup15q subjects
        T2 = T; % copy table
        
        T(dsidx,:) = []; % remove Dup15q subjects
        % We must redefine group variable after removing Dup15q subjects
        T.Group = categorical([repmat({'AS'},nAS,1); repmat({'TD'},nTD,1);...
            repmat({'AS'},nAS,1); repmat({'TD'},nTD,1)]);
        
end
%%
% Dup15q participant sex
% 5003 - Male
% 5004 - Male
% 5007 - Female
% 5012 - Female
% 5013 - Female
% 5014 - Male
% 5046 - Male
% 5077 - Female
% 5078 - Female
% 5080 - Male
% 5088 - Male

DSsex = [0 0 1 1 1 0 0 1 1 0 0 ];

TT = table([AS15qcn';TD15qcn';DS15qcn';AS15qcn';TD15qcn';DS15qcn'],...
    [ASsex';TDsex';DSsex';ASsex';TDsex';DSsex'],...
    [ASWakeDatLen'; TDWakeDatLen'; DSWakeDatLen'; ASSleepDatLen'; TDSleepDatLen'; DSSleepDatLen'],...
    'VariableNames',{'CopyNumber','IsFemale','Usable_data_length'})
    
Tsup = [T(:,2:3) TT(:,1:2) T(:,8) TT(:,3) T(:,9:end)] % supplementa table


%% remove missing data from table
missing = isnan(T.CTW);

% make sure that if we remove the entry from wake, we also remove the sleep
% entry for the same subjects, and vice versa

missing1 = missing(1:length(missing)/2);
missing2 = missing(length(missing)/2+1:end);

missing = missing | [missing2; missing1];

T(missing,:) = [];

TDidx = T.Group == categorical({'TD'});
DSidx = T.Group == categorical({'DS'});
ASidx = T.Group == categorical({'AS'});

Tsup(missing,:) = [];

writetable(Tsup,'./SupplementalData.csv')


%% Z-score each predictor

T2 = T; % copy table so we still have the original values
frstvar = find(strcmp(vartype,'IV'),1,'last')+1;
for icol = frstvar:size(T,2)
    T{:,icol} = zscore(T{:,icol});
end


%% Report demographics

fprintf('%i participants with AS\n',sum(ASidx)/2)
fprintf('%i participants with TD\n',sum(TDidx)/2)
fprintf('%i participants with DS\n',sum(DSidx)/2)

AS_idx_wake = ASWakeDatLen ~= 0;
AS_idx_sleep = ASSleepDatLen ~= 0;
TD_idx_wake = TDWakeDatLen ~= 0;
TD_idx_sleep = TDSleepDatLen ~= 0;
DS_idx_wake = DSWakeDatLen ~= 0;
DS_idx_sleep = DSSleepDatLen ~= 0;

assert(all(AS_idx_wake == AS_idx_sleep),'Different AS datasets rejected for sleep and wake')
assert(all(TD_idx_wake == TD_idx_sleep),'Different TD datasets rejected for sleep and wake')
assert(all(DS_idx_wake == DS_idx_sleep),'Different DS datasets rejected for sleep and wake')

ASWakeDatLen = ASWakeDatLen(AS_idx_wake);
ASSleepDatLen = ASSleepDatLen(AS_idx_sleep);
TDWakeDatLen = TDWakeDatLen(TD_idx_wake);
TDSleepDatLen = TDSleepDatLen(TD_idx_sleep);
DSWakeDatLen = DSWakeDatLen(DS_idx_wake);
DSSleepDatLen = DSSleepDatLen(DS_idx_sleep);

ASages = ASages(AS_idx_wake);
ASsids = ASsids(AS_idx_wake);
TDages = TDages(TD_idx_wake);
TDsids = TDsids(TD_idx_wake);
DSages = DSages(DS_idx_wake);
DSsids = DSsids(DS_idx_wake);

fprintf('AS wake data = %2.1f +/- %2.1f minutes (mean +/- std)\n',...
    mean(ASWakeDatLen),std(ASWakeDatLen))
fprintf('AS sleep data = %2.1f +/- %2.1f minutes (mean +/- std)\n',...
    mean(ASSleepDatLen),std(ASSleepDatLen))
fprintf('AS wake data = %2.1f +/- %2.1f minutes (nested mean +/- nested std)\n',...
    nestedmean(ASWakeDatLen',ASsids'),nestedstd(ASWakeDatLen',ASsids'))
fprintf('AS sleep data = %2.1f +/- %2.1f minutes (nested mean +/- nested std)\n',...
    nestedmean(ASSleepDatLen',ASsids'),nestedstd(ASSleepDatLen',ASsids'))
fprintf('TD wake data = %2.1f +/- %2.1f minutes (mean +/- std)\n',...
    mean(TDWakeDatLen),std(TDWakeDatLen))
fprintf('TD sleep data = %2.1f +/- %2.1f minutes ( mean +/-  std)\n',...
    mean(TDSleepDatLen),std(TDSleepDatLen))
fprintf('DS wake data = %2.1f +/- %2.1f minutes ( mean +/-  std)\n',...
    mean(DSWakeDatLen),std(DSWakeDatLen))
fprintf('DS sleep data = %2.1f +/- %2.1f minutes ( mean +/-  std)\n',...
    mean(DSSleepDatLen),std(DSSleepDatLen))

fprintf('AS ages = %2.1f +/- %2.1f (mean +/- std)\n',mean(ASages), std(ASages))
fprintf('AS ages = %2.1f +/- %2.1f (nested mean +/- nested std)\n',nestedmean(ASages',ASsids'), nestedstd(ASages',ASsids'))
fprintf('TD ages = %2.1f +/- %2.1f (mean +/- std)\n',mean(TDages),std(TDages))
fprintf('DS ages = %2.1f +/- %2.1f (mean +/- std)\n',mean(DSages),std(DSages))

fprintf('%i AS datasets with deletion genotype and %i AS datasets with non-deletion genotype\n',...
    sum(AS15qcn==1),sum(AS15qcn==2))

[a,b] = unique(T.Subject(ASidx));
fprintf('%i AS subjects with deletion genotype and %i AS subjects with non-deletion genotype\n',...
    sum(AS15qcn(b)==1),sum(AS15qcn(b)==2))

fprintf('%i DS subjects with interstitial dups and %i DS subjects with isodicentric duplications\n\n',...
    sum(DS15qcn==3),sum(DS15qcn==4))

fprintf('AS sex ratio = %i:%i, M:F\n',sum(~ASsex),sum(ASsex))
fprintf('AS sex ratio (unique subjects only) = %i:%i, M:F\n',sum(~ASsex(b)),sum(ASsex(b)))
fprintf('TD sex ratio = %i:%i, M:F\n',sum(~TDsex),sum(TDsex))
fprintf('DS sex ratio = %i:%i, M:F\n',6,5) % hardcode

tmp = T.Subject(ASidx&T.Sleep2,:);
unitmp = unique(tmp);
howmany = nan(1,length(unitmp));
for isub = 1:length(unitmp)
    howmany(isub) = sum(tmp==unitmp(isub));
end

for i = 1:max(howmany)
    fprintf('     %i AS subjects with %i datasets\n',sum(howmany==i),i)
end


%% Do PCA on features

switch job
    case 'FindFeatures'
        
        % Add PCs to table based on wake - sleep AS data
        
        feats = unique(vartype);
        EEGvars = find(~contains(vartype,'IV'));
        Tpc = table(); % new table
        PC_coefs = cell(2,length(feats)-1);
        
        % Wake minus sleep
        wms = T{~logical(T.Sleep2),EEGvars} - T{logical(T.Sleep2),EEGvars};
        Td = table();
        for icol = 1:size(wms,2)
            eval(sprintf('Td.%s = wms(:,icol);',varnames{EEGvars(icol)}));
        end
        
        
        for ift = 2:length(feats)
            % do PCA just on AS data
            [COEFF, ~, LATENT] = pca(Td{ASidx(1:length(ASidx)/2),contains(vartype(EEGvars),feats{ift})});
            VE = LATENT./sum(LATENT).*100;
            PC_coefs{1,ift-1} = feats{ift};
            PC_coefs{2,ift-1} = COEFF;
            PC_coefs{3,ift-1} = VE;
            score = (T{:,contains(vartype,feats{ift})}-...
                nanmean(T{:,contains(vartype,feats{ift})}))*COEFF;
            eval(sprintf('Tpc.PC1_%s_VE%i = score(:,1);',replace(feats{ift},' ','_'),round(VE(1))));
            eval(sprintf('Tpc.PC2_%s_VE%i = score(:,2);',replace(feats{ift},' ','_'),round(VE(2))));
            eval(sprintf('Tpc.PC3_%s_VE%i = score(:,3);',replace(feats{ift},' ','_'),round(VE(3))));
            eval(sprintf('Tpc.PC4_%s_VE%i = score(:,4);',replace(feats{ift},' ','_'),round(VE(3))));
        end
        save PCtable Tpc PC_coefs
    case 'MachineLearning'
        load PCtable PC_coefs
        for ift = 1:size(PC_coefs,2)
            cVar = 0; % cumulative variance explained
            pccnt = 0;
            while cVar < varthresh % explain at least 75% of the variance
                pccnt = pccnt + 1;
                score = (T{:,contains(vartype,PC_coefs{1,ift})} - nanmean(T{:,contains(vartype,PC_coefs{1,ift})}))*PC_coefs{2,ift};
                eval(sprintf('T.PC%i_%s = score(:,%i);',pccnt,replace(PC_coefs{1,ift},' ','_'),pccnt));
                VaEx = PC_coefs{3,ift};
                cVar = cVar + VaEx(pccnt);
                % print out which variables contribute most to this PC
                v = varnames(contains(vartype,PC_coefs{1,ift}));
                contp = v(PC_coefs{2,ift}(:,1)>=0.5);
                contn = v(PC_coefs{2,ift}(:,1)<=0.5);
                for ict = 1:length(contp)
                    if ict == 1
                        fprintf('     Largest positive contributions to %s PC%i from %s',PC_coefs{1,ift},pccnt,contp{ict})
                    else
                        fprintf(' and %s',contp{ict})
                    end
                end
                fprintf('\n')
                for ict = 1:length(contn)
                    if ict == 1
                        fprintf('     Largest negative contributions to %s PC%i from %s',PC_coefs{1,ift},pccnt,contn{ict})
                    else
                        fprintf(' and %s',contn{ict})
                    end
                end
                fprintf('\n')
            end  
            fprintf('%i total PCs used for %s, Var. explained = %2.2f\n',pccnt,PC_coefs{1,ift},cVar)
        end
end

 %%

switch job
    case 'FindFeatures'
        
        %% Table and stats with only Angelman data
        
        Tas = T(T.Group == categorical({'AS'}),:);
        
        ASsleepF = nan(1,length(varnames));
        ASsleepES = nan(3,length(varnames));
        ASsleepP = nan(1,length(varnames));
        
        predictors = '1 + Sleep + (1|Subject)';
        exclude = 'Sleep Sleep2 Group Group2 Conscious Conscious2 Subject Age CopyNumber';
        
        for ivar = 1:length(varnames)
            if ~contains(exclude,varnames(ivar)) % if the response isn't already a predictor
                mdlstr = sprintf('%s ~  %s',varnames{ivar},predictors)
                lme = fitlme(Tas,mdlstr)
                lme_anova = anova(lme); % anova test on fixed effect terms
                lme_anova.FStat
                
                sIDX = find(strcmpi(lme.Coefficients.Name,'Sleep_Yes'));
                
                
                ASsleepF(ivar) = lme_anova.FStat(2);
                ASsleepP(ivar) = lme_anova.pValue(2);
                
                % Beta coef as effect size
                ASsleepES(1,ivar) = lme.Coefficients.Estimate(2);
                ASsleepES(2,ivar) = lme.Coefficients.Lower(2);
                ASsleepES(3,ivar) = lme.Coefficients.Upper(2);
            end
        end
        %%
        
        switch plotting
            case true
                
                % Main effect of sleep p-values
                
                myfigure2
                stem(-log10(ASsleepP),'filled','linewidth',2)
                xticks(1:length(varnames))
                xticklabels(varnames)
                xlim([8.5 length(varnames)+0.5])
                xtickangle(45)
                xlabel('EEG features')
                ylabel('AS sleep -log_{10}(p-value)')
                title('AS Biomarkers of conscious state','fontsize',18)
                
                box off
                set(gca,'linewidth',3)
                xAX = get(gca,'XAxis');
                set(xAX,'FontSize', 20)
                set(xAX,'color','k')
                yAX = get(gca,'YAxis');
                set(yAX,'FontSize', 20)
                set(yAX,'color','k')
                set(gca, 'TickDir', 'out')
                set(gcf,'color','w')
                set(gca,'Layer','top')
                axis normal
                print('-dpng','./Figures/SleepLMMMainEffectsPvalue.png')
                print('-dsvg','./Figures/SleepLMMMainEffectsPvalue.svg')
                
                %% Main effect of sleep F-value
                
                myfigure2
                stem(log10(ASsleepF),'filled','linewidth',2)
                xticks(1:length(varnames))
                xticklabels(varnames)
                xlim([8.5 length(varnames)+0.5])
                xtickangle(45)
                xlabel('EEG features')
                ylabel('AS sleep log_{10}(F-value)')
                title('AS Biomarkers of conscious state','fontsize',18)
                
                box off
                set(gca,'linewidth',3)
                xAX = get(gca,'XAxis');
                set(xAX,'FontSize', 20)
                set(xAX,'color','k')
                yAX = get(gca,'YAxis');
                set(yAX,'FontSize', 20)
                set(yAX,'color','k')
                set(gca, 'TickDir', 'out')
                set(gcf,'color','w')
                set(gca,'Layer','top')
                axis normal
                print('-dpng','./Figures/SleepLMMMainEffectsFvalue.png')
                print('-dsvg','./Figures/SleepLMMMainEffectsFvalue.svg')
                
                %% Main effect of sleep beta coefs
                
                myfigure2
                stem(ASsleepES(1,:),'filled','linewidth',2)
                pltidx = find(abs(ASsleepES(1,:))>=0.5);
                scatter(pltidx,ASsleepES(1,pltidx),130,	'k','filled','d')
                scatter(1:length(ASsleepES),ASsleepES(2,:),100,'r','filled','s')
                scatter(1:length(ASsleepES),ASsleepES(3,:),100,'r','filled','s')
                makefighandsome
                plot(1:100,ones(1,100).*0.5,'k:')
                plot(1:100,ones(1,100).*-0.5,'k:')
                xticks(1:length(varnames))
                xticklabels(varnames)
                              xlim([8.5 length(varnames)+0.5])
                ylim([-2 2])
                xtickangle(45)
                xlabel('EEG features')
                ylabel('AS sleep beta effect sizes')
                title('AS Biomarkers of conscious state','fontsize',18)
                
                print('-dpng','./Figures/SleepLMMMainEffectsBetas.png')
                print('-dsvg','./Figures/SleepLMMMainEffectsBetas.svg')
                
                Tbetas = table(varnames',ASsleepES(1,:)','VariableNames',{'Var','beta'})
                writetable(Tbetas,'Tbetas.csv')
                
            case false
                %Do nothing
        end
        
        save AS_sleep_wake_stats ASsleepP ASsleepF ASsleepES
        
    case 'MachineLearning'
        
        
        %% Regularized logistic regression
        load AS_sleep_wake_stats
        
        % Notes from JEFF
        % Do CV for regularization ONLY using training data, and when you do your ...
        % split for CV, make sure that both samples from the same subjects are in ...
        % either the test or the validation set. Don't spread between both.
               
        rng(58980423)
        Nboot = 1e4; % bootstrapped resamples
        % results is group x [AUC, accuracy, PPV, RECALL, SPEC] x feature type
        % x [real value, CI lower bound, CI upper bound]
        results = nan(3,5,5,3);
        ccnt = 0; % counter for classifications
                
        % ground truth (is this person conscious? frame is this way so we are talking about a consciousness detector)
        TDgt = T.Conscious(TDidx) == categorical({'Yes'});
        DSgt = T.Conscious(DSidx) == categorical({'Yes'});
        ASgt = T.Conscious(ASidx) == categorical({'Yes'});
               
        % Bootstrapping allocation
        feats = unique(vartype);
        
        % group x resample x feature x measure (AUC, accuracy, sensitivity,
        % specificity) 
                 
        % Use n for each bootstrapped resample,
        ASboot = sum(ASidx);
        TDboot = sum(TDidx);
        DSboot = sum(DSidx);
        
        ASidx2 = find(ASidx);
        TDidx2 = find(TDidx);
        DSidx2 = find(DSidx);

        %%% TRAIN ON ANGELMAN %%%
        fprintf('ANGELMAN TRAINING\n')
        
        % ROC curve preallocation
        xAS = cell(5,1);
        xTD = cell(5,1);
        xDS = cell(5,1);
        yAS = cell(5,1);
        yTD = cell(5,1);
        yDS = cell(5,1);
        
        % Choose regularization parameter using 10-fold cross validation on
        % training data
        kfold = 10; % number of CV folds
        S = T.Subject(ASidx); % unique subjects

        us = unique(S);
        cvid = 1 + mod((1:length(us))',kfold);
        indices = randperm(length(us)); % randperm 
        cvid = cvid(indices);

        mytest = nan(sum(ASidx),kfold);
        mytraining = nan(sum(ASidx),kfold);
        for icv = 1:kfold
            sidx = cvid == icv; % subjects belonging to test data for this CV
            mytest(:,icv) = ismember(S,us(sidx));
            mytraining(:,icv) = ~mytest(:,icv); % the ones not in the test data are training data
        end
        assert(all(sum(mytest,2) == ones(size(mytest,1),1)),'Each dataset must be used only once')
        assert(all(mytest + mytraining,[1 2]),'Some data are used both as test and training data')

        pcv = nan(length(feats)-1,1);
        u1cv = nan(length(feats)-1,1);
        u2cv = nan(length(feats)-1,1);
        AUC1cv = nan(length(feats)-1,3);
        AUC2cv = nan(length(feats)-1,3);
        
        for ift = 2:length(feats)
     
            switch method
                case 'LMM'
                    feat  = feats{ift};
                    fprintf('\n%s FEATURES ONLY using %s method\n',upper(feat),method)
                    
                    fts = find(ismember(vartype,feat));
                    assert(~isempty(fts),'No features found belonging to this category')
                    ranks = find(abs(ASsleepES(1,ismember(vartype,feat))) >= betathresh); % use features with abs(beta) >= 0.5
                    
                    tmp = T(:,fts(ranks));
                    predALL = T{:,fts(ranks)};
                    predAS = T{ASidx,fts(ranks)};
                    predTD = T{TDidx,fts(ranks)};
                    predDS = T{DSidx,fts(ranks)};
                    fprintf('     Using these features:')
                    for jft = 1:length(ranks)
                        fprintf(' %s ',tmp.Properties.VariableNames{jft})
                    end
                    fprintf('\n')
                case 'PCA'
                    feat  = feats{ift};
                    fprintf('\n%s FEATURES ONLY using %s method\n',upper(feat),method)
                    
                    fts = find(ismember(vartype,feat));
                    assert(~isempty(fts),'No features found belonging to this category')
                    feat  = replace(feats{ift},' ','_');
                    pcidx = contains(T.Properties.VariableNames,'PC') & contains(T.Properties.VariableNames,feat);
                    
                    predAS = T{ASidx,pcidx};
                    predTD = T{TDidx,pcidx};
                    predDS = T{DSidx,pcidx};
            end
                       
            [B0,FitInfo,CV] = mylassoglm(predAS,T.Conscious(ASidx),'binomial','Link','logit','CV',kfold,'Balanced',cat(3,mytraining,mytest));
            
            % Check that data from the same subjects are on same side of the CV partion
            IS = nan(1,kfold);
            for ifold = 1:kfold
                St = S(logical(CV.training(:,ifold))); % training subjects
                Sv = S(~logical(CV.training(:,ifold))); % test subjects
                IS(ifold) = length(intersect(St,Sv));
            end
            
            assert(max(IS)==0,'At least one fold had data from same subject on different sides of CV partition')
            
            % Examine the cross-validation plot to see the effect of the Lambda regularization parameter.
            lassoPlot(B0,FitInfo,'plottype','CV');
            legend('show') % Show legend
            
            [~,useme] = min(FitInfo.Deviance); % find optimal regularization
            lambda = FitInfo.Lambda(useme); % regularization parameter
            LAMBDA(ift-1) = lambda;
            
            % Do logistic regression just on AS data using the chosen lambda value
            [B,LassoFit] = lassoglm(predAS,T.Conscious(ASidx),'binomial','Link','logit','Lambda',lambda);
            
            coef = [LassoFit.Intercept; B];
            
            % apply model fit with optimal regularization to validation sets
            ASfit = glmval(coef,predAS,'logit','Constant','on');
            TDfit = glmval(coef,predTD,'logit','Constant','on');
            DSfit = glmval(coef,predDS,'logit','Constant','on');
            
            % Compare AS training performance with 10-fold cross validation
            k2fold = kfold; % use same number of CV folds as for the hyperparameter fitting
            
            cvid2 = 1 + mod((1:length(us))',k2fold);
            indices = randperm(length(us)); % randperm
            cvid2 = cvid2(indices);
            
            mytest2 = nan(sum(ASidx),k2fold);
            mytraining2 = nan(sum(ASidx),k2fold);
            for icv = 1:k2fold
                sidx = cvid2 == icv; % subjects belonging to test data for this CV
                mytest2(:,icv) = ismember(S,us(sidx));
                mytraining2(:,icv) = ~mytest2(:,icv); % the ones not in the test data are training data
            end
            
            assert(all(sum(mytest2,2) == ones(size(mytest2,1),1)),'Each dataset must be used only once')
            assert(all(mytest2 + mytraining2,[1 2]),'Some data are used both as test and training data')
            
            [B,LassoFit] = lassoglm(predAS,T.Conscious(ASidx),'binomial','Link','logit','Lambda',lambda);
            coef = [LassoFit.Intercept; B];
            ASfitcvx = [];
            ASgtcvx = [];
            theGT = logical(T.Conscious2(ASidx));
            for icv = 1:k2fold
                [Bx,LassoFitx] = lassoglm(predAS(logical(mytraining2(:,icv)),:),...
                    theGT(logical(mytraining2(:,icv))),'binomial','Link','logit','Lambda',lambda);
                tmpcoef = [LassoFitx.Intercept; Bx];
                ASfitcvx = [ASfitcvx; glmval(tmpcoef,predAS(logical(mytest2(:,icv)),:),'logit','Constant','on')];
                ASgtcvx = [ASgtcvx; ASgt(logical(mytest2(:,icv)))];
            end
            
            [~,~,~,AS_AUCcvx] = perfcurve(ASgtcvx,ASfitcvx,1,'XCrit','FPR',...
                'YCrit','TPR','Nboot',Nboot,'Alpha',0.05);
                        
            %%% ANGELMAN SYNDROME %%%
            % Test to see if the curve has the correct number of points, 
            % N + 1, where N is the number of data sests
            
            [X_AS,Y_AS,T_AS,AUC_AS,opt_AS,bs_AS] = myperfcurve(ASgt,ASfit,1,...
                'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Alpha',0.05);
            
            % Compare training AUC to k-fold CV AUC
            assert(size(ASgt,1)==size(ASgtcvx,1),'Not all subjects were tested with crossvalidation!?')
            assert(size(ASfit,1)==size(ASfitcvx,1),'Not all subjects were tested with crossvalidation!?')
            Nsamp = size(ASgt,1)/2;
            [pvl,u1,u2] = AUCMannWhitney(AS_AUCcvx(1),Nsamp,Nsamp,'right',AUC_AS(1),Nsamp,Nsamp)
            
            if pvl < 0.05
                fprintf('%i-fold cross-validation AUC (%1.2f) is significantly different (p = %1.3f) from the training AUC (%1.2f) for %s\n',...
                    k2fold,AS_AUCcvx(1),pvl,AUC_AS(1),feats{ift})
            else
                fprintf('No sig. difference between CV (AUC = %1.2f) and training (AUC = %1.2f) AUCs (p = %1.3f) for %s\n',...
                    AS_AUCcvx(1),AUC_AS(1),pvl,feats{ift})
            end
            
            pcv(ift-1) = pvl;
            u1cv(ift-1) = u1;
            u2cv(ift-1) = u2;
            AUC1cv(ift-1,:) = AS_AUCcvx; % cross-validated
            AUC2cv(ift-1,:) = AUC_AS; % training
                        
            thresh = unique(T_AS(intersect(find(X_AS(:,1) == opt_AS(1)),find(Y_AS(:,1) == opt_AS(2)))));
            assert(length(thresh)==1,'More than one threshold found')
            TP = sum(ASgt & (ASfit >= thresh));
            FP = sum(~ASgt & (ASfit >= thresh));
            TN = sum(~ASgt & (ASfit < thresh));
            FN = sum(ASgt & (ASfit < thresh));
            sens_AS = TP/(TP+FN);
            spec_AS = TN/(TN+FP);
            ppv_AS = TP/(TP+FP);
            acc_AS = (TP+TN)/(TP+FP+TN+FN);

            % now do perfcurve again to get 95% CIs for accuracy,
            % sensitivity, specificity, etc.
            
            [FALL_AS,ACCU_AS] = perfcurve(ASgt,ASfit,1,'XCrit','FPR','YCrit','accu','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_AS(1));
            try % this one will sometimes fail if AUC is very small
                [FALL_AS,PPV_AS] = perfcurve(ASgt,ASfit,1,'XCrit','FPR','YCrit','PPV','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_AS(1));
            catch
                if ppv_AS == 0
                    PPV_AS = [0 0 0];
                else
                    error('Unknown bug computing PPV_AS')
                end
            end
            [FALL_AS,RECALL_AS] = perfcurve(ASgt,ASfit,1,'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_AS(1));
            SPEC_AS = 1-X_AS(find(X_AS(:,1)==opt_AS(1),1),:);
            
            SPEC_AS = SPEC_AS([1 3 2]); % the CI bounds are swapped, flip them
            
            % sanity checks
            assert(single(ACCU_AS(1))==single(acc_AS),'Accuracies don''t match')
            assert(single(PPV_AS(1))==single(ppv_AS),'PPVs don''t match')
            assert(single(RECALL_AS(1))==single(sens_AS),'sensitivities (recall) don''t match')
            assert(single(SPEC_AS(1))==single(spec_AS),'specificities don''t match')
            
            % extract results
            results(1,1,ift-1,:) = AUC_AS.*100;
            results(1,2,ift-1,:) = ACCU_AS.*100;
            results(1,3,ift-1,:) = PPV_AS.*100;
            results(1,4,ift-1,:) = RECALL_AS.*100;
            results(1,5,ift-1,:) = SPEC_AS.*100;
                            
            fprintf('AS training set accuracy = %1.2f%%, AUC = %1.2f%%, precision = %1.2f%%, recall = %1.2f%%\n',...
                results(1,2,ift-1,1),results(1,1,ift-1,1),results(1,3,ift-1,1),results(1,4,ift-1,1))
            
            %%% TD %%%
            % Test to see if the curve has the correct number of points, 
            % N + 1, where N is the number of data sests
                        
            [X_TD,Y_TD,T_TD,AUC_TD,opt_TD,bs_TD] = myperfcurve(TDgt,TDfit,1,...
                'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Alpha',0.05);
            
            thresh = unique(T_TD(intersect(find(X_TD(:,1) == opt_TD(1)),find(Y_TD(:,1) == opt_TD(2))),1));
            assert(length(thresh)==1,'More than one threshold found')
            TP = sum(TDgt & (TDfit >= thresh));
            FP = sum(~TDgt & (TDfit >= thresh));
            TN = sum(~TDgt & (TDfit < thresh));
            FN = sum(TDgt & (TDfit < thresh));
            sens_TD = TP/(TP+FN);
            spec_TD = TN/(TN+FP);
            ppv_TD = TP/(TP+FP);
            acc_TD = (TP+TN)/(TP+FP+TN+FN);

            % now do perfcurve again to get 95% CIs for accuracy,
            % sensitivity, specificity, etc.
            
            % I set up try-catch blocks below for TD ... for some reason,
            % bugs only happen with the TD data, probably because they
            % have "perfect" EEGs that don't match with AS and DS
             
            [FALL_TD,ACCU_TD] = perfcurve(TDgt,TDfit,1,'XCrit','FPR','YCrit','accu','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_TD(1));
            %if we had a wrong accuracy because of the weird operating point
            if single(ACCU_TD(1)) ~= single(acc_TD) && opt_TD(1) == 0 
                [FALL_TD,ACCU_TD] = perfcurve(TDgt,TDfit,1,'XCrit','FPR','YCrit','accu','Nboot',Nboot,'Alpha',0.05);
                ACCU_TD = ACCU_TD(find(ACCU_TD(:,1) == acc_TD,1),:);
                [FALL_TD,PPV_TD] = perfcurve(TDgt,TDfit,1,'XCrit','FPR','YCrit','PPV','Nboot',Nboot,'Alpha',0.05);
                PPV_TD = PPV_TD(find(PPV_TD(:,1) == ppv_TD,1),:);
                [FALL_TD,RECALL_TD] = perfcurve(TDgt,TDfit,1,'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Alpha',0.05);
                RECALL_TD = RECALL_TD(find(RECALL_TD(:,1) == sens_TD,1),:);
                SPEC_TD = 1-X_TD;
                SPEC_TD = SPEC_TD(find(SPEC_TD(:,1) == spec_TD,1),:);
                SPEC_TD = SPEC_TD([1 3 2]); % the CI bounds are swapped, flip them
            else
                try % this one will sometimes fail if AUC is very small
                    [FALL_TD,PPV_TD] = perfcurve(TDgt,TDfit,1,'XCrit','FPR','YCrit','PPV','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_TD(1));
                    assert(~any(isnan(PPV_TD)),'NaNs in PPV or its CI') % another specific bug I've encountered, only for DS
                catch
                    if ppv_TD == 0
                        PPV_TD = [0 0 0];
                    else
                        [FALL_TD,PPV_TD] = perfcurve(TDgt,TDfit,1,'XCrit','FPR','YCrit','PPV','Nboot',Nboot,'Alpha',0.05);
                        PPV_TD = PPV_TD(find(PPV_TD(:,1) == ppv_TD,1),:);
                    end
                end
                [FALL_TD,RECALL_TD] = perfcurve(TDgt,TDfit,1,'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_TD(1));
                SPEC_TD = 1-X_TD(find(X_TD(:,1)==opt_TD(1),1),:);
                SPEC_TD = SPEC_TD([1 3 2]); % the CI bounds are swapped, flip them
            end
            
            % sanity checks
            assert(single(ACCU_TD(1))==single(acc_TD),'Accuracies don''t match')
            assert(single(PPV_TD(1))==single(ppv_TD),'PPVs don''t match')
            assert(single(RECALL_TD(1))==single(sens_TD),'sensitivities (recall) don''t match')
            assert(single(SPEC_TD(1))==single(spec_TD),'specificities don''t match')
            
            % extract results
            results(2,1,ift-1,:) = AUC_TD.*100;
            results(2,2,ift-1,:) = ACCU_TD.*100;
            results(2,3,ift-1,:) = PPV_TD.*100;
            results(2,4,ift-1,:) = RECALL_TD.*100;
            results(2,5,ift-1,:) = SPEC_TD.*100;
                            
            fprintf('TD training set accuracy = %1.2f%%, AUC = %1.2f%%, precision = %1.2f%%, recall = %1.2f%%\n',...
                results(2,2,ift-1,1),results(2,1,ift-1,1),results(2,3,ift-1,1),results(2,4,ift-1,1))
            
            %%% Dup15q syndrome %%% 
            % Test to see if the curve has the correct number of points, 
            % N + 1, where N is the number of data sests
           
            [X_DS,Y_DS,T_DS,AUC_DS,opt_DS,bs_DS] = myperfcurve(DSgt,DSfit,1,...
                'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Alpha',0.05);
                                   
            thresh = unique(T_DS(intersect(find(X_DS(:,1) == opt_DS(1)),find(Y_DS(:,1) == opt_DS(2)))));
            assert(length(thresh)==1,'More than one threshold found')
            TP = sum(DSgt & (DSfit >= thresh));
            FP = sum(~DSgt & (DSfit >= thresh));
            TN = sum(~DSgt & (DSfit < thresh));
            FN = sum(DSgt & (DSfit < thresh));
            sens_DS = TP/(TP+FN);
            spec_DS = TN/(TN+FP);
            ppv_DS = TP/(TP+FP);
            acc_DS = (TP+TN)/(TP+FP+TN+FN);

            % now do perfcurve again to get 95% CIs for accuracy,
            % sensitivity, specificity, etc.
            
            [FALL_DS,ACCU_DS] = perfcurve(DSgt,DSfit,1,'XCrit','FPR','YCrit','accu','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_DS(1));
            try % this one will sometimes fail if AUC is very small
                [FALL_DS,PPV_DS] = perfcurve(DSgt,DSfit,1,'XCrit','FPR','YCrit','PPV','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_DS(1));
                assert(~any(isnan(PPV_DS)),'NaNs in PPV or its CI') % another specific bug I've encountered, only for DS
            catch
                if ppv_DS == 0
                    PPV_DS = [0 0 0];
                else
                    [FALL_DS,PPV_DS] = perfcurve(DSgt,DSfit,1,'XCrit','FPR','YCrit','PPV','Nboot',Nboot,'Alpha',0.05);
                    PPV_DS = PPV_DS(find(PPV_DS(:,1) == ppv_DS,1),:);
                end
            end
            [FALL_DS,RECALL_DS] = perfcurve(DSgt,DSfit,1,'XCrit','FPR','YCrit','TPR','Nboot',Nboot,'Alpha',0.05,'Xvals',opt_DS(1));
            SPEC_DS = 1-X_DS(find(X_DS(:,1)==opt_DS(1),1),:);
            SPEC_DS = SPEC_DS([1 3 2]); % the CI bounds are swapped, flip them
            
            if any(isnan(PPV_DS)), keyboard, end
            
            % sanity checks
            assert(single(ACCU_DS(1))==single(acc_DS),'Accuracies don''t match')
            assert(single(PPV_DS(1))==single(ppv_DS),'PPVs don''t match')
            assert(single(RECALL_DS(1))==single(sens_DS),'sensitivities (recall) don''t match')
            assert(single(SPEC_DS(1))==single(spec_DS),'specificities don''t match')
            
            % extract results
            results(3,1,ift-1,:) = AUC_DS.*100;
            results(3,2,ift-1,:) = ACCU_DS.*100;
            results(3,3,ift-1,:) = PPV_DS.*100;
            results(3,4,ift-1,:) = RECALL_DS.*100;
            results(3,5,ift-1,:) = SPEC_DS.*100;
                            
            fprintf('DS training set accuracy = %1.2f%%, AUC = %1.2f%%, precision = %1.2f%%, recall = %1.2f%%\n',...
                results(3,2,ift-1,1),results(3,1,ift-1,1),results(3,3,ift-1,1),results(3,4,ift-1,1))
            fprintf('\n')
            
            % plot histograms of boostrapped distributions
            myfigure
            hold on
            if contains(feats{ift},'Entropy') % dark colors for entropy
                [H_AS] = histogram(bs_AS,0:0.01:1,'facecolor','#0072BD');
            else % light colors for spectral
                [H_AS] = histogram(bs_AS,0:0.01:1,'facecolor','c');
            end
                
            title(sprintf('AS AUC %i bootstraps, %s %s',Nboot,feats{ift},method))
            xlabel('AUC (%)')
            ylabel('Bootstrapped resamples')
            set(gca,'YScale','log')
            axis([0 1 0 10000])
            xticks([0 0.2 0.4 0.6 0.8 1])
            xticklabels([0 20 40 60 80 100])
            makefighandsome
            outname = sprintf('AS_AUC_hist_%i_%s_%s',Nboot,feats{ift},method);
            print('-dsvg',sprintf('./Figures/%s',outname))
            print('-dpng',sprintf('./Figures/%s',outname))
            
            myfigure
            hold on
            if contains(feats{ift},'Entropy') % dark colors for entropy
                [H_TD] = histogram(bs_TD,0:0.01:1,'facecolor','#77AC30');
            else % light colors for spectral
                [H_TD] = histogram(bs_TD,0:0.01:1,'facecolor','g');
            end
                
            title(sprintf('NT AUC %i bootstraps, %s %s',Nboot,feats{ift},method))
            xlabel('AUC (%)')
            ylabel('Bootstrapped resamples')
            set(gca,'YScale','log')
            axis([0 1 0 10000])
            xticks([0 0.2 0.4 0.6 0.8 1])
            xticklabels([0 20 40 60 80 100])
            makefighandsome
            outname = sprintf('NT_AUC_hist_%i_%s_%s',Nboot,feats{ift},method);
            print('-dsvg',sprintf('./Figures/%s',outname))
            print('-dpng',sprintf('./Figures/%s',outname))
            
            myfigure
            hold on
            if contains(feats{ift},'Entropy') % dark colors for entropy
                [H_DS] = histogram(bs_DS,0:0.01:1,'facecolor','#A2142F');
            else % light colors for spectral
                [H_DS] = histogram(bs_DS,0:0.01:1,'facecolor','r');
            end
                
            title(sprintf('Dup15q AUC %i bootstraps, %s %s',Nboot,feats{ift},method))
            xlabel('AUC (%)')
            ylabel('Bootstrapped resamples')
            set(gca,'YScale','log')
            axis([0 1 0 10000])
            xticks([0 0.2 0.4 0.6 0.8 1])
            xticklabels([0 20 40 60 80 100])
            makefighandsome
            outname = sprintf('DS_AUC_hist_%i_%s_%s',Nboot,feats{ift},method);
            print('-dsvg',sprintf('./Figures/%s',outname))
            print('-dpng',sprintf('./Figures/%s',outname))
            
            % get histogram data
            eval(sprintf('AS_%s_boots = H_AS.Values''',feats{ift}))
            eval(sprintf('TD_%s_boots = H_TD.Values''',feats{ift}))
            eval(sprintf('DS_%s_boots = H_DS.Values''',feats{ift}))

            xAS{ift-1} = X_AS(:,1);
            xTD{ift-1} = X_TD(:,1);
            xDS{ift-1} = X_DS(:,1);
            yAS{ift-1} = Y_AS(:,1);
            yTD{ift-1} = Y_TD(:,1);
            yDS{ift-1} = Y_DS(:,1);
                    
            myfigure2
            plot(X_AS(:,1),Y_AS(:,1),'LineWidth',2)
            plot(X_TD(:,1),Y_TD(:,1),'LineWidth',2)
            plot(X_DS(:,1),Y_DS(:,1),'LineWidth',2)
            plot([-0.02 1],[-0.02 1],'k:')
            legend({sprintf('AS, AUC = %2.1f%%',AUC_AS(1)*100),sprintf('TD, AUC = %2.1f%%',AUC_TD(1)*100),...
                sprintf('Dup15q, AUC = %2.1f%%',AUC_DS(1)*100),'Chance performance'},'location',...
                'northeastoutside','fontsize',12,'autoupdate','off')
            legend box off
            xlabel('False positive rate')
            ylabel('True positive rate')
            title(sprintf('ROC curve, %s with %s feat. selection',feats{ift},method),'fontsize',18)
            axis([0 1 0 1])
            xticks(0:0.2:1)
            yticks(0:0.2:1)
            makefighandsome
            axis square
            pause(0.01)
            print('-dpng',sprintf('./Figures/ROC_curve_%s_%s',feats{ift},method))

        end
        
        assert(all(single(results(:,:,:,1)) >= single(results(:,:,:,2)) & ...
            single(results(:,:,:,1)) <= single(results(:,:,:,3)),[1 2 3]),...
            'True values are not contained within 95% CIs!')
        
        AUC_bins = [0.005:0.01:0.995]';
        Tfig7 = table(AUC_bins,AS_fcEntropy_boots,AS_fcSpectral_boots,AS_scEntropy_boots,AS_scSpectralA_boots,AS_scSpectralR_boots,...
            TD_fcEntropy_boots,TD_fcSpectral_boots,TD_scEntropy_boots,TD_scSpectralA_boots,TD_scSpectralR_boots,...
            DS_fcEntropy_boots,DS_fcSpectral_boots,DS_scEntropy_boots,DS_scSpectralA_boots,DS_scSpectralR_boots);
        writetable(Tfig7,sprintf('./Fig7_table_%s.csv',method))
        
        
        %% Plot ROC curves separately for each group
        
        ro = [4 3 1 5 2]; % reorder to get the colors we want for each group
        
        switch method
            case 'LMM'
                openfig('./Figures/TD_ROC_curve_PCA.fig')
                hold on
                for i = ro
                    switch i
                        case 1
                            pclr = [0.9290 0.6940 0.1250];
                        case 2
                            pclr = [0.4660 0.6740 0.1880];
                        case 3
                            pclr = [0.8500 0.3250 0.0980];
                        case 4
                            pclr = [0 0.4470 0.7410];
                        case 5
                            pclr = [0.4940 0.1840 0.5560];
                        otherwise
                            error('Index exceeds number of conditions')
                    end
                    plot(xTD{i}',yTD{i}',':','color',pclr,'LineWidth',2)        
                end
            case 'PCA'
                myfigure
                for i = ro, plot(xTD{i}',yTD{i}','LineWidth',2), end
        end
        
        plot([-0.02 1],[-0.02 1],'k:')
        xlabel('False positive rate (%)')
        ylabel('True positive rate (%)')
        title(sprintf('TD ROC curves with %s feat. selection',method),'fontsize',18)
        axis([-0.02 1 -0.02 1])
        xticks(0:0.2:1)
        xticklabels([0:20:100])
        yticks(0:0.2:1)
        yticklabels([0:20:100])
        makefigpretty
        print('-dpng',sprintf('./Figures/TD_ROC_curve_%s.png',method))
        print('-dsvg',sprintf('./Figures/TD_ROC_curve_%s.svg',method))
        savefig(sprintf('./Figures/TD_ROC_curve_%s.fig',method))
        
       
        switch method
            case 'LMM'
                openfig('./Figures/AS_ROC_curve_PCA.fig')
                hold on
                for i = ro
                    switch i
                        case 1
                            pclr = [0.9290 0.6940 0.1250];
                        case 2
                            pclr = [0.4660 0.6740 0.1880];
                        case 3
                            pclr = [0.8500 0.3250 0.0980];
                        case 4
                            pclr = [0 0.4470 0.7410];
                        case 5
                            pclr = [0.4940 0.1840 0.5560];
                        otherwise
                            error('Index exceeds number of conditions')
                    end
                    plot(xAS{i}',yAS{i}',':','color',pclr,'LineWidth',2)        
                end
            case 'PCA'
                myfigure
                for i = ro, plot(xAS{i}',yAS{i}','LineWidth',2), end
        end
        
        plot([-0.02 1],[-0.02 1],'k:')
        xlabel('False positive rate (%)')
        ylabel('True positive rate (%)')
        title(sprintf('AS ROC curves with %s feat. selection',method),'fontsize',18)
        axis([-0.02 1 -0.02 1])
        xticks(0:0.2:1)
        xticklabels([0:20:100])
        yticks(0:0.2:1)
        yticklabels([0:20:100])
        makefigpretty
        print('-dpng',sprintf('./Figures/AS_ROC_curve_%s.png',method))
        print('-dsvg',sprintf('./Figures/AS_ROC_curve_%s.svg',method))
        savefig(sprintf('./Figures/AS_ROC_curve_%s.fig',method))
        
       
        switch method
            case 'LMM'
                openfig('./Figures/Dup15q_ROC_curve_PCA.fig')
                hold on
                for i = ro
                    switch i
                        case 1
                            pclr = [0.9290 0.6940 0.1250];
                        case 2
                            pclr = [0.4660 0.6740 0.1880];
                        case 3
                            pclr = [0.8500 0.3250 0.0980];
                        case 4
                            pclr = [0 0.4470 0.7410];
                        case 5
                            pclr = [0.4940 0.1840 0.5560];
                        otherwise
                            error('Index exceeds number of conditions')
                    end
                    plot(xDS{i}',yDS{i}',':','color',pclr,'LineWidth',2)        
                end
            case 'PCA'
                myfigure
                for i = ro, plot(xDS{i}',yDS{i}','LineWidth',2), end
        end
        
        plot([-0.02 1],[-0.02 1],'k:')
        xlabel('False positive rate (%)')
        ylabel('True positive rate (%)')
        title(sprintf('Dup15q ROC curves with %s feat. selection',method),'fontsize',18)
        axis([-0.02 1 -0.02 1])
        xticks(0:0.2:1)
        xticklabels([0:20:100])
        yticks(0:0.2:1)
        yticklabels([0:20:100])
        makefigpretty
        print('-dpng',sprintf('./Figures/Dup15q_ROC_curve_%s.png',method))
        print('-dsvg',sprintf('./Figures/Dup15q_ROC_curve_%s.svg',method))
        savefig(sprintf('./Figures/Dup15q_ROC_curve_%s.fig',method))
        
        
        %% Test for differences in AUC between measures
        % feature category names that we will use in the paper
        featnames = feats(2:end);
        
        % initialize vars
        tststr = cell(9,1);
        groupz = cell(9,1);
        e2sp = nan(9,1); % entropy to spectral comparison p-vales
        AUC1 = nan(9,1);
        AUC2 = nan(9,1);
        MWU1 = nan(9,1);
        MWU2 = nan(9,1);
        SampN = nan(9,1);
        winner = cell(9,1);
        
        % Consider the CV AUC separately for AS
        e2sp_cv = nan(3,1); % entropy to spectral comparison p-vales
        AUC1_cv = nan(3,1);
        AUC2_cv = nan(3,1);
        MWU1_cv = nan(3,1);
        MWU2_cv = nan(3,1);
        SampN_cv = nan(3,1);
        winner_cv = cell(3,1);
        groupz_cv = cell(3,1);
        teststr_cv = cell(3,1);

        cnt = 0;
        cntcv = 0;
        g = {'AS';'TD';'DS'};
        for ift = 1:size(results,3)
            for jft = 1:size(results,3)
                for igrp = 1:3
                    if contains(featnames{ift},'Entropy') && contains(featnames{jft},'Spectral') && ...
                            ((contains(featnames{ift},'sc') && contains(featnames{jft},'sc')) || ...
                            (contains(featnames{ift},'fc') && contains(featnames{jft},'fc')))
                        cnt = cnt + 1;
                        groupz{cnt} = g{igrp};
                        tststr{cnt} = sprintf('%s vs %s',featnames{ift},featnames{jft});
                        N = sum(strcmp(cellstr(T.Group),g{igrp}))/2;
                        if strcmp(g{igrp},'AS')
                            % For AS only, also use the CV AUC
                            % Mann-Whitney
                            cntcv = cntcv + 1;
                            groupz_cv{cntcv} = g{igrp};
                            tststr_cv{cntcv} = sprintf('%s vs %s',featnames{ift},featnames{jft});
                            [pVal2,U1,U2] = AUCMannWhitney(AUC1cv(ift),N,N,'both',AUC1cv(jft),N,N);
                            e2sp_cv(cntcv) = pVal2;
                            AUC1_cv(cntcv) = AUC1cv(ift);
                            AUC2_cv(cntcv) = AUC1cv(jft);
                            MWU1_cv(cntcv) = U1;
                            MWU2_cv(cntcv) = U2;
                            SampN_cv(cnt) = N;
                            if U1 > U2
                                winner_cv{cntcv} = 'Entropy';
                            elseif U1 < U2
                                winner_cv{cntcv} = 'Spectral';
                            else
                                winner_cv{cntcv} = 'NA';
                            end
                        end
                            
                        % Mann-Whitney
                        [pVal2,U1,U2] = AUCMannWhitney(results(igrp,1,ift,1)./100,N,N,'both',results(igrp,1,jft,1)./100,N,N);
                        e2sp(cnt) = pVal2;
                        AUC1(cnt) = results(igrp,1,ift)./100;
                        AUC2(cnt) = results(igrp,1,jft)./100;
                        MWU1(cnt) = U1;
                        MWU2(cnt) = U2;
                        SampN(cnt) = N;
                        if U1 > U2
                            winner{cnt} = 'Entropy';
                        elseif U1 < U2
                            winner{cnt} = 'Spectral';
                        else
                            winner{cnt} = 'NA';
                        end
                    end
                end
            end
        end
                
        % Note: FDR correction is done separately with another script
        
        % Table for main manuscript
        Te2s = table(groupz,SampN,tststr,winner,e2sp,AUC1,AUC2,MWU1,MWU2,'VariableNames',...
            {'Cohort','N','Test','Larger AUC','MWU p-value','Entropy AUC','Spectral AUC',...
            'Entropy U','Spectral U'})
        writetable(Te2s,sprintf('./EntrVsSpecPvals_%s%i.csv',method,Nboot))
        
        
        Te2s_cv = table(groupz_cv,repmat(sum(ASgt),3,1),tststr_cv',winner_cv,...
            e2sp_cv,AUC1_cv,AUC2_cv,MWU1_cv,MWU2_cv,'VariableNames',...
            {'Cohort','N','Test','Larger AUC','MWU p-value','Entropy AUC','Spectral AUC',...
            'Entropy U','Spectral U'})
        writetable(Te2s_cv,sprintf('./EntrVsSpecPvals_10foldCV_%s%i.csv',method,Nboot))
        

        %% Save the output
        switch method
            case 'PCA'
                save EEG_ML_vars_PCA_MWU
            case 'LMM'
                save EEG_ML_vars_LMM_MWU
        end

        
        %% Create table of results
        
        allAUC = reshape(results(:,1,:,1),numel(results)/15,1)./100;
        allACC = reshape(results(:,2,:,1),numel(results)/15,1)./100;
        allPPV = reshape(results(:,3,:,1),numel(results)/15,1)./100; % same as precision
        allSENS = reshape(results(:,4,:,1),numel(results)/15,1)./100; % same as recall
        allSPEC = reshape(results(:,5,:,1),numel(results)/15,1)./100; % specificity
        
        features = [];
        for ift = 1:length(featnames)
            features = [features repmat(featnames(ift),1,3)];
        end
        
        groups = repmat({'Angelman';'Neurotypical';'Dup15q'},size(results,3),1);
        
        cls = repmat({'Training';'Validation';'Validation'},size(results,3),1);
        
        f2 = repmat(unique(feats(~contains(feats,'IV')))',size(unique(feats(~contains(feats,'IV'))),1),3)';
        F = reshape(f2,1,size(f2,1)*size(f2,2))';
        
        LAMBDA2 = repmat(LAMBDA',size(LAMBDA,1),3)';
        L = reshape(LAMBDA2,1,size(LAMBDA2,1)*size(LAMBDA2,2))';
        
        AUC_CIs = squeeze([results(1,1,:,2:3); results(2,1,:,2:3); results(3,1,:,2:3)]);
        AUC_CI1 = reshape(AUC_CIs(:,:,1),15,1);
        AUC_CI2 = reshape(AUC_CIs(:,:,2),15,1);
        
        ACC_CIs = squeeze([results(1,2,:,2:3); results(2,2,:,2:3); results(3,2,:,2:3)]);
        ACC_CI1 = reshape(ACC_CIs(:,:,1),15,1);
        ACC_CI2 = reshape(ACC_CIs(:,:,2),15,1);
        
        PPV_CIs = squeeze([results(1,3,:,2:3); results(2,3,:,2:3); results(3,3,:,2:3)]);
        PPV_CI1 = reshape(PPV_CIs(:,:,1),15,1);
        PPV_CI2 = reshape(PPV_CIs(:,:,2),15,1);
        
        SENS_CIs = squeeze([results(1,4,:,2:3); results(2,4,:,2:3); results(3,4,:,2:3)]);
        SENS_CI1 = reshape(SENS_CIs(:,:,1),15,1);
        SENS_CI2 = reshape(SENS_CIs(:,:,2),15,1);
        
        SPEC_CIs = squeeze([results(1,5,:,2:3); results(2,5,:,2:3); results(3,5,:,2:3)]);
        SPEC_CI1 = reshape(SPEC_CIs(:,:,1),15,1);
        SPEC_CI2 = reshape(SPEC_CIs(:,:,2),15,1);
        
        % Create strings with confidence intervals
        
        AUCstrs = cell(size(allAUC,1),1);
        ACCstrs = cell(size(allACC,1),1);
        PPVstrs = cell(size(allPPV,1),1);
        SENSstrs = cell(size(allSENS,1),1);
        SPECstrs = cell(size(allSPEC,1),1);
        
        % Get p-values for AUCs based on Mann-Whitney U
        MWp = nan(size(allAUC));
        MWU = nan(size(allAUC));
        for irow = 1:size(allAUC,1)
            switch groups{irow}
                case 'Angelman'
                    [P,U] = AUCMannWhitney(allAUC(irow),sum(ASidx)/2,sum(ASidx)/2,'right');
                    MWp(irow) = P;
                    MWU(irow) = U;
                case 'Neurotypical'
                    [P,U] = AUCMannWhitney(allAUC(irow),sum(TDidx)/2,sum(TDidx)/2,'right');
                    MWp(irow) = P;
                    MWU(irow) = U;
                case 'Dup15q'
                    [P,U] = AUCMannWhitney(allAUC(irow),sum(DSidx)/2,sum(DSidx)/2,'right');
                    MWp(irow) = P;
                    MWU(irow) = U;
            end
        end
        
        for irow = 1:size(allAUC,1)
            AUCstrs{irow,1} = sprintf('%2.1f%% (%2.1f%% - %2.1f%%)',allAUC(irow).*100,AUC_CI1(irow),AUC_CI2(irow));
            ACCstrs{irow,1} = sprintf('%2.1f%% (%2.1f%% - %2.1f%%)',allACC(irow).*100,ACC_CI1(irow),ACC_CI2(irow));
            PPVstrs{irow,1} = sprintf('%2.1f%% (%2.1f%% - %2.1f%%)',allPPV(irow).*100,PPV_CI1(irow),PPV_CI2(irow));
            SENSstrs{irow,1} = sprintf('%2.1f%% (%2.1f%% - %2.1f%%)',allSENS(irow).*100,SENS_CI1(irow),SENS_CI2(irow));
            SPECstrs{irow,1} = sprintf('%2.1f%% (%2.1f%% - %2.1f%%)',allSPEC(irow).*100,SPEC_CI1(irow),SPEC_CI2(irow));
        end
            
        
        % create big table (for supplement)
        Ta = table(F,cls,groups,L,MWp,MWU,AUCstrs,ACCstrs,PPVstrs,SENSstrs,SPECstrs,...
            'VariableNames',{'Features','Classifcation','Group','Lambda',...
            'MWU p-values','MW U (test-statistic)','AUC (95% CI)',...
            'Accuracy (95% CI)','Precision (95% CI)','Recall (95% CI)','Specificity (95% CI)'})
        writetable(Ta,sprintf('ClassificationResults%s_%i_resamples.csv',method,Nboot))
               
       % Just to be safe
        clear X_AS Y_AS T_AS AUC_AS opt_AS X_TD Y_TD T_TD AUC_TD opt_TD ...
            X_DS Y_DS T_DS AUC_DS opt_DS

        
        %% PCA coefficients
        
        figure('color','w')
        h=subplot(1,1,1);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 0.5]);
        stem([PC_coefs{2,2}(:,1); PC_coefs{2,3}(:,1); PC_coefs{2,1}(:,1); ...
            PC_coefs{2,4}(:,1); PC_coefs{2,5}(:,1)],'filled','linewidth',2)
        xticks(1:length(varnames))
        xticklabels(varnames(frstvar:end))
        ylim([-1 1])
        xtickangle(45)
        ylabel('PC_{1} loadings')
        
        box off
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', 20)
        set(xAX,'color','k')
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', 20)
        set(yAX,'color','k')
        set(gca, 'TickDir', 'out')
        set(gcf,'color','w')
        set(gca,'Layer','top')
        axis normal
        print('-dsvg','./Figures/PCA_coefficeints_stem.svg')
        
        %%
        for i = 1:size(PC_coefs,2)
            myfigure
            mypcolor(flipud(PC_coefs{2,i}))
            xticks(0.5:size(PC_coefs{2,i},2)+0.5)
            ve = round(PC_coefs{3,i},1);
            xl = cell(1,size(PC_coefs{2,i},1)+1);
            for j = 1:size(PC_coefs{2,i},1)
                xl{j+1} = sprintf('%2.1f%% \sigma^{2}',PC_coefs{3,i}(j));
            end
            xticklabels(xl)
            xtickangle(45)
            yticks(0.5:size(PC_coefs{2,i},1)+0.5)
            yticklabels([cell(1,1) fliplr(varnames(strcmp(vartype,PC_coefs{1,i})))])
            caxis([-1 1])
            title(PC_coefs{1,i})
            makefighandsome
            axis square
            axis([1 size(PC_coefs{2,i},2)+1 1 size(PC_coefs{2,i},1)+1])
            print('-dpng',sprintf('./Figures/PC%i_coefs.png',i))
            print('-dsvg',sprintf('./Figures/PC%i_coefs.svg',i))
        end

        
        %% Table with multivariate features (all channels)
        
        varnamesx = {'Sleep','Conscious','Group','Subject','Sleep2',...
            'Conscious2','Group2','Age','mMSE','LZc',...
            'CTW','PermEn8','PermEn16','PermEn32','PermEn64',...
            'PermEn128','slow','delta1','delta2',...
            'theta','alpha','beta','slowR','delta1R','delta2R',...
            'thetaR','alphaR','betaR'};
        
        % Multichannel permutation entropy

                
        ASPermEnwkTau8MC  = squeeze(ASWakePermEn(:,1,:));
        ASPermEnspTau8MC  = squeeze(ASSleepPermEn(:,1,:));
        DSPermEnwkTau8MC  = squeeze(DSWakePermEn(:,1,:));
        DSPermEnspTau8MC  = squeeze(DSSleepPermEn(:,1,:));
        TDPermEnwkTau8MC  = squeeze(TDWakePermEn(:,1,:));
        TDPermEnspTau8MC  = squeeze(TDSleepPermEn(:,1,:));
         
        ASPermEnwkTau16MC  = squeeze(ASWakePermEn(:,2,:));
        ASPermEnspTau16MC  = squeeze(ASSleepPermEn(:,2,:));
        DSPermEnwkTau16MC  = squeeze(DSWakePermEn(:,2,:));
        DSPermEnspTau16MC  = squeeze(DSSleepPermEn(:,2,:));
        TDPermEnwkTau16MC  = squeeze(TDWakePermEn(:,2,:));
        TDPermEnspTau16MC  = squeeze(TDSleepPermEn(:,2,:));
        
        ASPermEnwkTau32MC  = squeeze(ASWakePermEn(:,3,:));
        ASPermEnspTau32MC  = squeeze(ASSleepPermEn(:,3,:));
        DSPermEnwkTau32MC  = squeeze(DSWakePermEn(:,3,:));
        DSPermEnspTau32MC  = squeeze(DSSleepPermEn(:,3,:));
        TDPermEnwkTau32MC  = squeeze(TDWakePermEn(:,3,:));
        TDPermEnspTau32MC  = squeeze(TDSleepPermEn(:,3,:));
                
        ASPermEnwkTau64MC  = squeeze(ASWakePermEn(:,4,:));
        ASPermEnspTau64MC  = squeeze(ASSleepPermEn(:,4,:));
        DSPermEnwkTau64MC  = squeeze(DSWakePermEn(:,4,:));
        DSPermEnspTau64MC  = squeeze(DSSleepPermEn(:,4,:));
        TDPermEnwkTau64MC  = squeeze(TDWakePermEn(:,4,:));
        TDPermEnspTau64MC  = squeeze(TDSleepPermEn(:,4,:));
        
        ASPermEnwkTau128MC  = squeeze(ASWakePermEn(:,5,:));
        ASPermEnspTau128MC  = squeeze(ASSleepPermEn(:,5,:));
        DSPermEnwkTau128MC  = squeeze(DSWakePermEn(:,5,:));
        DSPermEnspTau128MC  = squeeze(DSSleepPermEn(:,5,:));
        TDPermEnwkTau128MC  = squeeze(TDWakePermEn(:,5,:));
        TDPermEnspTau128MC  = squeeze(TDSleepPermEn(:,5,:));
        
        
        Tx = table(issleep,isconscious,group,subjects,sleepBin,groupBin,consciousBin,...
            [ASages'; TDages'; DSages'; ASages'; TDages'; DSages'],...
            [ASWakeMSE'; TDWakeMSE'; DSWakeMSE'; ASSleepMSE'; TDSleepMSE'; DSSleepMSE'],...
            [ASWakeLZc'; TDWakeLZc'; DSWakeLZc'; ASSleepLZc'; TDSleepLZc'; DSSleepLZc'],...
            [ASWakeCTW'; TDWakeCTW'; DSWakeCTW'; ASSleepCTW'; TDSleepCTW'; DSSleepCTW'],...
            [ASPermEnwkTau8MC'; TDPermEnwkTau8MC'; DSPermEnwkTau8MC'; ASPermEnspTau8MC'; TDPermEnspTau8MC'; DSPermEnspTau8MC'],...
            [ASPermEnwkTau16MC'; TDPermEnwkTau16MC'; DSPermEnwkTau16MC'; ASPermEnspTau16MC'; TDPermEnspTau16MC'; DSPermEnspTau16MC'],...
            [ASPermEnwkTau32MC'; TDPermEnwkTau32MC'; DSPermEnwkTau32MC'; ASPermEnspTau32MC'; TDPermEnspTau32MC'; DSPermEnspTau32MC'],...
            [ASPermEnwkTau64MC'; TDPermEnwkTau64MC'; DSPermEnwkTau64MC'; ASPermEnspTau64MC'; TDPermEnspTau64MC'; DSPermEnspTau64MC'],...
            [ASPermEnwkTau128MC'; TDPermEnwkTau128MC'; DSPermEnwkTau128MC'; ASPermEnspTau128MC'; TDPermEnspTau128MC'; DSPermEnspTau128MC'],...
            [slow_wk_AS'; slow_wk_TD'; slow_wk_DS'; slow_wk_AS'; slow_sp_TD'; slow_sp_DS'],...
            [delta1_wk_AS'; delta1_wk_TD'; delta1_wk_DS'; delta1_sp_AS'; delta1_sp_TD'; delta1_sp_DS'],...
            [delta2_wk_AS'; delta2_wk_TD'; delta2_wk_DS'; delta2_sp_AS'; delta2_sp_TD'; delta2_sp_DS'],...
            [theta_wk_AS'; theta_wk_TD'; theta_wk_DS'; theta_sp_AS'; theta_sp_TD'; theta_sp_DS'],...
            [alpha_wk_AS'; alpha_wk_TD'; alpha_wk_DS'; alpha_sp_AS'; alpha_sp_TD'; alpha_sp_DS'],...
            [beta_wk_AS'; beta_wk_TD'; beta_wk_DS'; beta_sp_AS'; beta_sp_TD'; beta_sp_DS'],...
            [slow_rel_wk_AS'; slow_rel_wk_TD'; slow_rel_wk_DS'; slow_rel_sp_AS';slow_rel_sp_TD'; slow_rel_sp_DS'],...
            [delta1_rel_wk_AS'; delta1_rel_wk_TD'; delta1_rel_wk_DS'; delta1_rel_sp_AS';delta1_rel_sp_TD'; delta1_rel_sp_DS'],...
            [delta2_rel_wk_AS'; delta2_rel_wk_TD'; delta2_rel_wk_DS'; delta2_rel_sp_AS';delta2_rel_sp_TD'; delta2_rel_sp_DS'],...
            [theta_rel_wk_AS'; theta_rel_wk_TD'; theta_rel_wk_DS'; theta_rel_sp_AS';theta_rel_sp_TD'; theta_rel_sp_DS'],...
            [alpha_rel_wk_AS'; alpha_rel_wk_TD'; alpha_rel_wk_DS'; alpha_rel_sp_AS';alpha_rel_sp_TD'; alpha_rel_sp_DS'],...
            [beta_rel_wk_AS'; beta_rel_wk_TD'; beta_rel_wk_DS'; beta_rel_sp_AS';beta_rel_sp_TD'; beta_rel_sp_DS'],'VariableNames',varnamesx);
       
            % remove rows with NaNs
            Tx(any(isnan(Tx{:,max(find(strcmp(vartype,'IV')))+1:end}')),:) = []
        
            writetable(Tx,'MultichannelEEGFeatures.csv')
            
            %% AUC topographies
            
            Tx = readtable('./MultichannelEEGfeatures.csv');
            TxVars = properties(Tx);
            Tvars = properties(T);
                        
            % indvidual channel results: cohort x feature x channel
            chresults = nan(3,length(feats)-3,19);
            
            % variables belonging to each feature category
            % add underscores after absolute power names so that relative
            % power doesn't also get selected
            whichvars = { {varnames{frstvar:frstvar+7}}, {'slow_','delta1_',...
                'delta2_','theta_','alpha_','beta_'}, {varnames{33:38}} };
            
            for ich = 1:nchan
                for ift = 4:length(feats)
                switch method
                    case 'LMM'
                        feat  = feats{ift};
                        fprintf('\n%s FEATURES ONLY using %s method\n',upper(feat),method)
                        fts = find(ismember(vartype,feat));
                        ranks = find(abs(ASsleepES(1,ismember(vartype,feat))) >= betathresh); % use features with abs(beta) >= 0.5
                        selected = Tvars(fts(ranks));
                        pickmeidx = [];
                        for isl = 1:length(selected)
                            tmp = find(contains(TxVars,selected{isl}),ich);
                            TxVars(tmp(end))
                            assert(contains(TxVars(tmp(end)),num2str(ich)),'Wrong channel')
                            pickmeidx = [pickmeidx tmp(end)];
                        end
                        predALL = Tx{:,pickmeidx};
                        predAS = Tx{ASidx,pickmeidx };
                        predTD = Tx{TDidx,pickmeidx };
                        predDS = Tx{DSidx,pickmeidx };
                        fprintf('\n')
                    case 'PCA'
                        feat  = feats{ift};
                        fprintf('\n%s FEATURES ONLY using %s method\n',upper(feat),method)
                        howmany = find(cumsum(PC_coefs{3,1}) >= 90,1); % number of PCs
                        loadings = PC_coefs{2,contains(PC_coefs(1,:),feat)}(:,1:howmany);
                        tmp = find(contains(TxVars,whichvars{ift-3}));
                        TxVars{tmp(ich:19:end)}
                        eegvars = Tx{:,tmp(ich:19:end)};
                        chpcs = eegvars*loadings;
                        
                        predAS = chpcs(ASidx,:);
                        predTD = chpcs(TDidx,:);
                        predDS = chpcs(DSidx,:);
                end

                [B0,FitInfo,CV] = mylassoglm(predAS,categorical(Tx.Conscious(ASidx)),...
                    'binomial','Link','logit','CV',kfold,'Balanced',cat(3,mytraining,mytest));

                % Check that data from the same subjects are on same side of the CV partion
                IS = nan(1,kfold);
                for ifold = 1:kfold
                    St = S(logical(CV.training(:,ifold))); % training subjects
                    Sv = S(~logical(CV.training(:,ifold))); % test subjects
                    IS(ifold) = length(intersect(St,Sv));
                end

                assert(max(IS)==0,'At least one fold had data from same subject on different sides of CV partition')

                % Examine the cross-validation plot to see the effect of the Lambda regularization parameter.
                lassoPlot(B0,FitInfo,'plottype','CV');
                legend('show') % Show legend

                [~,useme] = min(FitInfo.Deviance); % find optimal regularization
                lambda = FitInfo.Lambda(useme); % regularization parameter
                LAMBDA(ift-1) = lambda;

                % Do logistic regression just on AS data using the chosen lambda value
                [B,LassoFit] = lassoglm(predAS,categorical(T.Conscious(ASidx)),'binomial','Link','logit','Lambda',lambda);

                coef = [LassoFit.Intercept; B];

                % apply model fit with optimal regularization to validation sets
                ASfit = glmval(coef,predAS,'logit','Constant','on');
                TDfit = glmval(coef,predTD,'logit','Constant','on');
                DSfit = glmval(coef,predDS,'logit','Constant','on');

                %%% ANGELMAN SYNDROME %%%
                % Test to see if the curve has the correct number of points, 
                % N + 1, where N is the number of data sests

                [X_AS,Y_AS,T_AS,AUC_AS,opt_AS] = perfcurve(ASgt,ASfit,1,...
                    'XCrit','FPR','YCrit','TPR');

                % extract results
                chresults(1,ift-3,ich) = AUC_AS.*100;

                %%% TD %%%
                % Test to see if the curve has the correct number of points, 
                % N + 1, where N is the number of data sests

                [X_TD,Y_TD,T_TD,AUC_TD,opt_TD] = perfcurve(TDgt,TDfit,1,...
                    'XCrit','FPR','YCrit','TPR');

                % extract results
                chresults(2,ift-3,ich) = AUC_TD.*100;

                %%% Dup15q syndrome %%% 
                % Test to see if the curve has the correct number of points, 
                % N + 1, where N is the number of data sests

                [X_DS,Y_DS,T_DS,AUC_DS,opt_DS] = perfcurve(DSgt,DSfit,1,...
                    'XCrit','FPR','YCrit','TPR');

                % extract results
                chresults(3,ift-3,ich) = AUC_DS.*100;

                end
            end
            
            for ift = 1:3
                myfigure
                plot_topo_AS(squeeze(chresults(1,ift,:)))
                caxis([0 100])
                colormap parula
                mycolorbar
                title(sprintf('AS AUC %s %s',feats{ift+3},method),'fontsize',18)
                print('-dpng',sprintf('./Figures/AS_AUC_%s_%s.png',feats{ift+3},method))
                
                myfigure
                plot_topo_AS(squeeze(chresults(2,ift,:)))
                caxis([0 100])
                colormap parula
                mycolorbar
                title(sprintf('TD AUC %s %s',feats{ift+3},method),'fontsize',18)
                print('-dpng',sprintf('./Figures/TD_AUC_%s_%s.png',feats{ift+3},method))
                
                myfigure
                plot_topo_AS(squeeze(chresults(3,ift,:)))
                caxis([0 100])
                colormap parula
                mycolorbar
                title(sprintf('DS AUC %s %s',feats{ift+3},method),'fontsize',18)
                print('-dpng',sprintf('./Figures/DS_AUC_%s_%s.png',feats{ift+3},method))
            end
            
            %% Make connectivity table
             
            % Make table
            
            wSMI8 = [squeeze(ASWakewSMIrs(:,1,:)) squeeze(TDWakewSMIrs(:,1,:)) ...
                squeeze(DSWakewSMIrs(:,1,:)) squeeze(ASSleepwSMIrs(:,1,:)) ...
                squeeze(TDSleepwSMIrs(:,1,:)) squeeze(DSSleepwSMIrs(:,1,:))];
            
            wSMI16 = [squeeze(ASWakewSMIrs(:,2,:)) squeeze(TDWakewSMIrs(:,2,:)) ...
                squeeze(DSWakewSMIrs(:,2,:)) squeeze(ASSleepwSMIrs(:,2,:)) ...
                squeeze(TDSleepwSMIrs(:,2,:)) squeeze(DSSleepwSMIrs(:,2,:))];
            
            wSMI32 = [squeeze(ASWakewSMIrs(:,3,:)) squeeze(TDWakewSMIrs(:,3,:)) ...
                squeeze(DSWakewSMIrs(:,3,:)) squeeze(ASSleepwSMIrs(:,3,:)) ...
                squeeze(TDSleepwSMIrs(:,3,:)) squeeze(DSSleepwSMIrs(:,3,:))];
            
            wSMI64 = [squeeze(ASWakewSMIrs(:,4,:)) squeeze(TDWakewSMIrs(:,4,:)) ...
                squeeze(DSWakewSMIrs(:,4,:)) squeeze(ASSleepwSMIrs(:,4,:)) ...
                squeeze(TDSleepwSMIrs(:,4,:)) squeeze(DSSleepwSMIrs(:,4,:))];
            
            wSMI128 = [squeeze(ASWakewSMIrs(:,5,:)) squeeze(TDWakewSMIrs(:,5,:)) ...
                squeeze(DSWakewSMIrs(:,5,:)) squeeze(ASSleepwSMIrs(:,5,:)) ...
                squeeze(TDSleepwSMIrs(:,5,:)) squeeze(DSSleepwSMIrs(:,5,:))];
            
            dwPLI1 = [squeeze(ASWakedwPLIrs(:,1,:)) squeeze(TDWakedwPLIrs(:,1,:)) ...
                squeeze(DSWakedwPLIrs(:,1,:)) squeeze(ASSleepdwPLIrs(:,1,:)) ...
                squeeze(TDSleepdwPLIrs(:,1,:)) squeeze(DSSleepdwPLIrs(:,1,:))];
            
            dwPLI2 = [squeeze(ASWakedwPLIrs(:,2,:)) squeeze(TDWakedwPLIrs(:,2,:)) ...
                squeeze(DSWakedwPLIrs(:,2,:)) squeeze(ASSleepdwPLIrs(:,2,:)) ...
                squeeze(TDSleepdwPLIrs(:,2,:)) squeeze(DSSleepdwPLIrs(:,2,:))];
            
            dwPLI3 = [squeeze(ASWakedwPLIrs(:,3,:)) squeeze(TDWakedwPLIrs(:,3,:)) ...
                squeeze(DSWakedwPLIrs(:,3,:)) squeeze(ASSleepdwPLIrs(:,3,:)) ...
                squeeze(TDSleepdwPLIrs(:,3,:)) squeeze(DSSleepdwPLIrs(:,3,:))];
            
            dwPLI4 = [squeeze(ASWakedwPLIrs(:,4,:)) squeeze(TDWakedwPLIrs(:,4,:)) ...
                squeeze(DSWakedwPLIrs(:,4,:)) squeeze(ASSleepdwPLIrs(:,4,:)) ...
                squeeze(TDSleepdwPLIrs(:,4,:)) squeeze(DSSleepdwPLIrs(:,4,:))];
            
            dwPLI5 = [squeeze(ASWakedwPLIrs(:,5,:)) squeeze(TDWakedwPLIrs(:,5,:)) ...
                squeeze(DSWakedwPLIrs(:,5,:)) squeeze(ASSleepdwPLIrs(:,5,:)) ...
                squeeze(TDSleepdwPLIrs(:,5,:)) squeeze(DSSleepdwPLIrs(:,5,:))];
                   
            dwPLI6 = [squeeze(ASWakedwPLIrs(:,6,:)) squeeze(TDWakedwPLIrs(:,6,:)) ...
                squeeze(DSWakedwPLIrs(:,6,:)) squeeze(ASSleepdwPLIrs(:,6,:)) ...
                squeeze(TDSleepdwPLIrs(:,6,:)) squeeze(DSSleepdwPLIrs(:,6,:))];
            
            Tfc = table(wSMI8',wSMI16',wSMI32',wSMI64',wSMI128',dwPLI1',dwPLI2',...
                dwPLI3',dwPLI4',dwPLI5',dwPLI6','VariableNames',{'wSMI8','wSMI16',...
                'wSMI32','wSMI64','wSMI128','dwPLIslow','dwPLIdelta1','dwPLIdelta2',...
                'dwPLItheta','dwPLIalpha','dwPLIbeta'});
            
            writetable(Tfc,'./FunctionalConnectivity.csv') % write table
            % we need to read the saved file so that each channel combo is its own column
            
            % AUC connectivity 
            
            Tfc = readtable('./FunctionalConnectivity.csv');
            % Change variable names
            for ivar = 1:length(Tfc.Properties.VariableNames)
                tmp = Tfc.Properties.VariableNames{ivar}; % read variable name
                varnum = str2double(tmp(find(tmp=='_') + 1:end));
                ROW = ceil(varnum/nchan);
                COL = mod(varnum,nchan);
                if COL == 0, COL = 19; end
                Tfc.Properties.VariableNames{ivar} = sprintf('%s_%s_%s',...
                    tmp(1:find(tmp=='_')-1),lay.label{ROW},lay.label{COL});
            end
                         
            Tvars = properties(T);
            TfcVars = Tfc.Properties.VariableNames;
                       
            % add IV columns if they're not already there
            if ~contains(Tfc.Properties.VariableNames,'Conscious')
                % remove rows that are just nans
                Tfc(all(isnan(Tfc{:,:}),2),:) = [];
                addme = find(contains(vartype,'IV'));
                for icol = 1:length(addme)
                    eval(sprintf('Tfc.%s = T.%s;',T.Properties.VariableNames{icol},...
                        T.Properties.VariableNames{icol}))
                end
            end
                        
            % indvidual channel results: cohort x feature x channel
            cxresults = nan(3,2,nchan,nchan);
            
            % variables belonging to each feature category
            % add underscores after absolute power names so that relative
            % power doesn't also get selected
            whichvars = { {varnames{frstvar:frstvar+7}}, {'slow_','delta1_',...
                'delta2_','theta_','alpha_','beta_'}, {varnames{33:38}} };
            
            for ich = 1:nchan
                for jch = 1:nchan
                    if ich > jch && ED(ich,jch) > 80
                        for ift = 2:3
                        switch method
                            case 'LMM'
                                feat  = feats{ift};
                                fprintf('\n%s features only using %s method\n',feat,method)
                                fts = find(ismember(vartype,feat));
                                ranks = find(abs(ASsleepES(1,ismember(vartype,feat))) >= betathresh); % use features with abs(beta) >= 0.5
                                selected = Tvars(fts(ranks))
                                if ED(ich,jch) < 130 % short-range
                                    selected = selected(contains(selected,'SR'));
                                else
                                    selected = selected(contains(selected,'LR'));
                                end
                                
                                pickmeidx = [];
                                for isl = 1:length(selected)
                                    tmp = find(contains(TfcVars,selected{isl}(3:end)) & ...
                                        contains(TfcVars,lay.label{ich}) & contains(TfcVars,lay.label{jch}) );
                                    TfcVars(tmp)
                                    assert(length(tmp)==2,'More than 2 results found')
                                    assert(strcmp(sort(TfcVars{tmp(1)}),sort(TfcVars{tmp(2)})),...
                                        'Results are not just the same variables with ith and jth channels flipped')
                                    pickmeidx = [pickmeidx tmp(1)];
                                end
                                
                                predALL = Tfc{:,pickmeidx};
                                predAS = Tfc{ASidx,pickmeidx };
                                predTD = Tfc{TDidx,pickmeidx };
                                predDS = Tfc{DSidx,pickmeidx };
                                fprintf('\n')
                            case 'PCA'
                                feat  = feats{ift};
                                fprintf('\n%s features only using %s method\n',feat,method)
                                fts = find(ismember(vartype,feat));
                                loadings = PC_coefs{2,ift-1};
                                % Just the PCs that were used
                                nb_pcs = find(cumsum(PC_coefs{3,ift-1}) >= 90,1);
                                L = loadings(:,1:nb_pcs);
                                % Averaged short-range and long-range
                                L2 = mean(cat(3,L,flipud(L)),3);
                                L2 = L2(1:size(L2,1)/2,:); 
                                predALL = nan(size(Tfc,1),size(L2,1),nb_pcs);
                                for ipc = 1:size(L2,2) % for each component                               
                                    switch feat
                                        case 'fcEntropy'
                                            cn = {'wSMI8','wSMI16','wSMI32','wSMI64','wSMI128'};
                                        case 'fcSpectral'
                                            cn = {'dwPLIslow','dwPLIdelta1',...
                                                'dwPLIdelta2','dwPLItheta','dwPLIalpha','dwPLIbeta'};
                                    end
                                    
                                    for ivar = 1:size(L2,1) % for each variable
                                        tmp = find(contains(TfcVars,cn{ivar}) & ...
                                            contains(TfcVars,lay.label{ich}) & contains(TfcVars,lay.label{jch}) );
                                        TfcVars(tmp)
                                        assert(length(tmp)==2,'More than 2 results found')
                                        assert(strcmp(sort(TfcVars{tmp(1)}),sort(TfcVars{tmp(2)})),...
                                            'Results are not just the same variables with ith and jth channels flipped')
                                        predALL(:,ivar,ipc) = Tfc{:,tmp(1)}.* L2(ivar,ipc);                                       
                                    end
                                end
                                
                                predALL = squeeze(mean(predALL,2)); % average across the variables that contribute to the PC
                                predAS = predALL(ASidx,:);
                                predTD = predALL(TDidx,:);
                                predDS = predALL(DSidx,:);
                        end
                        
                        [B0,FitInfo,CV] = mylassoglm(predAS,categorical(Tfc.Conscious(ASidx)),...
                            'binomial','Link','logit','CV',kfold,'Balanced',cat(3,mytraining,mytest));
 
                        % Check that data from the same subjects are on same side of the CV partion
                        IS = nan(1,kfold);
                        for ifold = 1:kfold
                            St = S(logical(CV.training(:,ifold))); % training subjects
                            Sv = S(~logical(CV.training(:,ifold))); % test subjects
                            IS(ifold) = length(intersect(St,Sv));
                        end

                        assert(max(IS)==0,'At least one fold had data from same subject on different sides of CV partition')

                        [~,useme] = min(FitInfo.Deviance); % find optimal regularization
                        lambda = FitInfo.Lambda(useme); % regularization parameter
                        LAMBDA(ift-1) = lambda;

                        % Do logistic regression just on AS data using the chosen lambda value
                        [B,LassoFit] = lassoglm(predAS,categorical(T.Conscious(ASidx)),'binomial','Link','logit','Lambda',lambda);

                        coef = [LassoFit.Intercept; B];
                        
                        % apply model fit with optimal regularization to validation sets
                        ASfit = glmval(coef,predAS,'logit','Constant','on');
                        TDfit = glmval(coef,predTD,'logit','Constant','on');
                        DSfit = glmval(coef,predDS,'logit','Constant','on');

                        %%% ANGELMAN SYNDROME %%%
                        % Test to see if the curve has the correct number of points, 
                        % N + 1, where N is the number of data sests

                        [X_AS,Y_AS,T_AS,AUC_AS,opt_AS] = perfcurve(ASgt,ASfit,1,...
                            'XCrit','FPR','YCrit','TPR');

                        % extract results
                        cxresults(1,ift-1,ich,jch) = AUC_AS.*100;

                        %%% TD %%%
                        % Test to see if the curve has the correct number of points, 
                        % N + 1, where N is the number of data sests

                        [X_TD,Y_TD,T_TD,AUC_TD,opt_TD] = perfcurve(TDgt,TDfit,1,...
                            'XCrit','FPR','YCrit','TPR');

                        % extract results
                        cxresults(2,ift-1,ich,jch) = AUC_TD.*100;

                        %%% Dup15q syndrome %%% 
                        % Test to see if the curve has the correct number of points, 
                        % N + 1, where N is the number of data sests

                        [X_DS,Y_DS,T_DS,AUC_DS,opt_DS] = perfcurve(DSgt,DSfit,1,...
                            'XCrit','FPR','YCrit','TPR');

                        % extract results
                        cxresults(3,ift-1,ich,jch) = AUC_DS.*100;
                        end
                    end
                end
            end
           
            save(sprintf('FC_channels_AUCs_%s',method'),'cxresults')
            
            %%
            
            myfigure2
            a = squeeze(cxresults(1,1,:,:));
            b = squeeze(cxresults(1,2,:,:));
            mypcolor(nanmean(cat(3,a,fliplr(flipud(b))),3))
            xticks([1:nchan]+0.5)
            xticklabels(lay.label)
            yticks([1:nchan]+0.5)
            yticklabels(flipud(lay.label))
            caxis([0 100])
            %mycolorbar
            axis([1 20 1 20])
            makefighandsome, axis square
            title(sprintf('AS AUCs %s',method))
            print('-dpng',sprintf('./Figures/ASConnMatAUC_%s.png',method))
            print('-dsvg',sprintf('./Figures/ASConnMatAUC_%s.svg',method))
            
            myfigure2
            a = squeeze(cxresults(2,1,:,:));
            b = squeeze(cxresults(2,2,:,:));
            mypcolor(nanmean(cat(3,a,fliplr(flipud(b))),3))
            xticks([1:nchan]+0.5)
            xticklabels(lay.label)
            yticks([1:nchan]+0.5)
            yticklabels(flipud(lay.label))
            caxis([0 100])
            %mycolorbar
            axis([1 20 1 20])
            makefighandsome, axis square
            title(sprintf('TD AUCs %s',method))
            print('-dpng',sprintf('./Figures/TDConnMatAUC_%s.png',method))
            print('-dsvg',sprintf('./Figures/TDConnMatAUC_%s.svg',method))
            
            myfigure2
            a = squeeze(cxresults(3,1,:,:));
            b = squeeze(cxresults(3,2,:,:));
            mypcolor(nanmean(cat(3,a,fliplr(flipud(b))),3))
            xticks([1:nchan]+0.5)
            xticklabels(lay.label)
            yticks([1:nchan]+0.5)
            yticklabels(flipud(lay.label))
            caxis([0 100])
            %mycolorbar
            axis([1 20 1 20])
            makefighandsome, axis square
            title(sprintf('DS AUCs %s',method))
            print('-dpng',sprintf('./Figures/DSConnMatAUC_%s.png',method))
            print('-dsvg',sprintf('./Figures/DSConnMatAUC_%s.svg',method))
            
        %% Look at univariate AUCs 
        
        rng(111) % set random number generator
        npc = 1;
        varnames2 = T.Properties.VariableNames(max(find(strcmp(vartype,'IV')))+1:end);
        vartype2 = [vartype T.Properties.VariableNames(contains(T.Properties.VariableNames,'PC'))];
        AUC = nan(3,length(varnames2));
        ACC = nan(3,length(varnames2));
        SENS = nan(3,length(varnames2));
        SPEC = nan(3,length(varnames2));
        PPV = nan(3,length(varnames2));
        lamb = nan(1,length(varnames2));
               
        for ivar = frstvar:length(T.Properties.VariableNames)
            ivar
            
            eval(sprintf('xvar = T.%s;',T.Properties.VariableNames{ivar}))
            
            % use UNREGULARIZED logistic regression
            glm = fitglm(xvar(ASidx),T.Conscious(ASidx),'Distribution','binomial','Link','logit')
            coef = glm.Coefficients{:,1};
            
            % apply model fit to validation sets
            ASfit = glmval(coef,xvar(ASidx),'logit','Constant','on');
            TDfit = glmval(coef,xvar(TDidx),'logit','Constant','on');
            DSfit = glmval(coef,xvar(DSidx),'logit','Constant','on');
            
            
            [X_AS,Y_AS,T_AS,AUC_AS,opt_AS] = perfcurve(ASgt,ASfit,1);
            thresh = T_AS(intersect(find(X_AS == opt_AS(1)),find(Y_AS == opt_AS(2))));
            acc_AS = sum(ASgt == (ASfit >= thresh))/length(ASfit)*100;
            fprintf('     AS training set accuracy = %1.2f%%, AUC = %1.2f%%\n',acc_AS,AUC_AS*100)
            TP = sum(ASgt & (ASfit >= thresh)); % true positives
            FP = sum(~ASgt & (ASfit >= thresh)); % false positives
            TN = sum(~ASgt & (ASfit < thresh)); % true negatives
            FN = sum(ASgt & (ASfit < thresh)); % false negatives
            assert(acc_AS==((TP+TN)/(TP+FP+TN+FN))*100,'Accuracy is wrong')
            
            if TP == 0 && FN == 0
                sens_AS = 0;
            else
                sens_AS = TP/(TP+FN);
            end
            
            if TN == 0 && FP == 0
                spec_AS = 0;
            else
                spec_AS = TN/(TN+FP);
            end
            
            if TP == 0 && FP == 0
                ppv_AS = 0;
            else
                ppv_AS = TP/(TP+FP);
            end
            
            [X_TD,Y_TD,T_TD,AUC_TD,opt_TD] = perfcurve(TDgt,TDfit,1);
            thresh = T_TD(intersect(find(X_TD == opt_TD(1)),find(Y_TD == opt_TD(2))));
            acc_TD = sum(TDgt == (TDfit >= thresh))/length(TDfit)*100;
            fprintf('    TD validation set accuracy = %1.2f%%, AUC = %1.2f%%\n',acc_TD,AUC_TD*100)
            TP = sum(TDgt & (TDfit >= thresh)); % true positives
            FP = sum(~TDgt & (TDfit >= thresh)); % false positives
            TN = sum(~TDgt & (TDfit < thresh)); % true negatives
            FN = sum(TDgt & (TDfit < thresh)); % false negatives
            assert(acc_TD==((TP+TN)/(TP+FP+TN+FN))*100,'Accuracy is wrong')
            
            if TP == 0 && FN == 0
                sens_TD = 0;
            else
                sens_TD = TP/(TP+FN);
            end
            
            if TN == 0 && FP == 0
                spec_TD = 0;
            else
                spec_TD = TN/(TN+FP);
            end
            
            if TP == 0 && FP == 0
                ppv_TD = 0;
            else
                ppv_TD = TP/(TP+FP);
            end
            
            [X_DS,Y_DS,T_DS,AUC_DS,opt_DS] = perfcurve(DSgt,DSfit,1);
            thresh = T_DS(intersect(find(X_DS == opt_DS(1)),find(Y_DS == opt_DS(2))));
            acc_DS = sum(DSgt == (DSfit >= thresh))/length(DSfit)*100;
            fprintf('    Dup15q validation set accuracy = %1.2f%%, AUC = %1.2f%%\n',acc_DS,AUC_DS*100)
            fprintf('\n')
            TP = sum(DSgt & (DSfit >= thresh)); % true positives
            FP = sum(~DSgt & (DSfit >= thresh)); % false positives
            TN = sum(~DSgt & (DSfit < thresh)); % true negatives
            FN = sum(DSgt & (DSfit < thresh)); % false negatives
            assert(acc_DS==((TP+TN)/(TP+FP+TN+FN))*100,'Accuracy is wrong')
            
            if TP == 0 && FN == 0
                sens_DS = 0;
            else
                sens_DS = TP/(TP+FN);
            end
            
            if TN == 0 && FP == 0
                spec_DS = 0;
            else
                spec_DS = TN/(TN+FP);
            end
            
            if TP == 0 && FP == 0
                ppv_DS = 0;
            else
                ppv_DS = TP/(TP+FP);
            end
            
            % Area under the ROC curve for each group
            AUC(1,ivar) = AUC_AS;
            AUC(2,ivar) = AUC_TD;
            AUC(3,ivar) = AUC_DS;
            
            ACC(1,ivar) = acc_AS/100;
            ACC(2,ivar) = acc_TD/100;
            ACC(3,ivar) = acc_DS/100;
            
            % Sensitivity for each group
            SENS(1,ivar) = sens_AS;
            SENS(2,ivar) = sens_TD;
            SENS(3,ivar) = sens_DS;
            
            % Specificity for each group
            SPEC(1,ivar) = spec_AS;
            SPEC(2,ivar) = spec_TD;
            SPEC(3,ivar) = spec_DS;
            
            % Precision for each group
            PPV(1,ivar) = ppv_AS;
            PPV(2,ivar) = ppv_TD;
            PPV(3,ivar) = ppv_DS; 
        end
        

        Tall = table(T.Properties.VariableNames(frstvar:end)',AUC(1,frstvar:end)',AUC(2,frstvar:end)',AUC(3,frstvar:end)',...
            ACC(1,frstvar:end)',ACC(2,frstvar:end)',ACC(3,frstvar:end)',...
            PPV(1,frstvar:end)',PPV(2,frstvar:end)',PPV(3,frstvar:end)',...
            SENS(1,frstvar:end)',SENS(2,frstvar:end)',SENS(3,frstvar:end)',...
            SPEC(1,frstvar:end)',SPEC(2,frstvar:end)',SPEC(3,frstvar:end)','VariableNames',...
            {'Vars','AS AUC','NT AUC','DS AUC',...
            'AS ACC','NT ACC','DS ACC',...
            'AS PPV','NT PPV','DS PPV',...
            'AS TPR','NT TPR','DS TPR',...
            'AS TNR','NT TNR','DS TNR'})
        
       
        
        
        %% Show classification stats for all variables
        tmp = properties(Tall);
        vartype3 = vartype2(frstvar:end);
        reorder = [find(contains(vartype3,'scEntropy')) find(contains(vartype3,'fcEntropy')) ...
            find(contains(vartype3,'scSpectralA')) find(contains(vartype3,'scSpectralR')) ...
            find(contains(vartype3,'fcSpectral'))];
        writetable(Tall(reorder,:),'./AllVariableScores.csv')
        tmp2 = Tall.Vars(reorder);
        varlabels = replace(tmp2,'_',' ');
        varlabels = replace(varlabels,'slow','s');
        varlabels = replace(varlabels,'delta','?');
        varlabels = replace(varlabels,'theta','?');
        varlabels = replace(varlabels,'alpha','?-?');
        varlabels = replace(varlabels,'beta','?');
        
        
        myfigure2
        % do it this way because pcolor function needs buffer on edges
        plttbl = ones(size(Tall,2),size(Tall,1)+1);
        plttbl(1:size(Tall,2)-1,1:size(Tall,1)) = flipud(Tall{reorder,2:end}');
        pcolor(plttbl)
        colormap parula
        caxis([0 1])
        axis([1 size(Tall,1)+1 1 size(Tall,2)])
        xticks(1.5:1:size(Tall,1)+1)
        xticklabels(varlabels)
        xtickangle(45)
        yticks(1.5:1:size(Tall,2))
        yticklabels(flipud(tmp(2:size(Tall,2))))
        box off
        set(gca,'linewidth',3)
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', 16)
        set(xAX,'color','k')
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', 20)
        set(yAX,'color','k')
        set(gca, 'TickDir', 'out')
        set(gcf,'color','w')
        set(gca,'Layer','top')
        axis normal
        print('-dpng','./Figures/AllVariableScores.png')
        print('-dsvg','./Figures/AllVariableScores.svg')
        
         %% Plot precision vs recall
                
        for igrp = 1:3
            myfigure2
            scatter(PPV(igrp,contains(vartype2,'scEntropy')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'scEntropy')&~idxLMM&~idxPC),270,CXclr,'fill',...
                'markeredgecolor','k','LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'fcEntropy')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'fcEntropy')&~idxLMM&~idxPC),270,MIclr,'fill',...
                'markeredgecolor','k','LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'fcSpectral')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'fcSpectral')&~idxLMM&~idxPC),270,PIclr,'fill',...
                'markeredgecolor','k','LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'scSpectralA')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'scSpectralA')&~idxLMM&~idxPC),270,ASclr,'fill',...
                'markeredgecolor','k','LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'scSpectralR')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'scSpectralR')&~idxLMM&~idxPC),270,RSclr,'fill',...
                'markeredgecolor','k','LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            
            buf = -0.01;

            % features selected by LMMs
            scatter(PPV(igrp,contains(vartype2,'scEntropy')&idxLMM),SENS(igrp,contains(vartype2,...
                'scEntropy')&idxLMM),300,'d','fill','markerfacecolor',CXclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'fcEntropy')&idxLMM),SENS(igrp,contains(vartype2,...
                'fcEntropy')&idxLMM),300,MIclr,'d','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'fcSpectral')&idxLMM),SENS(igrp,contains(vartype2,...
            'fcSpectral')&idxLMM),300,PIclr,'d','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'scSpectralA')&idxLMM),SENS(igrp,contains(vartype2,...
                'scSpectralA')&idxLMM),300,ASclr,'d','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'scSpectralR')&idxLMM),SENS(igrp,contains(vartype2,...
                'scSpectralR')&idxLMM),300,RSclr,'d','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            
            % First principal components
            scatter(PPV(igrp,contains(vartype2,'scEntropy')&idxPC),SENS(igrp,contains(vartype2,...
                'scEntropy')&idxPC),400,CXclr,'p','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'fcEntropy')&idxPC),SENS(igrp,contains(vartype2,...
                'fcEntropy')&idxPC),400,MIclr,'p','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'fcSpectral')&idxPC),SENS(igrp,contains(vartype2,...
                'fcSpectral')&idxPC),400,PIclr,'p','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'scSpectralA')&idxPC),SENS(igrp,contains(vartype2,...
                'scSpectralA')&idxPC),400,ASclr,'p','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'scSpectralR')&idxPC),SENS(igrp,contains(vartype2,...
                'scSpectralR')&idxPC),400,RSclr,'p','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            
            %legend boxoff
            
            title(g{igrp},'fontsize',18)
            xlabel('Precision (%)')
            ylabel('Recall (%)')
            xticks(0:0.2:1)
            xticklabels(0:20:100)
            yticks(0:0.2:1)
            yticklabels(0:20:100)
            axis([buf 1 buf 1])
            makefigpretty
            
            DT = delaunayTriangulation(PPV(igrp,contains(vartype2,'scEntropy'))',...
                SENS(igrp,contains(vartype2,'scEntropy'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),CXclr,'facealpha',0.15,'edgealpha',0)
            
            DT = delaunayTriangulation(PPV(igrp,contains(vartype2,'fcEntropy'))',...
                SENS(igrp,contains(vartype2,'fcEntropy'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),MIclr,'facealpha',0.15,'edgealpha',0)
            
            DT = delaunayTriangulation(PPV(igrp,contains(vartype2,'scSpectralA'))',...
                SENS(igrp,contains(vartype2,'scSpectralA'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),ASclr,'facealpha',0.15,'edgealpha',0)
            
            DT = delaunayTriangulation(PPV(igrp,contains(vartype2,'scSpectralR'))',...
                SENS(igrp,contains(vartype2,'scSpectralR'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),RSclr,'facealpha',0.15,'edgealpha',0)
            
            DT = delaunayTriangulation(PPV(igrp,contains(vartype2,'fcSpectral'))',...
                SENS(igrp,contains(vartype2,'fcSpectral'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),PIclr,'facealpha',0.15,'edgealpha',0)
                      
            axis([0 1 0 1])
            axis square
            
            print('-dsvg',sprintf('./figures/PPV_vs_recall_scatter_plot_%s.svg',g{igrp}))
            print('-dpng',sprintf('./figures/PPV_vs_recall_scatter_plot_%s.png',g{igrp}))
        end
        %% Plot precision vs recall upper quadrant
                
        
        ASclr = [0, 0.4470, 0.7410];
        CXclr = [0.8500, 0.3250, 0.0980];
        MIclr = [0.9290, 0.6940, 0.1250];
        RSclr = [0.4940, 0.1840, 0.5560];
        PIclr = [0.4660, 0.6740, 0.1880];
          
        
        hb1 = 320;
        hb2 = 200;
        hb3 = 250;
        
        for igrp = 1:3
            myfigure2
            scatter(PPV(igrp,contains(vartype2,'scEntropy')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'scEntropy')&~idxLMM&~idxPC),hb1,CXclr,'fill','markerfacecolor',CXclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'fcEntropy')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'fcEntropy')&~idxLMM&~idxPC),hb1,MIclr,'fill','markerfacecolor',MIclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'fcSpectral')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'fcSpectral')&~idxLMM&~idxPC),hb1,PIclr,'fill','markerfacecolor',PIclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'scSpectralA')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'scSpectralA')&~idxLMM&~idxPC),hb1,ASclr,'fill','markerfacecolor',ASclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'scSpectralR')&~idxLMM&~idxPC),...
                SENS(igrp,contains(vartype2,'scSpectralR')&~idxLMM&~idxPC),hb1,RSclr,'fill','markerfacecolor',RSclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            
            buf = -0.01;

            % features selected by LMMs
            scatter(PPV(igrp,contains(vartype2,'scEntropy')&idxLMM),SENS(igrp,contains(vartype2,...
            'scEntropy')&idxLMM),hb2,CXclr,'d','fill','markerfacecolor',CXclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'fcEntropy')&idxLMM),SENS(igrp,contains(vartype2,...
                'fcEntropy')&idxLMM),hb2,MIclr,'d','fill','markerfacecolor',MIclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'fcSpectral')&idxLMM),SENS(igrp,contains(vartype2,...
                'fcSpectral')&idxLMM),hb2,PIclr,'d','fill','markerfacecolor',PIclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'scSpectralA')&idxLMM),SENS(igrp,contains(vartype2,...
                'scSpectralA')&idxLMM),hb2,ASclr,'d','fill','markerfacecolor',ASclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'scSpectralR')&idxLMM),SENS(igrp,contains(vartype2,...
                'scSpectralR')&idxLMM),hb2,RSclr,'d','fill','markerfacecolor',RSclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            
            % First principal components
            scatter(PPV(igrp,contains(vartype2,'scEntropy')&idxPC),SENS(igrp,contains(vartype2,...
                'scEntropy')&idxPC),hb3,CXclr,'p','fill','markerfacecolor',CXclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'fcEntropy')&idxPC),SENS(igrp,contains(vartype2,...
                'fcEntropy')&idxPC),hb3,MIclr,'p','fill','markerfacecolor',MIclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'fcSpectral')&idxPC),SENS(igrp,contains(vartype2,...
                'fcSpectral')&idxPC),hb3,PIclr,'p','fill','markerfacecolor',PIclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'scSpectralA')&idxPC),SENS(igrp,contains(vartype2,...
                'scSpectralA')&idxPC),hb3,ASclr,'p','fill','markerfacecolor',ASclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(PPV(igrp,contains(vartype2,'scSpectralR')&idxPC),SENS(igrp,contains(vartype2,...
                'scSpectralR')&idxPC),hb3,RSclr,'p','fill','markerfacecolor',RSclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            
            %legend boxoff
            title(g{igrp},'fontsize',18)
            xlabel('Precision (%)')
            ylabel('Recall (%)')
            xticks(0:0.1:1)
            xticklabels(0:10:100)
            yticks(0:0.1:1)
            yticklabels(0:10:100)
            axis([0.5 1 0.5 1])
            makefigpretty
            
            for ivar = frstvar:length(T.Properties.VariableNames)
                if PPV(igrp,ivar) >= 0.5 && SENS(igrp,ivar) >= 0.5
                    idx = ivar-frstvar+1;  
                    if contains(vartype2{ivar},'scEntropy')
                        text(PPV(igrp,ivar)+jitter(1,ivar),SENS(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',CXclr,'FontSize',16)
                    elseif contains(vartype2{ivar},'fcEntropy')
                        text(PPV(igrp,ivar)+jitter(1,ivar),SENS(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',MIclr,'FontSize',16)
                    elseif contains(vartype2{ivar},'scSpectralA')
                        text(PPV(igrp,ivar)+jitter(1,ivar),SENS(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',ASclr,'FontSize',16)
                    elseif contains(vartype2{ivar},'scSpectralR')
                        text(PPV(igrp,ivar)+jitter(1,ivar),SENS(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',RSclr,'FontSize',16)
                    elseif contains(vartype2{ivar},'fcSpectral')
                        text(PPV(igrp,ivar)+jitter(1,ivar),SENS(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',PIclr,'FontSize',16)
                    else
                        error('Unknown variable type')
                    end
                end
            end
            
            print('-dsvg',sprintf('./Figures/PPV_vs_recall_scatter_plot_upper_quadrant_%s.svg',g{igrp}))
            print('-dpng',sprintf('./Figures/PPV_vs_recall_scatter_plot_upper_quadrant_%s.png',g{igrp}))
        end      

        
        %% Plot sensitivity vs specificity
                
        for igrp = 1:3
            myfigure2
            scatter(SENS(igrp,contains(vartype2,'scEntropy')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'scEntropy')&~idxLMM&~idxPC),270,CXclr,'fill',...
                'markeredgecolor','k','LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'fcEntropy')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'fcEntropy')&~idxLMM&~idxPC),270,MIclr,'fill',...
                'markeredgecolor','k','LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'fcSpectral')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'fcSpectral')&~idxLMM&~idxPC),270,PIclr,'fill',...
                'markeredgecolor','k','LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'scSpectralA')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'scSpectralA')&~idxLMM&~idxPC),270,ASclr,'fill',...
                'markeredgecolor','k','LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'scSpectralR')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'scSpectralR')&~idxLMM&~idxPC),270,RSclr,'fill',...
                'markeredgecolor','k','LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            
            buf = -0.01;
            plot([buf 1],[buf 1],'k--')
            

            % features selected by LMMs
            scatter(SENS(igrp,contains(vartype2,'scEntropy')&idxLMM),SPEC(igrp,contains(vartype2,...
                'scEntropy')&idxLMM),300,'d','fill','markerfacecolor',CXclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'fcEntropy')&idxLMM),SPEC(igrp,contains(vartype2,...
                'fcEntropy')&idxLMM),300,MIclr,'d','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'fcSpectral')&idxLMM),SPEC(igrp,contains(vartype2,...
            'fcSpectral')&idxLMM),300,PIclr,'d','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'scSpectralA')&idxLMM),SPEC(igrp,contains(vartype2,...
                'scSpectralA')&idxLMM),300,ASclr,'d','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'scSpectralR')&idxLMM),SPEC(igrp,contains(vartype2,...
                'scSpectralR')&idxLMM),300,RSclr,'d','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            
            % First principal components
            scatter(SENS(igrp,contains(vartype2,'scEntropy')&idxPC),SPEC(igrp,contains(vartype2,...
                'scEntropy')&idxPC),400,CXclr,'p','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'fcEntropy')&idxPC),SPEC(igrp,contains(vartype2,...
                'fcEntropy')&idxPC),400,MIclr,'p','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'fcSpectral')&idxPC),SPEC(igrp,contains(vartype2,...
                'fcSpectral')&idxPC),400,PIclr,'p','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'scSpectralA')&idxPC),SPEC(igrp,contains(vartype2,...
                'scSpectralA')&idxPC),400,ASclr,'p','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'scSpectralR')&idxPC),SPEC(igrp,contains(vartype2,...
                'scSpectralR')&idxPC),400,RSclr,'p','fill','markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            
            %legend boxoff
            
            title(g{igrp},'fontsize',18)
            xlabel('Sensitivity (%)')
            ylabel('Specificity (%)')
            xticks(0:0.2:1)
            xticklabels(0:20:100)
            yticks(0:0.2:1)
            yticklabels(0:20:100)
            axis([buf 1 buf 1])
            makefigpretty
            
            DT = delaunayTriangulation(SENS(igrp,contains(vartype2,'scEntropy'))',...
                SPEC(igrp,contains(vartype2,'scEntropy'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),CXclr,'facealpha',0.15,'edgealpha',0)
            
            DT = delaunayTriangulation(SENS(igrp,contains(vartype2,'fcEntropy'))',...
                SPEC(igrp,contains(vartype2,'fcEntropy'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),MIclr,'facealpha',0.15,'edgealpha',0)
            
            DT = delaunayTriangulation(SENS(igrp,contains(vartype2,'scSpectralA'))',...
                SPEC(igrp,contains(vartype2,'scSpectralA'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),ASclr,'facealpha',0.15,'edgealpha',0)
            
            DT = delaunayTriangulation(SENS(igrp,contains(vartype2,'scSpectralR'))',...
                SPEC(igrp,contains(vartype2,'scSpectralR'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),RSclr,'facealpha',0.15,'edgealpha',0)
            
            DT = delaunayTriangulation(SENS(igrp,contains(vartype2,'fcSpectral'))',...
                SPEC(igrp,contains(vartype2,'fcSpectral'))');
            k = convexHull(DT);
            fill(DT.Points(k,1),DT.Points(k,2),PIclr,'facealpha',0.15,'edgealpha',0)
                      
            axis([0 1 0 1])
            axis square
            
            print('-dsvg',sprintf('./figures/sens_vs_spec_scatter_plot_%s.svg',g{igrp}))
            print('-dpng',sprintf('./figures/sens_vs_spec_scatter_plot_%s.png',g{igrp}))
        end
        %% Plot sensitivity vs specificity upper quadrant
                
        
        ASclr = [0, 0.4470, 0.7410];
        CXclr = [0.8500, 0.3250, 0.0980];
        MIclr = [0.9290, 0.6940, 0.1250];
        RSclr = [0.4940, 0.1840, 0.5560];
        PIclr = [0.4660, 0.6740, 0.1880];
          
        
        hb1 = 320;
        hb2 = 200;
        hb3 = 250;
        
        for igrp = 1:3
            myfigure2
            scatter(SENS(igrp,contains(vartype2,'scEntropy')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'scEntropy')&~idxLMM&~idxPC),hb1,CXclr,'fill','markerfacecolor',CXclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'fcEntropy')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'fcEntropy')&~idxLMM&~idxPC),hb1,MIclr,'fill','markerfacecolor',MIclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'fcSpectral')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'fcSpectral')&~idxLMM&~idxPC),hb1,PIclr,'fill','markerfacecolor',PIclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'scSpectralA')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'scSpectralA')&~idxLMM&~idxPC),hb1,ASclr,'fill','markerfacecolor',ASclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'scSpectralR')&~idxLMM&~idxPC),...
                SPEC(igrp,contains(vartype2,'scSpectralR')&~idxLMM&~idxPC),hb1,RSclr,'fill','markerfacecolor',RSclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            
            buf = -0.01;

            % features selected by LMMs
            scatter(SENS(igrp,contains(vartype2,'scEntropy')&idxLMM),SPEC(igrp,contains(vartype2,...
            'scEntropy')&idxLMM),hb2,CXclr,'d','fill','markerfacecolor',CXclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'fcEntropy')&idxLMM),SPEC(igrp,contains(vartype2,...
                'fcEntropy')&idxLMM),hb2,MIclr,'d','fill','markerfacecolor',MIclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'fcSpectral')&idxLMM),SPEC(igrp,contains(vartype2,...
                'fcSpectral')&idxLMM),hb2,PIclr,'d','fill','markerfacecolor',PIclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'scSpectralA')&idxLMM),SPEC(igrp,contains(vartype2,...
                'scSpectralA')&idxLMM),hb2,ASclr,'d','fill','markerfacecolor',ASclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'scSpectralR')&idxLMM),SPEC(igrp,contains(vartype2,...
                'scSpectralR')&idxLMM),hb2,RSclr,'d','fill','markerfacecolor',RSclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            
            % First principal components
            scatter(SENS(igrp,contains(vartype2,'scEntropy')&idxPC),SPEC(igrp,contains(vartype2,...
                'scEntropy')&idxPC),hb3,CXclr,'p','fill','markerfacecolor',CXclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'fcEntropy')&idxPC),SPEC(igrp,contains(vartype2,...
                'fcEntropy')&idxPC),hb3,MIclr,'p','fill','markerfacecolor',MIclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'fcSpectral')&idxPC),SPEC(igrp,contains(vartype2,...
                'fcSpectral')&idxPC),hb3,PIclr,'p','fill','markerfacecolor',PIclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'scSpectralA')&idxPC),SPEC(igrp,contains(vartype2,...
                'scSpectralA')&idxPC),hb3,ASclr,'p','fill','markerfacecolor',ASclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            scatter(SENS(igrp,contains(vartype2,'scSpectralR')&idxPC),SPEC(igrp,contains(vartype2,...
                'scSpectralR')&idxPC),hb3,RSclr,'p','fill','markerfacecolor',RSclr,'markeredgecolor','k',...
                'LineWidth',1.5,'markerfacealpha',0.5,'markeredgealpha',0.5)
            
            %legend boxoff
            title(g{igrp},'fontsize',18)
            xlabel('Sensitivity (%)')
            ylabel('Specificity (%)')
            xticks(0:0.1:1)
            xticklabels(0:10:100)
            yticks(0:0.1:1)
            yticklabels(0:10:100)
            axis([0.5 1 0.5 1])
            makefigpretty
            
            for ivar = frstvar:length(T.Properties.VariableNames)
                if SENS(igrp,ivar) >= 0.5 && SPEC(igrp,ivar) >= 0.5
                    idx = ivar-frstvar+1;  
                    if contains(vartype2{ivar},'scEntropy')
                        text(SENS(igrp,ivar)+jitter(1,ivar),SPEC(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',CXclr,'FontSize',16)
                    elseif contains(vartype2{ivar},'fcEntropy')
                        text(SENS(igrp,ivar)+jitter(1,ivar),SPEC(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',MIclr,'FontSize',16)
                    elseif contains(vartype2{ivar},'scSpectralA')
                        text(SENS(igrp,ivar)+jitter(1,ivar),SPEC(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',ASclr,'FontSize',16)
                    elseif contains(vartype2{ivar},'scSpectralR')
                        text(SENS(igrp,ivar)+jitter(1,ivar),SPEC(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',RSclr,'FontSize',16)
                    elseif contains(vartype2{ivar},'fcSpectral')
                        text(SENS(igrp,ivar)+jitter(1,ivar),SPEC(igrp,ivar)+jitter(2,ivar)+0.01,varnames3{idx},'color',PIclr,'FontSize',16)
                    else
                        error('Unknown variable type')
                    end
                end
            end
                                
            % vertical and horizontal lines for AUC thresholds
            tcks = 0:0.1:1;

            for i = 1:length(tcks)      
                plot(ones(1,100).*tcks(i),linspace(buf,1,100),'k:')
                plot(linspace(buf,1,100),ones(1,100).*tcks(i),'k:')
            end
            
            print('-dsvg',sprintf('./Figures/sens_vs_spec_scatter_plot_upper_quadrant_%s.svg',g{igrp}))
            print('-dpng',sprintf('./Figures/sens_vs_spec_scatter_plot_upper_quadrant_%s.png',g{igrp}))
        end      
        
        
        %% Save the output
        switch method
            case 'PCA'
                save EEG_ML_vars_PCA_MWU
            case 'LMM'
                save EEG_ML_vars_LMM_MWU
        end
        

end

% The End :) 
