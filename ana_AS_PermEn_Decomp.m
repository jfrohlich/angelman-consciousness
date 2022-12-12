%%% STATS 2021 analysis %%%
% Joel Frohlich
% University of California, Los Angeles (UCLA)
% Monti Lab, Psychology Department
%
% University of Tuebingen, Germany
% Institue for Neuromodulation and Neurotechnology 
%
% Last update: 12 Dec, 2022 (cleaned up code and comments)

% For each subject, the code performs a permutation entropy decomposition
% using the PermEnDecomposition function, and saves the results to a
% variable. The code then saves the results to a file using the save
% function, and repeats this process for the remaining subjects. Finally,
% the code loads the saved results and performs some additional analysis to
% test for differences between spectral and phase components of the data.


% Optional, clear all variables from workspace
clearvars

rng(45783) % seed added 02.26.19

OS = computer; % detect opperating system; disclaimer: the computer function only returns a limited number of values, and it is not guaranteed to return the correct value on all systems
dbstop if error
new_fsample = 125; % target sampling rate, see King et al. 2013 supplement
ft_defaults
nchan = 19;

switch OS
    case 'GLNXA64'
        pth = './postprocessed/AS/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
    case 'MACI64'
        pth = './postprocessed/AS/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
    case 'PCWIN64'
        pth = './postprocessed/AS/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
end


% these are the files we already know we want to use
% load('select_files','select_files')

sleep = true; % asleep or awake?
[sleepsec,sleep_subjects,sleep_ages] = get_sections(sleep,files,pth,new_fsample);

sleep = false; % asleep or awake?
[wakesec,wake_subjects,wake_ages] = get_sections(sleep,files,pth,new_fsample);

save DecompSubs sleep_subjects sleep_ages wake_subjects wake_ages
%%
assert(length(sleep_subjects)==length(wake_subjects),'Different number of wake and sleep subjects')

for isub = 1:length(sleep_subjects)
    assert(strcmp(sleep_subjects{isub},wake_subjects{isub}),'Subjects don''t match between conditions')
end

assert(length(wakesec)==length(sleepsec),'Different number of sleep and wake data');
all_decomps = cell(1,length(wakesec));
Nsur = 1000; % use this many surrogates
Nsegs = 20; % use only subjects with at least this many good segments
tic
%
start_ndx = 1;
for isub = length(wakesec)-2:-1:1 % loop backwards
    if length(wakesec{isub}) >= Nsegs && length(sleepsec{isub}) >= Nsegs
        W = wakesec{isub}(randperm(length(wakesec{isub}),Nsegs));
        S = sleepsec{isub}(randperm(length(sleepsec{isub}),Nsegs));
        OUT = PermEnDecomposition(W,wake_subjects{isub},S,sleep_subjects{isub},Nsur,new_fsample);
        all_decomps{isub} = OUT;
        fprintf('Finished subject %s\n',sleep_subjects{isub})
        toc
    else
        fprintf('     Skipping this one, not enough good segments ...\n')
        fprintf('%i%% overall complete, %i minutes elapsed\n',round(isub/length(wakesec)*100),round(toc/60))
    end
    save PermEnDecompAS_Ns1000 all_decomps Nsur Nsegs
end

save PermEnDecompAS_Ns1000 all_decomps Nsur Nsegs sleep_subjects wake_subjects

%% Test for differences between spectral and phase components 

clearvars
Nsur = 1000;
nchan = 19;
load DecompSubs
load(sprintf('PermEnDecompAS_Ns%i.mat',Nsur))

% don't use these subjects
bl = {'100577'};

P = nan(length(all_decomps),length(all_decomps{1}.lags),nchan);
S = nan(length(all_decomps),length(all_decomps{1}.lags),nchan);
X = nan(length(all_decomps),length(all_decomps{1}.lags),nchan);
L = nan(length(all_decomps),length(all_decomps{1}.lags));


rmv = []; % remove these indices

for isub = 1:length(all_decomps)
    if ~isempty(all_decomps{isub}) && ~strcmp(sleep_subjects{isub},bl)
        P(isub,:,:) = all_decomps{isub}.phase;
        S(isub,:,:) = all_decomps{isub}.spec;
        X(isub,:,:) = all_decomps{isub}.interact;
        L(isub,:,:) = all_decomps{isub}.lags;
    else
        rmv = [rmv isub];
    end
end

sleep_subjects(rmv)
sleep_ages(rmv)


P(rmv,:,:) = [];
S(rmv,:,:) = [];
X(rmv,:,:) = [];
L(rmv,:,:) = [];
sleep_subjects(rmv) = [];
wake_subjects(rmv) = [];

PX = P + X; % add phase + interaction

assert(length(sleep_subjects)==length(wake_subjects))
for i = 1:length(sleep_subjects)
    assert(strcmp(sleep_subjects{i},wake_subjects{i}))
end

alpha = 0.01;
assert(all(diff(L)==0,[1 2]),'Different subjects have different lags')
pval = nan(size(P,2), 3);
tval = nan(size(P,2),3);
cohd = nan(size(P,2),3);

% average across channels
Pall = P; % copy
Sall = S; % copy
Xall = X; % copy
PXall = PX; % copy
P = mean(P,3); % average across channels
S = mean(S,3); % average across channels
X = mean(X,3); % average across channels
PX = mean(PX,3); % average across channels

LAGS = L(1,:);

pstats = nan([size(Pall,2) size(Pall,3) 3]);
tstats = nan([size(Pall,2) size(Pall,3) 3]);

for itau = 1:size(P,2)

    PartPSI = [repmat(categorical({'Phase'}),size(P,1),1); ...
        repmat(categorical({'Spectral'}),size(S,1),1); ...
        repmat(categorical({'Interaction'}),size(X,1),1)];
    PartPS = [repmat(categorical({'Phase'}),size(P,1),1); ...
        repmat(categorical({'Spectral'}),size(S,1),1)];
    PartSI = [repmat(categorical({'Spectral'}),size(S,1),1); ...
        repmat(categorical({'Interaction'}),size(X,1),1)];
    PartPI = [repmat(categorical({'Phase'}),size(P,1),1); ...
        repmat(categorical({'Interaction'}),size(X,1),1)];
    PartSPX = [repmat(categorical({'Spectral'}),size(P,1),1); ...
        repmat(categorical({'Non-spectral'}),size(PX,1),1)];
    
    TdPSI = table([P(:,itau); S(:,itau); X(:,itau)],PartPSI,categorical([sleep_subjects sleep_subjects sleep_subjects])',...
        'VariableNames',{'Contribution','Part','Subject'});
    TdPS = table([P(:,itau); S(:,itau)],PartPS,categorical([sleep_subjects sleep_subjects])',...
        'VariableNames',{'Contribution','Part','Subject'});
    TdSI = table([S(:,itau); X(:,itau)],PartSI,categorical([sleep_subjects sleep_subjects])',...
        'VariableNames',{'Contribution','Part','Subject'});
    TdPI = table([P(:,itau); X(:,itau)],PartPI,categorical([sleep_subjects sleep_subjects])',...
        'VariableNames',{'Contribution','Part','Subject'});
    TdSPX = table([P(:,itau); PX(:,itau)],PartSPX,categorical([sleep_subjects wake_subjects])',...
        'VariableNames',{'Contribution','Part','Subject'});
    
    result    = fitlme(TdPSI,'Contribution ~ Part + (1|Subject)');
    resultPS  = fitlme(TdPS,'Contribution ~ Part + (1|Subject)');
    resultSI  = fitlme(TdSI,'Contribution ~ Part + (1|Subject)');
    resultPI  = fitlme(TdPI,'Contribution ~ Part + (1|Subject)');
    resultSPX = fitlme(TdSPX,'Contribution ~ Part + (1|Subject)');
    
    AOV = anova(result);
    comparisons = {'phase vs amplitude','amplitude vs interaction',...
        'phase vs interaction','amplitude vs non-amplitude'};
    
    pval(itau,1) = resultPS.Coefficients.pValue(end);
    pval(itau,2) = resultSI.Coefficients.pValue(end);
    pval(itau,3) = resultPI.Coefficients.pValue(end);
    pval(itau,4) = resultSPX.Coefficients.pValue(end);
    
    tval(itau,1) = resultPS.Coefficients.tStat(end);
    tval(itau,2) = resultSI.Coefficients.tStat(end);
    tval(itau,3) = resultPI.Coefficients.tStat(end);
    tval(itau,4) = resultSPX.Coefficients.tStat(end);
    
    cohd(itau,1) = cohen_d(TdPS{TdPS.Part==categorical(cellstr(repmat('Phase',size(TdPS,1),1))),1},...
        TdPS{TdPS.Part==categorical(cellstr(repmat('Spectral',size(TdPS,1),1))),1});
    cohd(itau,2) = cohen_d(TdSI{TdSI.Part==categorical(cellstr(repmat('Spectral',size(TdSI,1),1))),1},...
        TdSI{TdSI.Part==categorical(cellstr(repmat('Interaction',size(TdSI,1),1))),1});
    cohd(itau,3) = cohen_d(TdPI{TdPI.Part==categorical(cellstr(repmat('Phase',size(TdPI,1),1))),1},...
        TdPI{TdPI.Part==categorical(cellstr(repmat('Interaction',size(TdPI,1),1))),1});
    cohd(itau,4) = cohen_d(TdSPX{TdSPX.Part==categorical(cellstr(repmat('Spectral',size(TdSPX,1),1))),1},...
        TdSPX{TdSPX.Part==categorical(cellstr(repmat('Non-spectral',size(TdSPX,1),1))),1});
    
    cl = 0.04; % axis limit for topoplots
    
    myfigure
    plot_topo_AS_thick(squeeze(mean(Sall(:,itau,:))))
    caxis([-cl cl])
    %mycolorbar
    %title({sprintf('Spectral component, tau = %i',LAGS(itau)),''})
    print('-dpng',sprintf('./Figures/Spectral_comp_Tau%i.png',LAGS(itau)))

    myfigure
    plot_topo_AS_thick(squeeze(mean(PXall(:,itau,:))))
    caxis([-cl cl])
    %mycolorbar
    %title({sprintf('Non-spectral component, tau = %i',LAGS(itau)),''})
    print('-dpng',sprintf('./Figures/Nonspectral_comp_Tau%i.png',LAGS(itau)))
    %pause(0.1)
    %close all
    
end

%% Export table with data for Figure 8

load joel_lay.mat lay
AmplitudeComponent = reshape(Sall,size(Sall,1),size(Sall,2)*size(Sall,3))';
NonamplitudeComponent = reshape(PXall,size(PXall,1),size(PXall,2)*size(PXall,3))';
rows = cell(size(Sall,2)*size(Sall,3),1);
lags = 2.^[3:7];

rcnt = 0;
for ich = 1:size(Sall,3)
    for itau = 1:size(Sall,2)
        rcnt = rcnt + 1;
        rows{rcnt} = sprintf('%s_tau=%ims',lay.label{ich},lags(itau));
    end
end

ColHeaders = [{'Channel_and_lag'} sleep_subjects sleep_subjects];


Table_Fig8_Amp = table(rows,AmplitudeComponent);
Table_Fig8_Nonamp = table(rows,NonamplitudeComponent);

% Export these tables and combine them later in Excel
writetable(Table_Fig8_Amp,'./Table_Fig8_Amp.csv')
writetable(Table_Fig8_Nonamp,'./Table_Fig8_Nonamp.csv')
Theaders = table(ColHeaders);
writetable(Theaders,'./Table_Fig8_headers.csv')

%% Report stats
qval = mafdr(pval(:),'bhfdr',true);
qval = reshape(qval,5,4);

for itau = 1:5
    for icomp = 1:4
        if qval(itau,icomp) < 0.0005
            fprintf('Significant effect *** of %s for tau = %i, p = %1.2d (FDR corrected), t = %2.2f, d = %2.3f\n',...
                comparisons{icomp},LAGS(itau),qval(itau,icomp),tval(itau,icomp),cohd(itau,icomp))
        elseif qval(itau,icomp) < 0.005
            fprintf('Significant effect ** of %s for tau = %i, p = %1.2d (FDR corrected), t = %2.2f, d = %2.3f\n',...
                comparisons{icomp},LAGS(itau),qval(itau,icomp),tval(itau,icomp),cohd(itau,icomp))
        elseif qval(itau,icomp) < 0.05
            fprintf('Significant effect * of %s for tau = %i, p = %1.2d (FDR corrected), t = %2.2f, d = %2.3f\n',...
                comparisons{icomp},LAGS(itau),qval(itau,icomp),tval(itau,icomp),cohd(itau,icomp))
        end
            
    end
end

comparison = [repmat(comparisons(1),5,1); repmat(comparisons(2),5,1); repmat(comparisons(3),5,1); repmat(comparisons(4),5,1)];
taucol = repmat([8 16 32 64 128]',4,1);

% Generate supplemental table
Tdecomp = table(comparison,taucol,qval(:),0-tval(:),abs(cohd(:)),'VariableNames',...
    {'Comparison','Tau (ms)','p-value (FDR corrected)','t-stat','Cohen''s d'})

writetable(Tdecomp,'PermEnDecomp.csv')


%% Figure

qval = qval(:,end); % just keep last column

IL = [S;PX];
IL = reshape(IL,size(IL,1)/2,size(IL,2)*2);

myfigure2
myviolin(IL,'paired',false,'bw',0.01)
xticks(1:size(IL,2))
yticks(-0.2:0.1:0.2)
plot(linspace(0,size(IL,2)+1,100),zeros(1,100),'k:')
xlim([0 size(IL,2)+1])
makefighandsome
% Mark significant results with an asterisk
for itau = 1:size(qval,1)
    for icomp = 1:size(qval,2)
        where = (itau-1)*2+icomp;
        if pval(itau,icomp) < alpha
            text(where,0.1,sprintf('p = %1.1d *',qval(itau,icomp)),'fontsize',11)
        else
            %text(where,0.1,sprintf('p = %1.2f',pval(itau,icomp)),'fontsize',10)
        end
    end
end

% Uncomment below for Latex
% ylabel('Contribution','Interpreter','latex')
% set(gca,'XTickLabel',{sprintf('Phase, %stau$ = %i ms','$\',L(1,1)),sprintf('Power, %stau$ = %i ms','$\',L(1,1)),...
%     sprintf('Phase, %stau$ = %i ms','$\',L(1,2)),sprintf('Power, %stau$ = %i ms','$\',L(1,2)),...
%     sprintf('Phase, %stau$ = %i ms','$\',L(1,3)),sprintf('Power, %stau$ = %i ms','$\',L(1,3)),...
%     sprintf('Phase, %stau$ = %i ms','$\',L(1,4)),sprintf('Power, %stau$ = %i ms','$\',L(1,4)),...
%     sprintf('Phase, %stau$ = %i ms','$\',L(1,5)),sprintf('Power, %stau$ = %i ms','$\',L(1,5))},...
%     'TickLabelInterpreter', 'latex')
% xtickangle(45)
% title('Permutation entropy decomposition, channel Pz','fontsize',22,'Interpreter','latex')

ylabel('\Delta PermEn')
ylim([-0.2 0.2])
xticklabels({'Amplitude, 16-40 Hz','Non-amplitude, 16-40 Hz'...
    'Amplitude, 8-20 Hz','Non-amplitude, 8-20 Hz',...
    'Amplitude, 4-10 Hz','Non-amplitude, 4-10 Hz',...
    'Amplitude, 2-5 Hz','Non-amplitude, 2-5 Hz',...
    'Amplitude, 1-2.5 Hz','Non-amplitude, 1-2.5 Hz'})
% xticklabels({sprintf('Phase, tau = %i ms',L(1,1)),sprintf('Power, tau = %i ms',L(1,1)),...
%     sprintf('Phase, tau = %i ms',L(1,2)),sprintf('Power, tau = %i ms',L(1,2)),...
%     sprintf('Phase, tau = %i ms',L(1,3)),sprintf('Power, tau = %i ms',L(1,3)),...
%     sprintf('Phase, tau = %i ms',L(1,4)),sprintf('Power, tau = %i ms',L(1,4)),...
%     sprintf('Phase, tau = %i ms',L(1,5)),sprintf('Power, tau = %i ms',L(1,5))})
xtickangle(45)

title('Permutation entropy decomposition, channel average','fontsize',22)
print('-dpng','./Figures/PermEnDecomp.png')
print('-dsvg','./Figures/PermEnDecomp.svg')


%%

function[all_sections,subs,ages] = get_sections(sleep,files,pth,new_fsample)

ntrl = 1; % only one trial for AS data

all_sections = cell(1,length(files));
subs = cell(1,length(files));
ages = cell(1,length(files));
Tages = readtable('./AS_table_2019_JF.csv');

tic
for ifile = 1:length(files)
    %if logical(sum(ismember(select_files,files(ifile).name(1:end-4))))
    load(sprintf('%s%s',pth,files(ifile).name),'data')
    fprintf('Now loading %s \n',files(ifile).name)
    subs{ifile} = files(ifile).name(1:6);
    
    us = find(data.cfg.fileout=='_'); % find underscores
    slsh = find(data.cfg.fileout=='/'); % find slashes
    fd = str2double(data.cfg.fileout(us(1)+1:us(2)-1)); % file date
    if ~isfield(data,'ID')
        if isempty(slsh)
            data.ID = data.cfg.fileout(1:us(1)-1);
        else
            data.ID = data.cfg.fileout(slsh(end)+1:us(1)-1);
        end
    end
    try
        assert(isa(data.ID,'double'))
        rw = find(data.ID == Tages.ID & fd == Tages.TimeRecording);
    catch
        rw = find(str2double(data.ID) == Tages.ID & fd == Tages.TimeRecording);
    end
    ages{ifile} = Tages.AgeRecording(rw);
    
    % check if this file had useable sleep
    
    % Commented out because second trial is always sleep
    if ~isfield(data.cfg.dattype,'SLEEP_MONTI') || isempty(data.cfg.dattype.SLEEP_MONTI)
        fprintf('Skipping this one, no useable sleep ... \n')
        continue
    end
    
    %Get datamat var for freq transform
    datamat = data.trial{ntrl};
    
    % re-average-ref datamat (after interpolation and ICA)
    datamat = datamat - nanmean(datamat,1);
    
    % Do some santiy checks
    assert(size(datamat,1) == 19,'channel(s) missing!')
    % allow for small number of NaNs at beginning and end of file
    assert(sum(isnan(datamat(1:end)))/numel(datamat) < 0.01, 'Many NaNs in data!')
    
    switch sleep
        case true
            %% convert awake data to NaNs
            keep = data.cfg.dattype.SLEEP_MONTI; % Data to retain for analysis (NOT converted to NaNs!)
            datamat2 = ones(size(datamat)); % ones indicate data that will be switched to NaNs
            
            for islp = 1:size(keep,1)
                idpt = keep(islp,1) - data.sampleinfo(1,1) + 1;
                datamat2(:,idpt:idpt+diff(keep(islp,:))) = 0;
            end
            
            wkedx = datamat2 == 1;
            
            datamat(wkedx) = NaN; % replace awake data with NaNs
            badart = [];
            
        case false
            % exclude the Monti sleep labels, plus also the "drowsy" sleep
            % labels (these two add together to the original sleep labels
            % from the Roche analysis)
            badart = data.cfg.dattype.SLEEP_MONTI; 
    end
    
    % Set bad data to NaNs for purpose of taking freq transform
    badart = [badart; data.cfg.dattype.bad; data.cfg.dattype.flash; data.cfg.dattype.drowsy];
    
    % add muscle to artifacts if it is a field
    if isfield(data.cfg.dattype,'muscle')
        badart = [badart; data.cfg.dattype.muscle];
    end
    
    badart(badart==0) = 1; % do this so no 0 indicies that give bugs
    
    for iart = 1:size(badart,1)
        idpt = badart(iart,1) - data.sampleinfo(ntrl,1) + 1;
        if idpt > 0 && idpt+diff(badart(iart,:)) <= size(datamat,2)
            datamat(:,idpt:idpt+diff(badart(iart,:))) = nan; 
        end
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
    
    
    % downsample to 125 Hz prior to wSMI, this sampling rate is justified
    % by the supplment for King et al. 2013
    datamat = single(resample(double(datamat)',new_fsample,data.fsample))';
    
    %% partion into n sec windows and compute wSMI
    
    nwin = 5;
    win_len = new_fsample*nwin; % n sec windows
    
    gcount = 0;
    
    % sections for decomp
    sections = cell(1,1);
 
    % windowing parameters
    % 1.0 --> no overlap; 0.5 --> 50 percent
    window_shift=0.5; % do 50% overlap with the matched data since it's less data
    n_shift = round(win_len*window_shift);
    for isection = 1:n_shift:size(datamat,2)-win_len+1
        section = double(datamat(:,isection:isection+win_len-1));
               
        if ~any(isnan(section(:))) % IF there are no NaNs in the window
            gcount = gcount + 1;
            sections{gcount} = section;
        end
    end
    all_sections{ifile} = sections;
    
end


end


