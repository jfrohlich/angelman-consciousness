% This Matlab code uses the FOOOF toolbox to analyze the power spectra of 
% the data in dataPeaks. It loads the data and the necessary toolbox and 
% sets some parameters for the analysis. Then, it creates some empty arrays 
% to store the results of the analysis.
% 
% The code then iterates over all of the subjects in the data and, for each 
% subject, extracts the power spectrum of the data and interpolates it onto 
% a new frequency range, foihd2. It then runs the fooof function on the 
% interpolated data to find the peak frequency and bandwidth of the 
% spectrum and saves the results to the appropriate arrays. Finally, it 
% saves the results to a file.

%     The code begins by loading the dataPeaks file, which contains power spectrum data for the sleep and wake data of each subject in the 3 different groups.
%     The settings variable is then defined as an empty array, and variables for the frequency range and frequencies for the power spectrum analysis are defined.
%     The code then creates empty arrays to store the peak frequency data for each subject in each group.
%     The code then loops through each subject in the TD group, performs the power spectrum analysis on their wake data, and finds the peak frequency. The peak frequency is then added to the TDwakePeaksExact array and to the appropriate bin in the TDwakePeaks histogram.
%     This process is repeated for the AS and DS groups with their respective arrays.
%     Next, the code plots the histograms of peak frequencies for the 3 groups and saves the figures to files.
%     Finally, the code saves the peak frequency data for each subject in each group to a .mat file.

load dataPeaks

try 
    pyenv("ExecutionMode","OutOfProcess")
catch
    fprintf('Python already added\n')
end

settings = [];
%settings.aperiodic_mode = 'knee';
frange = [0.5 32];
foihd2 = linspace(0.5,32,200);

% Note: power vectors must be linearly scaled and linearly spaced for FOOOF

N = size(TDWakePowCA,2);
assert(size(TDWakePowCA,2)==size(TDSleepPowCA,2))
TDwakePeaks = zeros(N,6);
TDsleepPeaks = zeros(N,6);
TDwakePeaksExact = zeros(N,1);
TDsleepPeaksExact = zeros(N,1);

N = size(ASWakePowCA,2);
assert(size(ASWakePowCA,2)==size(ASSleepPowCA,2))
ASwakePeaks = zeros(N,6);
ASsleepPeaks = zeros(N,6);
ASwakePeaksExact = zeros(N,1);
ASsleepPeaksExact = zeros(N,1);

N = size(DSWakePowCA,2);
assert(size(DSWakePowCA,2)==size(DSSleepPowCA,2))
DSwakePeaks = zeros(N,6);
DSsleepPeaks = zeros(N,6);
DSwakePeaksExact = zeros(N,1);
DSsleepPeaksExact = zeros(N,1);

% TD wake
for isub = 1:size(TDWakePowCA,2)
    fprintf('TD wake subject %i\n',isub)
    tmp = 10.^interp1(foi',TDWakePowCA(:,isub),foihd2);
    if any(isnan(tmp))% spectra with nans were rejected due to number of good windows
        fprintf('All nans, skipping this on ...\n')
        continue
    end
    out = fooof(foihd2,tmp,frange,settings);
    if isempty(out.peak_params)
        fprintf('No peak detected!\n')
        continue
    end
    [~,ifreq] = max(out.peak_params(:,2));
    f = out.peak_params(ifreq,1);
    fprintf('     Peak at %1.3f Hz\n',f)
    fprintf('          Peak power %1.3f\n',out.peak_params(ifreq,2))
    fprintf('          Peak bandwidth %1.3f\n',out.peak_params(ifreq,3))
    TDwakePeaksExact(isub) = f;
    if f < 1
        TDwakePeaks(isub,1) = TDwakePeaks(isub,1) + 1;
    elseif f < 2
        TDwakePeaks(isub,2) = TDwakePeaks(isub,2) + 1;
    elseif f < 4
        TDwakePeaks(isub,3) = TDwakePeaks(isub,3) + 1;
    elseif f < 8
        TDwakePeaks(isub,4) = TDwakePeaks(isub,4) + 1;
    elseif f < 16
        TDwakePeaks(isub,5) = TDwakePeaks(isub,5) + 1;
    elseif f < 32
        TDwakePeaks(isub,6) = TDwakePeaks(isub,6) + 1;
    else
        error('Unrecognized frequency')
    end
end

% AS wake
for isub = 1:size(ASWakePowCA,2)
    fprintf('AS wake subject %i\n',isub)
    tmp = 10.^interp1(foi,ASWakePowCA(:,isub),foihd2);
    if any(isnan(tmp))% spectra with nans were rejected due to number of good windows
        fprintf('All nans, skipping this on ...\n')
        continue
    end
    out = fooof(foihd2,tmp,frange,settings);
    if isempty(out.peak_params)
        fprintf('No peak detected!\n')
        continue
    end
    [~,ifreq] = max(out.peak_params(:,2));
    
    f = out.peak_params(ifreq,1);
    fprintf('     Peak at %1.3f Hz\n',f)
    fprintf('          Peak power %1.3f\n',out.peak_params(ifreq,2))
    fprintf('          Peak bandwidth %1.3f\n',out.peak_params(ifreq,3))
    ASwakePeaksExact(isub) = f;
    if f < 1
        ASwakePeaks(isub,1) = ASwakePeaks(isub,1) + 1;
    elseif f < 2
        ASwakePeaks(isub,2) = ASwakePeaks(isub,2) + 1;
    elseif f < 4
        ASwakePeaks(isub,3) = ASwakePeaks(isub,3) + 1;
    elseif f < 8
        ASwakePeaks(isub,4) = ASwakePeaks(isub,4) + 1;
    elseif f < 16
        ASwakePeaks(isub,5) = ASwakePeaks(isub,5) + 1;
    elseif f < 32
        ASwakePeaks(isub,6) = ASwakePeaks(isub,6) + 1;
    else
        error('Unrecognized frequency')
    end
end

%
% DS wake
for isub = 1:size(DSWakePowCA,2)
    fprintf('DS wake subject %i\n',isub)
    tmp = 10.^interp1(foi,DSWakePowCA(:,isub),foihd2);
    if any(isnan(tmp))% spectra with nans were rejected due to number of good windows
        fprintf('All nans, skipping this on ...\n')
        continue
    end
    out = fooof(foihd2,tmp,frange,settings);
    if isempty(out.peak_params)
        fprintf('No peak detected!\n')
        continue
    end
    [~,ifreq] = max(out.peak_params(:,2));
    
    f = out.peak_params(ifreq,1);
    fprintf('     Peak at %1.3f Hz\n',f)
    fprintf('          Peak power %1.3f\n',out.peak_params(ifreq,2))
    fprintf('          Peak bandwidth %1.3f\n',out.peak_params(ifreq,3))
    DSwakePeaksExact(isub) = f;
    if f < 1
        DSwakePeaks(isub,1) = DSwakePeaks(isub,1) + 1;
    elseif f < 2
        DSwakePeaks(isub,2) = DSwakePeaks(isub,2) + 1;
    elseif f < 4
        DSwakePeaks(isub,3) = DSwakePeaks(isub,3) + 1;
    elseif f < 8
        DSwakePeaks(isub,4) = DSwakePeaks(isub,4) + 1;
    elseif f < 16
        DSwakePeaks(isub,5) = DSwakePeaks(isub,5) + 1;
    elseif f < 32
        DSwakePeaks(isub,6) = DSwakePeaks(isub,6) + 1;
    else
        error('Unrecognized frequency')
    end
end

% TD sleep
for isub = 1:size(TDSleepPowCA,2)
    fprintf('TD sleep subject %i\n',isub)
    tmp = 10.^interp1(foi,TDSleepPowCA(:,isub),foihd2);
    if any(isnan(tmp))% spectra with nans were rejected due to number of good windows
        fprintf('All nans, skipping this on ...\n')
        continue
    end
    out = fooof(foihd2,tmp,frange,settings);
    if isempty(out.peak_params)
        fprintf('No peak detected!\n')
        continue
    end
    [~,ifreq] = max(out.peak_params(:,2));
    
    f = out.peak_params(ifreq,1);
    TDsleepPeaksExact(isub) = f;
    fprintf('     Peak at %1.3f Hz\n',f)
    fprintf('          Peak power %1.3f\n',out.peak_params(ifreq,2))
    fprintf('          Peak bandwidth %1.3f\n',out.peak_params(ifreq,3))
    if f < 1
        TDsleepPeaks(isub,1) = TDsleepPeaks(isub,1) + 1;
    elseif f < 2
        TDsleepPeaks(isub,2) = TDsleepPeaks(isub,2) + 1;
    elseif f < 4
        TDsleepPeaks(isub,3) = TDsleepPeaks(isub,3) + 1;
    elseif f < 8
        TDsleepPeaks(isub,4) = TDsleepPeaks(isub,4) + 1;
    elseif f < 16
        TDsleepPeaks(isub,5) = TDsleepPeaks(isub,5) + 1;
    elseif f < 32
        TDsleepPeaks(isub,6) = TDsleepPeaks(isub,6) + 1;
    else
        error('Unrecognized frequency')
    end
end
%
% AS sleep
for isub = 1:size(ASSleepPowCA,2)
    fprintf('AS sleep subject %i\n',isub)
    tmp = 10.^interp1(foi,ASSleepPowCA(:,isub),foihd2);
    if any(isnan(tmp))% spectra with nans were rejected due to number of good windows
        fprintf('All nans, skipping this on ...\n')
        continue
    end
    out = fooof(foihd2,tmp,frange,settings);
    if isempty(out.peak_params)
        fprintf('No peak detected!\n')
        continue
    end
    [~,ifreq] = max(out.peak_params(:,2));
    % don't allow peaks with a height < 1 and bandwidth >= 2
    %if ~(out.peak_params(ifreq,2) < 1 && out.peak_params(ifreq,3) >= 2)
    f = out.peak_params(ifreq,1);
    ASsleepPeaksExact(isub) = f;
    fprintf('     Peak at %1.3f Hz\n',f)
    fprintf('          Peak power %1.3f\n',out.peak_params(ifreq,2))
    fprintf('          Peak bandwidth %1.3f\n',out.peak_params(ifreq,3))
    if f < 1
        ASsleepPeaks(isub,1) = ASsleepPeaks(isub,1) + 1;
    elseif f < 2
        ASsleepPeaks(isub,2) = ASsleepPeaks(isub,2) + 1;
    elseif f < 4
        ASsleepPeaks(isub,3) = ASsleepPeaks(isub,3) + 1;
    elseif f < 8
        ASsleepPeaks(isub,4) = ASsleepPeaks(isub,4) + 1;
    elseif f < 16
        ASsleepPeaks(isub,5) = ASsleepPeaks(isub,5) + 1;
    elseif f < 32
        ASsleepPeaks(isub,6) = ASsleepPeaks(isub,6) + 1;
    else
        error('Unrecognized frequency')
    end
end


% DS sleep
for isub = 1:size(DSSleepPowCA,2)
    fprintf('DS sleep subject %i\n',isub)
    tmp = 10.^interp1(foi,DSSleepPowCA(:,isub),foihd2);
    if any(isnan(tmp))% spectra with nans were rejected due to number of good windows
        fprintf('All nans, skipping this on ...\n')
        continue
    end
    out = fooof(foihd2,tmp,frange,settings);
    if isempty(out.peak_params)
        fprintf('No peak detected!\n')
        continue
    end
    [~,ifreq] = max(out.peak_params(:,2));
    f = out.peak_params(ifreq,1);
    DSsleepPeaksExact(isub) = f;
    fprintf('     Peak at %1.3f Hz\n',f)
    fprintf('          Peak power %1.3f\n',out.peak_params(ifreq,2))
    fprintf('          Peak bandwidth %1.3f\n',out.peak_params(ifreq,3))
    if f < 1
        DSsleepPeaks(isub,1) = DSsleepPeaks(isub,1) + 1;
    elseif f < 2
        DSsleepPeaks(isub,2) = DSsleepPeaks(isub,2) + 1;
    elseif f < 4
        DSsleepPeaks(isub,3) = DSsleepPeaks(isub,3) + 1;
    elseif f < 8
        DSsleepPeaks(isub,4) = DSsleepPeaks(isub,4) + 1;
    elseif f < 16
        DSsleepPeaks(isub,5) = DSsleepPeaks(isub,5) + 1;
    elseif f < 32
        DSsleepPeaks(isub,6) = DSsleepPeaks(isub,6) + 1;
    else
        error('Unrecognized frequency')
    end
end

save('AllPeakData','ASwakePeaks','ASsleepPeaks','TDwakePeaks','TDsleepPeaks',...
    'DSwakePeaks','DSsleepPeaks')

%% Make histograms

bins = 15:30:345;

%% TD
Wake = sum(TDwakePeaks);
Sleep = sum(TDsleepPeaks);

theta = [];
here = 9:-1:4;
for i = 1:length(Wake)
    theta = [theta ones(1,Wake(i)).*bins(here(i))];
end
here = [10 11 12 1 2 3];
for i = 1:length(Sleep)
    theta = [theta ones(1,Sleep(i)).*bins(here(i))];
end

figure
polarhistogram(deg2rad(theta),deg2rad(0:30:360),'Normalization','probability')
thetaticks(bins)
thetaticklabels({'\theta','\alpha-\sigma','\beta','\beta','\alpha-\sigma',...
    '\theta','\delta2','\delta1','s','s','\delta1','\delta2'})
xAX = get(gca,'ThetaAxis');
set(xAX,'FontSize', 24)
title('Neurotypical')
print('-dsvg','./Figures/TDpeakHistogram.svg')

%% AS
Wake = sum(ASwakePeaks);
Sleep = sum(ASsleepPeaks);

theta = [];
here = 9:-1:4;
for i = 1:length(Wake)
    theta = [theta ones(1,Wake(i)).*bins(here(i))];
end
here = [10 11 12 1 2 3];
for i = 1:length(Sleep)
    theta = [theta ones(1,Sleep(i)).*bins(here(i))];
end

figure
polarhistogram(deg2rad(theta),deg2rad(0:30:360),'Normalization','probability')
thetaticks(bins)
thetaticklabels({'\theta','\alpha-\sigma','\beta','\beta','\alpha-\sigma',...
    '\theta','\delta2','\delta1','s','s','\delta1','\delta2'})
xAX = get(gca,'ThetaAxis');
set(xAX,'FontSize', 24)
title('Angelman')
print('-dsvg','./Figures/ASpeakHistogram.svg')

%% DS
Wake = sum(DSwakePeaks);
Sleep = sum(DSsleepPeaks);

theta = [];
here = 9:-1:4;
for i = 1:length(Wake)
    theta = [theta ones(1,Wake(i)).*bins(here(i))];
end
here = [10 11 12 1 2 3];
for i = 1:length(Sleep)
    theta = [theta ones(1,Sleep(i)).*bins(here(i))];
end

figure
polarhistogram(deg2rad(theta),deg2rad(0:30:360),'Normalization','probability')
thetaticks(bins)
thetaticklabels({'\theta','\alpha-\sigma','\beta','\beta','\alpha-\sigma',...
    '\theta','\delta2','\delta1','s','s','\delta1','\delta2'})
xAX = get(gca,'ThetaAxis');
set(xAX,'FontSize', 24)
title('Dup15q')
print('-dsvg','./Figures/DSpeakHistogram.svg')

%%
edges = -1:0.5:5;
myfigure
histogram(log2(TDwakePeaksExact),edges)
histogram(log2(TDsleepPeaksExact),edges)
legend({'Wake','NREMS'},'fontsize',18,'location','northwest')
legend boxoff
xticks(-1:5)
xticklabels(2.^(-1:5))
xlabel('Peak frequency (Hz)')
ylabel('Count')
title('NT peak frequency','fontsize',18)
makefighandsome
print('-dsvg','./Figures/TDpfhist.svg')

myfigure
histogram(log2(ASwakePeaksExact),edges)
histogram(log2(ASsleepPeaksExact),edges)
legend({'Wake','NREMS'},'fontsize',18,'location','northwest')
legend boxoff
xticks(-1:5)
xticklabels(2.^(-1:5))
xlabel('Peak frequency (Hz)')
ylabel('Count')
title('AS peak frequency','fontsize',18)
makefighandsome
print('-dsvg','./Figures/ASpfhist.svg')

myfigure
histogram(log2(DSwakePeaksExact),edges)
histogram(log2(DSsleepPeaksExact),edges)
legend({'Wake','NREMS'},'fontsize',18,'location','northwest')
legend boxoff
xticks(-1:5)
xticklabels(2.^(-1:5))
xlabel('Peak frequency (Hz)')
ylabel('Count')
title('Dup15q peak frequency','fontsize',18)
makefighandsome
print('-dsvg','./Figures/DSpfhist.svg')

%% Export data

load AllPeakData

Group = [repmat({'AS'},size(ASwakePeaks,1),1); repmat({'NT'},size(TDwakePeaks,1),1); ...
    repmat({'Dup15q'},size(DSwakePeaks,1),1)];

Condition = [repmat({'Wake'},length(Group),1); repmat({'Sleep'},length(Group),1)];

Group = [Group; Group];

slow = [ASwakePeaks(:,1); TDwakePeaks(:,1); DSwakePeaks(:,1); ...
    ASsleepPeaks(:,1); TDsleepPeaks(:,1); DSsleepPeaks(:,1)];

delta1 = [ASwakePeaks(:,2); TDwakePeaks(:,2); DSwakePeaks(:,2); ...
    ASsleepPeaks(:,2); TDsleepPeaks(:,2); DSsleepPeaks(:,2)];

delta2 = [ASwakePeaks(:,3); TDwakePeaks(:,3); DSwakePeaks(:,3); ...
    ASsleepPeaks(:,3); TDsleepPeaks(:,3); DSsleepPeaks(:,3)];

theta = [ASwakePeaks(:,4); TDwakePeaks(:,4); DSwakePeaks(:,4); ...
    ASsleepPeaks(:,4); TDsleepPeaks(:,4); DSsleepPeaks(:,4)];

alpha_sigma = [ASwakePeaks(:,5); TDwakePeaks(:,5); DSwakePeaks(:,5); ...
    ASsleepPeaks(:,5); TDsleepPeaks(:,5); DSsleepPeaks(:,5)];

beta = [ASwakePeaks(:,6); TDwakePeaks(:,6); DSwakePeaks(:,6); ...
    ASsleepPeaks(:,6); TDsleepPeaks(:,6); DSsleepPeaks(:,6)];

Tout = table(Group,Condition,slow,delta1,delta2,theta,alpha_sigma,beta)

writetable(Tout,'./Fig2ABC_data.csv')
