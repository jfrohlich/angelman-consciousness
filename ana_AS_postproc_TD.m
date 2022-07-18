%function[] = ana_AS_postproc_TD()
%%
clearvars
dbstop if error

td = dir('./TD_controls/Imported/Detrend_butterHP_firLP/Clean/*.mat');
tdpath = './TD_controls/Imported/Detrend_butterHP_firLP/Clean/';
stNDX = 1;

if stNDX ~= 1
    warning('Start index is set to %i\n',stNDX)
end

for ifl = stNDX:length(td)
    postproc(tdpath,td(ifl).name);
end

function[] = postproc(lddrc,fileStr)

% Joel Frohlich
% 05-25-2020
%
% lddrc is the directory to load the file from 
%
% fileStr is string with name of the file
%
% norm is string with the space in which the normalization should be
% performed
%
% smooth = smoothing factor (1/3 recommended, use fraction and not decimal
% approximation)
%
% srInt - freq to start the normalization
%
% spInt - freq to stop the normalization
% 
% shift - scale shift factor

load angelman_lay.mat lay

% Boolean var indicted if the interpolation happened
intp = 0;

DIRRESULT = './2021_analysis/postprocessed/TD/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end

load(strcat(lddrc,fileStr))
fprintf('\n Now loading %s \n',fileStr)
% If structure contains ICA field
if isfield(data.cfg,'ica')
    filt = data;
    % Apply ICA
    filt.label = upper(data.label(data.cfg.ica.userawch));
    idcompbad = setdiff(1:data.cfg.ica.numofic,data.cfg.ica.ctype.(data.cfg.ica.ctypemainlabel).num);
    % ica weights matrix
    weights = data.cfg.ica.a(:,idcompbad)*data.cfg.ica.w(idcompbad,:);

    ICA_cleaned = filt.trial{1}(data.cfg.ica.userawch,:) - ...
        weights*filt.trial{1}(data.cfg.ica.userawch,:);
    filt.trial{1} = ICA_cleaned;
    display('ICA subtracted')
else % if we need to fetch it from the old files
    try
        icafile = strcat(data.cfg.fileout,'_ICA');
    catch
        icafile = strcat(fileStr,'_ICA');
    end
    % rename current data to specify it is the (correctly) filtered one
    filt = data; 
    % load ICA data
    load(strcat('AllProcessedFiles_1_32_Hz/',icafile));
    assert(length(data.trial) == 1,'Wrong number of trials')
    % Apply ICA
    filt.label = upper(data.label(data.cfg.ica.userawch));
    idcompbad = setdiff(1:data.cfg.ica.numofic,data.cfg.ica.ctype.(data.cfg.ica.ctypemainlabel).num);
    % ica weights matrix
    weights = data.cfg.ica.a(:,idcompbad)*data.cfg.ica.w(idcompbad,:);
    
    % Apply to first trial (awake)
    ICA_cleaned = filt.trial{1}(data.cfg.ica.userawch,:) - ...
        weights*filt.trial{1}(data.cfg.ica.userawch,:);
    filt.trial{1} = ICA_cleaned;
    
    % Apply to second trial (sleep)
    ICA_cleaned = filt.trial{2}(data.cfg.ica.userawch,:) - ...
        weights*filt.trial{2}(data.cfg.ica.userawch,:);
    filt.trial{2} = ICA_cleaned;
    
    display('ICA subtracted')
end

load Angelman_lay
elec = lay.cfg.elec;

filt.cfg.raw = data.cfg.raw;

% fill in the missing channel as nan
filt2 = filt;
filt2 = rmfield(filt2,'trial');
filt2 = rmfield(filt2,'label');
for ielec=1:length(elec.label)
    idx = find(strcmpi(elec.label{ielec},filt.label));
    if ~isempty(idx)
        filt2.trial{1}(ielec,:) = filt.trial{1}(idx,:);
        filt2.trial{2}(ielec,:) = filt.trial{2}(idx,:);
        filt2.label{ielec} = elec.label{ielec};
    else
        filt2.trial{1}(ielec,:) = nan(1,size(filt.trial{1},2));
        filt2.trial{2}(ielec,:) = nan(1,size(filt.trial{2},2));
        filt2.label{ielec} = elec.label{ielec};
    end
end


try
    badch = setdiff(upper(lay.label),upper(data.cfg.raw.mainplot.allchlabel));
catch %if the line above doesn't work
    badch = data.cfg.raw.ctype.bad.label;
end

if ~isempty(badch)
    clear cfg
    cfg.badchannel=badch;
    cfg.missingchannel = [];
    % First trial (awake)
    out = jf_interpchan(filt2.trial{1},cfg);
    filt2.trial{1} = out;
    % Second trial (sleep)
    out = jf_interpchan(filt2.trial{2},cfg);
    filt2.trial{2} = out;
end

data = filt2;

%% redo average ref 
data.trial{1} = data.trial{1} - mean(data.trial{1},1);
data.trial{2} = data.trial{2} - mean(data.trial{2},1);


%% Save postprocessed data 

% find file suffix
NDX = strfind(fileStr,'_ICA');
if ~isempty(NDX)
    outstr = fileStr(1:NDX-1);
else
    outstr = fileStr(1:end-4);
end
data.fstr = outstr;
outname = strcat(DIRRESULT,outstr,'_pstprc_',date,'.mat');
save(outname,'data','-v7.3');

clear data

end


   