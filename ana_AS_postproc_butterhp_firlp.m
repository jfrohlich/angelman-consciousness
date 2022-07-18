%function[] = ana_AS_postproc_butterhp_firlp()
clearvars
dbstop if error

load angelman_lay.mat lay

AS = dir('./2021_analysis/butterHPfirLP_2021_wICA/NewICA/*.mat');
aspath = './2021_analysis/butterHPfirLP_2021_wICA/NewICA/';

stNDX = 1;

if stNDX ~= 1
    warning('Start index is set to %i\n',stNDX)
end

for ifl = stNDX:length(AS)
    postproc(aspath,AS(ifl).name);
end

function[] = postproc(aspath,fileStr)

load angelman_lay.mat lay

% Boolean var indicted if the interpolation happened
% intp = 0;

DIRRESULT = './2021_analysis/postprocessed/AS/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end

% we need this to get the sleep field, which might be missing from some
% fixed files
oldpath = './2021_analysis/butterHPfirLP_2021_wICA/';

assert(contains(fileStr,'fixed'),'This is not the file with ICA redone')
load(strcat(aspath,fileStr))
fprintf('\n Now loading %s \n',fileStr)

old = load(strcat(oldpath,sprintf('%s23-Jun-2021',fileStr(1:end-21))));

% combine my sleep markings (SLEEP_MONTI) with Mark's (SLEEP_NESPECA)
if isfield(data.cfg.dattype,'SLEEP_NESPECA')
    data.cfg.dattype.SLEEP_ALL = [data.cfg.dattype.SLEEP_NESPECA; data.cfg.dattype.SLEEP_MONTI];
else
    data.cfg.dattype.SLEEP_ALL = data.cfg.dattype.SLEEP_MONTI;
end

% add 30 s buffers around MONTI_ALL, to make sure that sleep or drowsy
% doesn't get counted as awake

buffer = [];

for irow = 1:size(data.cfg.dattype.SLEEP_ALL,1)
    x = data.cfg.dattype.SLEEP_ALL(irow,1);
    y = data.cfg.dattype.SLEEP_ALL(irow,2);
    buffer = [buffer; x-30*data.fsample x-1; y+1 y+30*data.fsample];
end

% add buffer to the "drowsy" sections
data.cfg.dattype.drowsy = [data.cfg.dattype.drowsy; buffer];

% check that there is nothing in the sleep markings from Roche (which
% really indicate "Not sleep") that we haven't already covered here with
% either the SLEEP_ALL or drowsy fields

D = []; % drowsy
S = []; % sleep
NW = []; % not wake
WN = []; % wake sections identified by nespeca
SN = []; % sleep sections identified by nespeca

for irow = 1:size(data.cfg.dattype.SLEEP_ALL,1)
    x = data.cfg.dattype.SLEEP_ALL(irow,1);
    y = data.cfg.dattype.SLEEP_ALL(irow,2);
    S = [S x:y];
end
clear x y

for irow = 1:size(data.cfg.dattype.drowsy,1)
    x = data.cfg.dattype.drowsy(irow,1);
    y = data.cfg.dattype.drowsy(irow,2);
    D = [D x:y];
end
clear x y

if isfield(data.cfg.dattype,'AWAKE_NESPECA')
    for irow = 1:size(data.cfg.dattype.AWAKE_NESPECA,1)
        x = data.cfg.dattype.AWAKE_NESPECA(irow,1);
        y = data.cfg.dattype.AWAKE_NESPECA(irow,2);
        WN = [WN x:y];
    end
end
clear x y

if isfield(data.cfg.dattype,'SLEEP_NESPECA')
    for irow = 1:size(data.cfg.dattype.SLEEP_NESPECA,1)
        x = data.cfg.dattype.SLEEP_NESPECA(irow,1);
        y = data.cfg.dattype.SLEEP_NESPECA(irow,2);
        SN = [SN x:y];
    end
end
clear x y

for irow = 1:size(old.data.cfg.dattype.sleep,1)
    x = old.data.cfg.dattype.sleep(irow,1);
    y = old.data.cfg.dattype.sleep(irow,2);
    NW = [NW x:y];
end
clear x y

% Sanity check 
assert(isempty(intersect(SN,WN)),'Nespeca labels disagree with themselves')

%%%%%
% check that Mark's awake labels don't overlap with anything we marked as
% sleep or drowsy
notsleep = intersect(WN,S);
notdrowsy = intersect(WN,D);

% It appears that the only file(s) with serious disagreement with Mark
% Nespeca's awake labels have a "sleep" section completely encompassed by
% wake--thus, the patch below should fix. For more complicated scenarios,
% an error will be generated. 
if length(notsleep) > data.fsample/2 % If labels disagree by more than a trivial amount (0.5 sec)
    if min(WN) < min(S) & max(WN) > max(S) % if the section Nespeca marked as awake starts before sleep and ends after sleep
        S = [];
        data.cfg.dattype.SLEEP_ALL = [];
        data.cfg.dattype.SLEEP_MONTI = [];
    else
        error('Labels disagree')
    end
end

if length(notdrowsy) > data.fsample/2 % If labels disagree by more than a trivial amount (0.5 sec)
    if min(WN) < min(D) & max(WN) > max(D) % if the section Nespeca marked as awake starts before drowsy and ends after drowsy
        D = [];
        data.cfg.dattype.drowsy = [];
    else
        error('Labels disagree')
    end
end
%%%%%

% check for anything we didn't cover with our new labels
uncovered = setdiff(NW,union(S,D));

if ~isempty(uncovered) 
    breakpoints = find(diff(uncovered)~=1);
    drowsy2 = [];
    
    for i = 1:length(breakpoints)+1
        if i == 1
            A = 1;
        else
            A = breakpoints(i-1)+1;
        end
        if i == length(breakpoints)+1
            B = length(uncovered);
        else
            B = breakpoints(i);
        end
        drowsy2 = [drowsy2; uncovered(A) uncovered(B)];
    end
    data.cfg.dattype.drowsy = [data.cfg.dattype.drowsy; drowsy2];
end

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
    disp('ICA subtracted')
else % if we need to fetch it from the old files
    error('No ICA!')
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
        filt2.label{ielec} = elec.label{ielec};
    else
        filt2.trial{1}(ielec,:) = nan(1,size(filt.trial{1},2));
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
    out = jf_interpchan(filt2.trial{1},cfg);
    filt2.trial{1} = out;
end

data = filt2;

%% redo average ref 

data.trial{1} = data.trial{1} - mean(data.trial{1},1);


%% Save postprocessed data 

% find file suffix
% NDX = strfind(fileStr,'_ICA');
% if ~isempty(NDX)
%     outstr = fileStr(1:NDX-1);
% else
%     outstr = fileStr(1:end-4);
% end

outstr = fileStr(1:21);
data.fstr = outstr;
outname = strcat(DIRRESULT,outstr,'pstprc_',date,'.mat');
save(outname,'data','-v7.3');

clear data

end


   