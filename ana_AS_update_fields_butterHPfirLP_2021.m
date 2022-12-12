function[] = ana_AS_update_fields_butterHPfirLP_2021()

%%
dbstop if error

load angelman_lay.mat lay

DIRRESULT = './Monti/2021_analysis/butterHPfirLP_2021_wICA/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end

ICAFileList = ReadMatFileNames('.\Monti\postprocessed\labeled\');

% Get preprocessed data with Butterworth filter
butterFileList   = ReadMatFileNames('.\Monti\AS_butterHP_firLP\'); % get list of files filtered with new filtering

anno1 = dir('./Monti/Philpot_sleep_annotations/*.xls');
anno2 = dir('./Monti/Nespeca_sleep_annotations/*.xls');
anno = [anno1(:).' anno2(:).']; % combine Philpot and Nespeca annotations

% Remove duplicates
rmidx = [];
for i = 1:length(anno)
    if contains(anno(i).name,'._')
        rmidx = [rmidx i];
    end
end
anno(rmidx) = []; % remove extraneous elements

start_ndx = 22;
if start_ndx ~= 1, warning('start_ndx is set equal to %i',start_ndx), end

for ifile = start_ndx:length(butterFileList) % for each file with new filtering
    
    fprintf('\n Updating %s (file %i out of %i) ... \n',butterFileList{ifile},...
        ifile,length(butterFileList))
    filt = load(butterFileList{ifile},'data');
    filt = filt.data;
    A = strfind(butterFileList{ifile},'AS_') + 3; A = A(end);
    B = strfind(butterFileList{ifile},'.mat') - 1;
    filename = butterFileList{ifile}(A:B);
    try
        altname = sprintf('%s_%im',filename(1:6),filt.cfg.age); % in case ICA file uses age and not date
    catch
        altname = sprintf('%s_xyzxyzm',filename(1:6)); % if we can't find the age (placeholder)
    end
    
    T = []; % initialize annotation table
    for iann = 1:length(anno)
        if contains(anno(iann).name,filename)  
            try
                Ttmp = readtable(sprintf('./Monti/Philpot_sleep_annotations/%s',anno(iann).name));
            catch
                Ttmp = readtable(sprintf('./Monti/Nespeca_sleep_annotations/%s',anno(iann).name));
            end
            T = [T; Ttmp]; % combine, e.g., if there are annotations from both Philpot and Nespeca
        end
    end
    
    % replace 0 indicies with 1s
    if ~isempty(T), T.Time_Start_sec_(T.Time_Start_sec_==0) = 1; end
    
    % find ICA file
    idx1 = find(contains(ICAFileList,filename) | contains(ICAFileList,altname));
    if isempty(idx1), idx1 = nan; end
    try
        data2 = load(ICAFileList{idx1},'data');
    catch
        fprintf('\nCan''t load ICA file, skipping this one ...\n')
        continue
    end
    
    data2 = data2.data;
    % copy ICA (no subtraction yet! Still need to run again seeded
    % without old output).
    filt.cfg.ica = data2.cfg.ica;
    data = filt;
    if ~isfield(data,'label')
        assert(length(filt.label) == size(data.trial{1},1),'Channels missing from labels')
        data.label = filt.label;
    end
    
    % Copy over bad channels
    data.cfg.raw.ctype = data2.cfg.raw.ctype;
    
    
    % Get old sleep labels
    if isfield(data2.cfg.dattype,'SLEEP_MONTI')
        data.cfg.dattype.SLEEP_MONTI = data2.cfg.dattype.SLEEP_MONTI;
    else
        data.cfg.dattype.SLEEP_MONTI = [];
    end
    
    % also copy over old sleep field, in includes things like drowsy
    if isfield(data2.cfg.dattype,'sleep')
        data.cfg.dattype.sleep = data2.cfg.dattype.sleep;
    else
        data.cfg.dattype.sleep = [];
    end
    
    if isempty(data.cfg.dattype.SLEEP_MONTI)
    %if isempty(data.cfg.dattype.SLEEP_ALL)
        continue % don't save this file if there's no sleep
    end
    
    % also update bad, muscle, and light flash fields
    data.cfg.dattype.bad = data2.cfg.dattype.bad;
    data.cfg.dattype.muscle = data2.cfg.dattype.muscle;
    data.cfg.dattype.flash = data2.cfg.dattype.flash;
    
    % update channels, in case some are marked bad
    data.cfg.raw = data2.cfg.raw;
    
    % sample info
    data.sampleinfo = [1 length(data.trial{1})];
    data.cfg.trl.trl = [1 length(data.trial{1}) 0];
    
    % make sure all fields are single precision
    data.trial{1} = single(data.trial{1});
    %data.dat_hp = single(data.dat_hp);
    data.time{1} = single(data.time{1});
    
    %% Save new data
    
    data.fstr = filename;
    outname = strcat(DIRRESULT,filename,'_wICA_',date,'.mat');
    fprintf('\n     Saving %s \n',outname)
    save(outname,'data','-v7.3');
    
    clear data data2 filt idx1 idx2
    
end


end
