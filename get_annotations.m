% compares our sleep annotations with those done by a trained neurologist
% through Ben Philpot's lab

% NOTE: some subjects that Philphot lists as only having one visit might
% have multiple visits if one or more occur after late 2017 or early 2018.
% Den Bakker et al. 2018 was received 22 November 2017 and accepted 11 
% April 2018, so it should be safe to assume there was no date from 2018 on
% when this was written. 

clearvars
anno1 = dir('./Philpot_sleep_annotations/*.xls');
anno2 = dir('./Nespeca_sleep_annotations/*.xls');
anno = [anno1(:).' anno2(:).']; % combine Philpot and Nespeca annotations

EEG = dir('./postprocessed/labeled/*.mat');
overlap = [];
dbstop if error

% Remove duplicates

rmidx = [];
for i = 1:length(anno)
    if contains(anno(i).name,'._')
        rmidx = [rmidx i];
    end
end
    
anno(rmidx) = []; % remove extraneous elements

Nespeca_IDs = [101132 102786 106234 106635 100382 110429 110882 120206 124283 125120 128280 200007];

skip = [100580 105499]; % unclear which visits go with which annotations 
SIDs = []; % store SIDs 
Pvals = []; % from Wilcoxin rank-sum test
lost = cell(1,1); % subjects not found;
p_spln = cell(1,1); % store data length of sleep, use SIDs are key

PhasbutnotM = 0;
MhasbutnotP = 0;
Ponly = cell(1,1);
Monly = cell(1,1);

lcnt = 0; % longitudinal data
scnt = 0; % single visit data
xcnt = 0; % file not found
fcnt = 0; % files
ncnt = 0; % files for table (to guide postprocessing)

Tsw = table(); % create new table
istheresleep = []; 

all_slp_subs = cell(1,1);
has_philpot_sleep = cell(1,1);
has_philpot_sleepwake = cell(1,1);
hscnt = 0;
hswcnt = 0;


for ifile = 1:length(anno)
    ifile
    %%% DEBUG %%%
    if ifile == 81
        keyboard
    end
    %%%%%%
    if strcmp(anno(ifile).name(1:2),'._')
        anno(ifile).name = anno(ifile).name(3:end);
    end
    try
        T = readtable(sprintf('./Philpot_sleep_annotations/%s',anno(ifile).name));
    catch
        % If it's not from Philpot, it must be from Nespeca
        T = readtable(sprintf('./Nespeca_sleep_annotations/%s',anno(ifile).name));
    end
    
    if iscell(T.Type) % If there's anything to read in the table (check this first other we get error when table is blank)
        if any(contains(T.Type,'S')) && any(contains(T.Type,'W')) % if wake and sleep contained in same file
            hswcnt = hswcnt + 1;
            has_philpot_sleepwake{hswcnt} = anno(ifile).name;
        end
        
        if any(contains(T.Type,'S'))
            hscnt = hscnt + 1;
            has_philpot_sleep{hscnt} = anno(ifile).name;
        end
    end
    
    if length(anno(ifile).name) == 20 % these are the subjects w/out longitudinal data
        if strcmp(anno(ifile).name(1:2),'._')
            anno(ifile).name = anno(ifile).name(3:end);
        end
                 
        % get SID
        start = find(ismember(anno(ifile).name,'_')) + 1;
        stop = find(ismember(anno(ifile).name,'.')) - 1;
        SID1 = str2double(anno(ifile).name(start:stop));
                   
        if any(ismember(skip,SID1))
            fprintf('\nSkipping this one, unclear which visit goes with which annotations\n')
            clear subject SID1 SID2 date1 date2
            continue
        end
        % find file with matching SID
        for ieeg = 1:length(EEG)
            NDX = find(ismember(EEG(ieeg).name,'_'),2);
            start = NDX(1)+1; stop = NDX(2)-1;
            SID2 = str2double(EEG(ieeg).name(start:stop));
            if SID1 == SID2
                subject = EEG(ieeg).name;
                scnt = scnt + 1;                                    
                % if this is a file with sleep
                if strcmp(has_philpot_sleep{hscnt},anno(ifile).name) 
                    fcnt = fcnt + 1;
                    all_slp_subs{fcnt} = subject;
                end
                
                % Add rows to table
                if ~contains(all_slp_subs{fcnt},'SEDATED') % if this subject wasn't sedated
                    ncnt = ncnt + 1;
                    usidx = find(all_slp_subs{fcnt} == '_',3);
                    Tsw.file{ncnt} = all_slp_subs{fcnt}(1:usidx(3)-1);
                    %%% DEBUG %%%
                    if ncnt > 1 && strcmp(Tsw.file{ncnt},Tsw.file{ncnt-1})
                        %keyboard
                        continue % skip repeat subject
                    end
                    %%%%%%
                    SLP = [];
                    WK = [];
                    DRWY = [];
                    for irow = 1:size(T,1) % look through annotation table
                        if strcmp(T.Type(irow),'S')
                            SLP = [SLP; T.Time_Start_sec_(irow) ...
                                T.Time_End_sec_(irow)];
                        elseif strcmp(T.Type(irow),'W')
                            WK = [WK; T.Time_Start_sec_(irow) ...
                                T.Time_End_sec_(irow)];
                        elseif strcmp(T.Type(irow),'D')
                            DRWY = [DRWY; T.Time_Start_sec_(irow) ...
                                T.Time_End_sec_(irow)];
                        end
                    end
                                       
                    istheresleep(ncnt) = ~isempty(SLP);
                    Tsw.Sleep{ncnt} = SLP;
                    Tsw.Wake{ncnt} = WK;
                    Tsw.Drowsy{ncnt} = DRWY;
                end
                
                if isempty(SIDs) || ~any(ismember(SIDs,SID1))
                    SIDs = [SIDs SID1];
                else
                    fprintf('Skipping this repeat of a subject without longitudinal data')
                    continue
                end
                break
            end
        end
    elseif length(anno(ifile).name) >= 29 % subjects with longitudinal data
        if strcmp(anno(ifile).name(1:2),'._')
            anno(ifile).name = anno(ifile).name(3:end);
        end
        
        % get SID
        NDX = find(ismember(anno(ifile).name,'_'),2);
        start = NDX(1) + 1;
        stop = start + 5;
        SID1 = str2double(anno(ifile).name(start:stop));
        if any(ismember(skip,SID1))
            fprintf('\nSkipping this one, unclear which visit goes with which annotations\n')
            continue
        end
        clear start stop
        % get date
        y = str2double(anno(ifile).name(NDX(2)+1:NDX(2)+4));
        m = str2double(anno(ifile).name(NDX(2)+5:NDX(2)+6));
        d = str2double(anno(ifile).name(NDX(2)+7:NDX(2)+8));
        date1 = datetime(y,m,d);
        clear NDX y m d
        % find file with matching SID
        for ieeg = 1:length(EEG)
            NDX = find(ismember(EEG(ieeg).name,'_'),2);
            start = NDX(1)+1; 
            stop = NDX(2)-1;
            SID2 = str2double(EEG(ieeg).name(start:stop));
            % get date
            y = str2double(EEG(ieeg).name(NDX(2)+1:NDX(2)+4));
            m = str2double(EEG(ieeg).name(NDX(2)+5:NDX(2)+6));
            d = str2double(EEG(ieeg).name(NDX(2)+7:NDX(2)+8));
            
            date2 = datetime(y,m,d);
            if SID1 == SID2 && eq(date1,date2)
                subject = EEG(ieeg).name;
                
                % if this is a file with sleep
                if strcmp(has_philpot_sleep{hscnt},anno(ifile).name) 
                    fcnt = fcnt + 1;
                    all_slp_subs{fcnt} = subject;
                end
                
                % Add rows to table
                if ~contains(all_slp_subs{fcnt},'SEDATED')% if this subject wasn't sedated
                    ncnt = ncnt + 1;
                    usidx = find(all_slp_subs{fcnt} == '_',3);
                    Tsw.file{ncnt} = all_slp_subs{fcnt}(1:usidx(3)-1);
                    %%% DEBUG %%%
                    if ncnt > 1 && strcmp(Tsw.file{ncnt},Tsw.file{ncnt-1})
                        %keyboard
                        continue % skip repeat subject
                    end
                    %%%%%%
                
%                 % if this is a file with both sleep 
%                 if ~contains(all_slp_subs{fcnt},'SEDATED') %&& ... % if this subject wasn't sedated
%                         %strcmp(has_philpot_sleep{hscnt},anno(ifile).name)% and it has both sleep and wake
%                     fcnt = fcnt + 1;
%                     all_slp_subs{fcnt} = subject
%                     
%                     ncnt = ncnt + 1;
%                     usidx = find(all_slp_subs{fcnt} == '_',3);
%                     Tsw.file{ncnt} = all_slp_subs{fcnt}(1:usidx(3)-1);
%                     %%% DEBUG %%%
%                     if ncnt > 1 && strcmp(Tsw.file{ncnt},Tsw.file{ncnt-1})
%                         keyboard
%                         continue % skip repeat subject
%                     end
%                     %%%%%%
                    SLP = [];
                    WK = [];
                    DRWY = [];
                    for irow = 1:size(T,1) % look through annotation table
                        if strcmp(T.Type(irow),'S')
                            SLP = [SLP; T.Time_Start_sec_(irow) ...
                                T.Time_End_sec_(irow)];
                        elseif strcmp(T.Type(irow),'W')
                            WK = [WK; T.Time_Start_sec_(irow) ...
                                T.Time_End_sec_(irow)];
                        elseif strcmp(T.Type(irow),'D')
                            DRWY = [DRWY; T.Time_Start_sec_(irow) ...
                                T.Time_End_sec_(irow)];
                        end
                    end
                    
                    istheresleep(ncnt) = ~isempty(SLP);
                    Tsw.Sleep{ncnt} = SLP;
                    Tsw.Wake{ncnt} = WK;
                    Tsw.Drowsy{ncnt} = DRWY;
                end
                
                lcnt = lcnt + 1;
                if ~ismember(SIDs,SID1)
                    SIDs = [SIDs SID1];
                end
                break
            end
        end
    else
        error('Invalid file name')
    end
    if exist('subject','var')
        load(sprintf('%s%s','./postprocessed/labeled/',subject),'data');
        SLEEP_PHILPOT = [];
        philpot = [];
        % check if table is empty
        if ~iscell(T.Time_Start_sec_) && sum(isnan(T.Time_Start_sec_)) == size(T,1)
            fprintf('\nAnnotation table from Philpot is blank, skipping this one ...\n')
            clear subject SID1 SID2 date1 date2
            continue
        end
        for irow = 1:size(T,1)
            if strcmp(T.Type(irow),'S')
                SLEEP_PHILPOT = [SLEEP_PHILPOT; T.Time_Start_sec_(irow) ...
                    T.Time_End_sec_(irow)];
                philpot = [philpot  T.Time_Start_sec_(irow):T.Time_End_sec_(irow)];
            end
        end
        SLEEP_MONTI = round(data.cfg.dattype.SLEEP_MONTI./data.fsample);
        monti = [];
        % Extract samples from Monti data
        for irow = 1:size(SLEEP_MONTI,1)
            monti = [monti SLEEP_MONTI(irow,1):SLEEP_MONTI(irow,2)];
        end
        if isempty(monti) && isempty(philpot)
            fprintf('Skipping this one, they''re both empty\n')
            clear subject SID1 SID2 date1 date2
            continue
        end
        
        % Compute percent overlap in sleep sections
        shared = length(intersect(monti,philpot));
        PM = shared/length(monti)*100; % percent overlap ref'ed to Monti
        PP = shared/length(philpot)*100; % percent overlap ref'ed to Philpot
        overlap = [overlap; PM PP];
        
        if ~isempty(monti) && ~isempty(philpot)
            % Wilcoxin rank sum test
            P = ranksum(monti,philpot,'alpha',0.05);
            Pvals = [Pvals P];
        end

        if isempty(monti) && ~isempty(philpot)
            PhasbutnotM = PhasbutnotM + 1;
            if exist('date1','var')
                fstr = sprintf('%i %s',SID1,datestr(date1));
            else
                fstr = num2str(SID1);
            end
            Ponly{PhasbutnotM} = fstr;
        elseif ~isempty(monti) && isempty(philpot)
            MhasbutnotP = MhasbutnotP + 1;
            if exist('date1','var')
                fstr = sprintf('%i %s',SID1,datestr(date1));
            else
                fstr = num2str(SID1);
            end
            Monly{MhasbutnotP} = fstr;
        end
%         % store length of sleep from Philpot annotations
%         I = find(ismember(SIDs,SID1));
%         if I > length(p_spln) | isempty(I) | isempty(p_spln{I}) % if it's a new subject
%             p_spln{I} = [length(philpot)];
%         else                                        % if it's a subject already in the cell array
%             p_spln{I}(length(p_spln{I})+1) = length(philpot);
%         end
%         % fprintf('\nClearing old vars...\n')
        clear subject SID1 SID2 date1 date2
    else
        xcnt = xcnt + 1;
        if exist('date1','var') && ~isnat(date1)
            fprintf('\nNo EEG for subject with Philpot annotation: %i %s\n',SID1,datestr(date1))
            lost{xcnt} = strcat(num2str(SID1),'_',datestr(date1));
        else
            fprintf('\nNo EEG for subject with Philpot annotation: %i \n',SID1)
            lost{xcnt} = num2str(SID1);
        end
        continue
    end
end

% %% Report results
% 
% % find subjects with at least 30 s sleep in Philpot annotations
% gcnt = 0;
% hassleep = 0;
% for i = 1:length(p_spln)
%     m = min(p_spln{i});
%     if m >= 30
%         gcnt = gcnt + 1;
%     end
%     if ~isempty(p_spln{i}) && m > 0
%         hassleep = hassleep + 1;
%     end
% end
% 
% sleepboth = sum(nanmean(overlap,2) > 0); % sleep in both annotations
% 
% overlap2 = overlap;
% overlap2(overlap2==0) = NaN; % replace 0s with NaNs ... this excludes pairing where one lab had sleep and the othere did not
% 
% 
% figure
% histogram(overlap2(:,1),10,'binwidth',10) % percent overlap ref'ed to Monti
% title('Monti sleep shared with Philpot','fontsize',18)
% xlabel('Percent shared')
% ylabel('Files')
% makefigpretty
% 
% figure
% histogram(overlap2(:,2),10,'binwidth',10) % percent overlap ref'ed to Philpot
% title('Philpot sleep shared with Monti','fontsize',18)
% xlabel('Percent shared')
% ylabel('Files')
% makefigpretty
% 
% assert(length(SIDs) == length(unique(SIDs)),'Repeat subjects in SID list')
% 
% Q = mafdr(Pvals,'BHFDR',true);
% H = sum(Q < 0.05);
% 
% fprintf('\n%i files total\n',length(SIDs))
% fprintf('\n%i shared subjects with sleep in Philpot annotations\n',hassleep)
% fprintf('\n%i shared subjects that have sleep according to both annotations\n',sleepboth)
% fprintf('\n%i shared subjects with at least 30 s sleep (according to Philpot annotations)\n',gcnt)
% fprintf('\nMean percentage of Monti sleep shared with Philpot: %2.2f +/- %2.2f percent\n',nanmean(overlap2(:,1)),nanstd(overlap2(:,1)))
% fprintf('\nMean percentage of Philpot sleep shared with Monti: %2.2f +/- %2.2f percent\n',nanmean(overlap2(:,2)),nanstd(overlap2(:,2)))
% fprintf('\nMedian percentage of Monti sleep shared with Philpot: %2.2f +/- %2.2f percentpercent\n',nanmedian(overlap2(:,1)),nanstd(overlap2(:,1)))
% fprintf('\nMedian percentage of Philpot sleep shared with Monti: %2.2f +/- %2.2f percentpercent\n',nanmedian(overlap2(:,2)),nanstd(overlap2(:,2)))
% fprintf('\n%i out of %i sleep labeling significantly differ between labs (Wilcoxin rank-sum, FDR corrected)\n',H,length(Q))
% fprintf('\n%i EEGs have sleep annotations from Philpot but not from Monti\n',PhasbutnotM)
% fprintf('\n%i EEGs have sleep annotations from Monti but not from Philpot\n',MhasbutnotP)
% 
% %% How many files have sleep?
% 
% hsID = [];
% 
% for i = 1:length(has_philpot_sleep)
%     sleepfile = has_philpot_sleep{i};
%     idx = find(has_philpot_sleep{1}=='_',1);
%     hsID = [hsID str2num(sleepfile(idx+1:idx+6))];
% end
% 
% fprintf('\n%i files from %i subjects with sleep from Philpot\n',...
%     length(has_philpot_sleep),length(unique(hsID)))
% 
% fprintf('\n%i total subjects with sleep from Philpot and Nespeca combinded\n',...
%     length(unique(union(hsID,Nespeca_IDs))))
% 
% unique(union(hsID,Nespeca_IDs))'

%% Save the names of files with paired sleep-wake

fprintf('%i files with sleep\n',sum(istheresleep))

writetable(Tsw,'./Nespeca_sleep_annotations/Philpot_Nespeca_sleep.csv')
 
save PhilpotNespecaSleep Tsw % save the table
        
        



