clearvars
rng(45783) % seed added 02.26.19
OS = computer; % detect opperating system
dbstop if error
ft_defaults

switch OS
    case 'GLNXA64'
        pth = './2021_analysis/postprocessed/TD/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
    case 'MACI64'
        pth = './2021_analysis/postprocessed/TD/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
    case 'PCWIN64'
        pth = './2021_analysis/postprocessed/TD/';
        files = dir(sprintf('%s%s',pth,'*.mat'));
end

sleep = true; % asleep or awake?
compute_wPLI(sleep,files,pth)

sleep = false; % asleep or awake?
compute_wPLI(sleep,files,pth)


%%

function[] = compute_wPLI(sleep,files,pth)

switch sleep
    case true
        DIRRESULT = './2021_analysis/TD_output/wPLI/sleep/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end   
        ntrl = 2;
    case false       
        DIRRESULT = './2021_analysis/TD_output/wPLI/wake/'; if ~exist(DIRRESULT), mkdir(DIRRESULT), end  
        ntrl = 1; 
end

tic
start_ndx = 1;

if start_ndx ~= 1, warning('Starting index is set to %i',start_ndx), end

for ifile = start_ndx:length(files)
    %if logical(sum(ismember(select_files,files(ifile).name(1:end-4))))
    load(sprintf('%s%s',pth,files(ifile).name),'data')
    fprintf('Now loading %s \n',files(ifile).name)
    
    % check if this file had useable sleep
    
    % The second trial is always sleep, so this might not be necessary
    if ~isfield(data.cfg.dattype,'SLEEP_MONTI') || isempty(data.cfg.dattype.SLEEP_MONTI)
        fprintf('Skipping this one, no useable sleep ... \n')
        continue
    end
    
    %Get datamat
    datamat = data.trial{ntrl};
       
    % Do some santiy checks
    assert(size(datamat,1) == 19,'channel(s) missing!')
    % allow for small number of NaNs at beginning and end of file
    assert(sum(isnan(datamat(1:end)))/numel(datamat) < 0.01, 'Many NaNs in data!')
    
    % Set bad data to NaNs for purpose of taking freq transform
    badart = [data.cfg.dattype.bad; data.cfg.dattype.flash ];
    
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
    
    percent = c/length(datamat);
    durgood = c/data.fsample;
    

    
    %% Compute wPLI (as alternative to wSMI)
    
    n_win =5*data.fsample; % window size (make sure this is at least twice as slow as the slowest oscillation
    overlap = true; % allow windows to overlap by 50%  
    
    % pad with nans to create to a matrix that is evenly divisible for
    % reshaping into 3D array
    datamat = [datamat nan(size(datamat,1),n_win-mod(size(datamat,2),n_win))];
    datamat_3d = reshape(datamat,size(datamat,1),n_win,size(datamat,2)/n_win);
    
    switch overlap % overlapping windows?
        case true   
            % create second datamat to get 50% overlapping windows
            datamat2 = datamat(:,(n_win/2)+1:end-(n_win/2));
            datamat2_3d = reshape(datamat2,size(datamat2,1),n_win,size(datamat2,2)/n_win);

            % check that the second matrix has appropriate dimensions
            assert(size(datamat_3d,1)==size(datamat2_3d,1))
            assert(size(datamat_3d,2)==size(datamat2_3d,2))
            assert(size(datamat_3d,3)==size(datamat2_3d,3)+1)

            % interleave the two 3D arrays so we acheive 50% window overlap

            data_3d = nan(size(datamat_3d,1),size(datamat_3d,2),size(datamat_3d,3)*2-1);
            for ilay = 1:size(data_3d,3)
                if mod(ilay,2) == 1 % if it's an odd numbered layer
                    data_3d(:,:,ilay) = datamat_3d(:,:,ceil(ilay/2));
                else % if it's an even numbered layer
                    data_3d(:,:,ilay) = datamat2_3d(:,:,ilay/2);
                    % check that they really overlap
                        % if this and the last layer don't have nans
                    if ~any(isnan(data_3d(:,:,ilay)),[1 2]) && ~any(isnan(data_3d(:,:,ilay-1)),[1 2])
                        assert(all(data_3d(:,1:n_win/2,ilay) == data_3d(:,n_win/2+1:n_win,ilay-1),[1 2]))
                    end
                end
            end
        case false
            data_3d = datamat_3d;
    end
    
    [wpli,dwpli,foi] = ro_wPLI(data_3d,data.fsample);
   
    wPLIout.prc_used_raw = percent; % percent of data used
    wPLIout.dur_used_raw = durgood; % seconds of data used
    wPLIout.win_len      = n_win/data.fsample; % window size in seconds
    wPLIout.overlap      = overlap; % overlapping windows? 
    wPLIout.fsample      = data.fsample;
    wPLIout.wPLI         = wpli; % wPLI all frequencies
    wPLIout.dwPLI        = dwpli; % dwPLI all frequencies
    
    % wPLI by frequency
    wPLIout.wPLI_slow   = mean(wpli(:,:,foi < 1),3);
    wPLIout.wPLI_delta1 = mean(wpli(:,:,foi >= 1 & foi < 2),3);
    wPLIout.wPLI_delta2 = mean(wpli(:,:,foi >= 2 & foi < 4),3);
    wPLIout.wPLI_theta = mean(wpli(:,:,foi >= 4 & foi < 8),3);
    wPLIout.wPLI_alpha = mean(wpli(:,:,foi >= 8 & foi < 16),3);
    wPLIout.wPLI_beta = mean(wpli(:,:,foi >= 16 & foi < 32),3);    
    
    % dwPLI by frequency
    wPLIout.dwPLI_slow   = mean(dwpli(:,:,foi < 1),3);
    wPLIout.dwPLI_delta1 = mean(dwpli(:,:,foi >= 1 & foi < 2),3);
    wPLIout.dwPLI_delta2 = mean(dwpli(:,:,foi >= 2 & foi < 4),3);
    wPLIout.dwPLI_theta = mean(dwpli(:,:,foi >= 4 & foi < 8),3);
    wPLIout.dwPLI_alpha = mean(dwpli(:,:,foi >= 8 & foi < 16),3);
    wPLIout.dwPLI_beta = mean(dwpli(:,:,foi >= 16 & foi < 32),3);

    
    
    outstr = data.fstr; % cutoff '.mat' suffix
    
    if sleep
        outname = sprintf('%s%s%s%s',DIRRESULT,outstr,'sleep_wPLI_fieldtrip','.mat');
    else 
        outname = sprintf('%s%s%s%s',DIRRESULT,outstr,'wake_wPLI_fieldtrip','.mat');
    end
    
    fprintf('     Now saving...\n')
    save(outname,'wPLIout','-v7.3');
    fprintf('\n')
    clear data datamat wPLIout
end

end