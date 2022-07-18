function [cfgout] = ro_LZCv(dat,cfg)

% Only compute the vanilla Lempel-Ziv, all other fields will be left empty
%
% Input:
% - dat ... [channels x samples], works with NaN sections
% - cfg ... struct, see script itself for default parameters
% Output:
% - cfgout ... data structure with LZC values and n_valid parameters
%
% Implementation follows: Mediano et al. 2021 https://www.biorxiv.org/content/10.1101/2020.11.01.356071v1

% See Grandy et al., 2016 for justification of concatenating shorter snipets (as may be the case for data with intermitted nan sesctions). Grandy, T.H., Garrett, D.D., Schmiedek, F., and Werkle-Bergner, M. (2016). On the estimation of brain signal entropy from sparse neuroimaging data. Scientific Reports 6, 23073.
% For justification of *multiscale* Lempel-Ziv, see Ib??ez-Molina, Antonio J., et al. "Multism Lempel?Ziv complexity for EEG measures." Clinical Neurophysiology 126.3 (2015): 541-548.
%
% 
%
%

%parpool([2 16]) % parallel processing for time scales 

[n_chan,~]=size(dat);
dbstop if error

if ~isfield(cfg,'verbose'),            cfg.verbose=0;                   end % if true, plot progress in command line



%if cfg.verbose, fprintf('Compute Lempel Ziv Complexity (%i scales)\n',length(cfg.smoothing)), end

% preallocation
vLZC = nan(n_chan,1); 
gLZC = nan(n_chan,length(cfg.foi)); 
CTW = nan(n_chan,1); % CTW entropy rate (LZc alternative), from Mediano et al. 
CSER = nan(n_chan,1); % State space entropy rate (LZc alternative), from Mediano et al. 

n_valid_vLZC = nan(size(vLZC));
n_valid_gLZC = nan(size(gLZC));
n_valid_CTW = nan(size(CTW));
n_valid_CSER = nan(size(CSER));

vB = nan(size(vLZC));
gB = nan(size(vLZC));


for ichan=1:n_chan
    if cfg.verbose, tic, fprintf('Channel %i of %i ..',ichan,n_chan); end
    
    x=dat(ichan,:);
    idx_up = find(diff(isnan([1,x]))==1);
    idx_down = find(diff(isnan([1,x]))==-1);
    if length(idx_down)<length(idx_up), idx_down(end+1)=length(x); end
    for icnt=length(idx_up):-1:1
        if idx_up(icnt)~=idx_down(icnt)
            x(idx_up(icnt)+1:idx_down(icnt)-1)=[]; % shrink NaN sections to just one NaN (speeds things up)
        end
    end
    
    assert(size(x,1)==1,'Tdw variable contains more than 1 row')
    if ~isempty(x) % check for NaNs
        discard = isnan(x);
        x(discard) = []; % remove NaNs
    end
    
    if isempty(x) || length(x) < cfg.smooth_lo(end)
        vLZC(ichan) = NaN;
        n_valid_vLZC(ichan) = NaN;
        gLZC(ichan,:) = NaN;

        n_valid_gLZC(ichan,:) = NaN;

    else
        % "Vanilla" LZ (using median of whole signal)
        n_valid_vLZC(ichan)=length(x); % number of data points
        vB(ichan) = n_valid_vLZC(ichan)/log2(n_valid_vLZC(ichan));
        vLZC(ichan) = LZ76(x>=median(x))/vB(ichan);
        n_valid_CTW(ichan)=length(x); % number of data points
        CTW(ichan) = CTWEntropyRate(x); % CTW entropy rate, from Mediano et al. 
        n_valid_CSER(ichan)=length(x); % number of data points
        %CSER(ichan) = StateSpaceEntropyRate(x, cfg.new_srate); % state space entropy rate from Mediano et al. 
        
        for ism=1:length(cfg.foi)
            % make a copy of the signal
            refsig = x; % this will become the origina "reference" signal
            if isempty(refsig) || length(refsig) < cfg.smooth_lo(end) % if this data was all NaNs or shorter than the largest smoothing window
                gLZC(ichan,ism) = NaN;
                n_valid_gLZC(ichan,ism) = NaN;

            else
                % make binary sequence using moving-median filter
                % Ibanez-Molina et al 2015
                LPsig = medfilt1(refsig,cfg.smooth_hi(ism),'truncate'); % attenuate high frequencies
                Tdw   = medfilt1(refsig,cfg.smooth_lo(ism),'truncate');  % generate threshold for binarization
                
                % trim the beginning and end of the time series where median is
                % computed based on truncated data
                assert(length(Tdw)==length(LPsig),'Moving threshold doesn''t match signal length')
                start = (cfg.smooth_lo(ism)-1)/2+1;
                stop  = length(Tdw) - (cfg.smooth_lo(ism)-1)/2;
                LPsig = LPsig(start:stop);
                Tdw = Tdw(start:stop);
                
                % do one last check on data length, since we truncated the time
                % series
                
                if length(Tdw) <= 2 % Needs to be at least two elements (since diff operation removes one)
                    gLZC(ichan,ism) = NaN;
                    n_valid_gLZC(ichan,ism) = NaN;

                else
                    n_valid_gLZC(ichan,ism)=length(LPsig); % number of data points
                    
                    % normalization factors (to control for data length)
                    gB(ichan,ism) = n_valid_gLZC(ichan,ism)/log2(n_valid_gLZC(ichan,ism));
                    
                    % Compute gLZC below, credit Mediano et al for LZ algoirthm
                    % and Yeh and Shi 2018 for gLZC
                    gLZC(ichan,ism) =  LZ76(LPsig>=Tdw)/gB(ichan,ism);
                end
                clear Tdw refsig LPsig % just to be safe
            end
        end
     end
     if cfg.verbose, fprintf('. done (%.2f sec)\n',toc), end
 end

cfg.n_valid_vLZC = n_valid_vLZC;
cfg.n_valid_gLZC = n_valid_gLZC;
cfg.n_valid_CTW = n_valid_CTW;
cfg.n_valid_CSER = n_valid_CSER;

cfg.vLZC = vLZC;
cfg.gLZC = gLZC;
cfg.CTW = CTW;
cfg.CSER = CSER;

cfgout = cfg;

end


