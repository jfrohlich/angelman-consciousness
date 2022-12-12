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

%gcp
%parfor ichan=1:n_chan
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


