function[wSMI_out] = ro_PermEn(X,cfg)

% This function only computes permutation entropy (no connectivity
% measures like wSMI)

if ~isfield(cfg,'taus')
    cfg.taus = 32;
end

if ~isfield(cfg,'k')
    cfg.kernel = 3;
end

if ~isfield(cfg,'chan_sel')
    cfg.chan_sel = 1:size(X,1); % include all channels
end

if ~isfield(cfg,'data_sel')
    cfg.data_sel = 1:size(X,2); % include all samples
end

n = length(cfg.chan_sel);
assert(isfield(cfg,'sf'),'Sampling rate missing')

%%% this computes the symbolic time series the output is:
%%% Sym: symbolic time series
%%% Count the number of symbols in each channel/trial

try
    [sym ,count] = S_Transf(X,cfg);
catch
    % If the signal can't be filtered for some reason
    for itau = 1:length(cfg.taus)
        wSMI_out.PE{itau} = nan(n,1);
    end
    wSMI_out.cfg = cfg;
    return
end

for itau = 1:length(sym)
    %%%% penalization symbols
    % 1 <-> 1 & 1 <-> 6
    % 2 <-> 2 & 2 <-> 5
    % 3 <-> 3 & 3 <-> 4
    % 4 <-> 4 & 4 <-> 3
    % 5 <-> 5 & 5 <-> 2
    % 6 <-> 6 & 6 <-> 1


    %%% permutation entropy calculation
    PE            = count{itau}.*log(count{itau}); 
    PE(isnan(PE)) = 0;
    PE            = -sum(PE,2);

    wSMI_out.PE = PE;
    wSMI_out.cfg = cfg;

end