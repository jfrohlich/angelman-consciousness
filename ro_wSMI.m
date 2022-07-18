function[wSMI_out] = ro_wSMI(X,cfg)


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
        wSMI_out.SMI{itau} = nan(n,n);
        wSMI_out.wSMI{itau} = nan(n,n);
        wSMI_out.PE{itau} = nan(n,1);
    end
    wSMI_out.cfg = cfg;
    return
end

SMI2 = cell(1,1);
wSMI2 = cell(1,1);
ind = find(tril(ones(n,n),-1)); % lower triangle indices

for itau = 1:length(sym)
    [SMI{itau}, wSMI{itau}, PE{itau}] = SMI_and_wSMI(sym{itau},count{itau},cfg.kernel); % Sitt et al. code
    SMI2{itau} = nan(n,n,size(SMI{itau},3));
    wSMI2{itau} = nan(n,n,size(wSMI{itau},3));
    for ilay = 1:size(SMI{itau},3)
        tmp = zeros(n);
        tmp(ind) = SMI{itau};
        SMI2{itau} = tmp + tmp';
        tmp = zeros(n);
        tmp(ind) = wSMI{itau};
        wSMI2{itau}(:,:,ilay) = tmp + tmp';
%         %%% DEBUG %%%
%         checkme = wSMI2{itau}(:,:,ilay);
%         if all(isnan(checkme(:)))
%             keyboard
%         end
%         %%% END DEBUG %%%
    end
end

wSMI_out.SMI = SMI2;
wSMI_out.wSMI = wSMI2;
wSMI_out.PE = PE;
wSMI_out.cfg = cfg;

end