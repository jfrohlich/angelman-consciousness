function[PEout] = CalcPermEn(dat,tau,chan)

% wrapper function for computing permutation entropy
% tau should be a scalar lag values in ms

load angelman_lay lay
%ft_defaults
assert(length(tau)==1,'Only one tau value should be used (must be scalar)')
cfg = [];
cfg.sf = 125;
cfg.verbose = false;
cfg.kernel = 3;
cfg.taus = round(tau*(cfg.sf/1000)); % covert tau from ms to samples (e.g., 32 ms becomes 8 samples at fs = 250 Hz)
cfg.taus_ms = tau;
cfg.chan_sel = find(contains(lay.label,chan)); 
% for some reason this crashes if we only select one channel, so we need to
% add a second channel
% if cfg.chan_sel < length(lay.label)-1;
%     cfg.chan_sel = [cfg.chan_sel cfg.chan_sel + 1]; 
% else
%     cfg.chan_sel = [cfg.chan_sel cfg.chan_sel - 1];
% end

cfg.data_sel = 1:size(dat,2); % use all samples
cfg_out = ro_PermEn(dat,cfg);

idx = strcmp(lay.label(cfg.chan_sel),chan); % Just return the channel of interest
PEout = cfg_out.PE(idx);
%PEout = cfg_out.PE;

end

