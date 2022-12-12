function [pow,pow_full,pow_median,pow_var,n,unit,foi_target,foi_delta,foi,...
    cfg,w_lin2log] = ro_freq_wavelet_JF(dat,cfg)

% Morlet wavelet transform of continous electrophysiological data
%
% Input
% - dat ... [channels x samples], set invalid data sections to NaN assuming the input signal was in uV
% - cfg ... struct, see script itself for default parameters, cfg.fsample
%           has to be set.
% Output
% 
% Pow = output power spectrum. The total power in uV^2 for channel ich over a 
%   frequency range [f1 f2] can be obtained with idf = foi> f1 & foi<f2 and
%       if cfg.norm is 'linear',   trapz(foi(idf),pow(ich,idf))
%       if cfg.norm is 'log2',   trapz(log2(foi(idf)),pow(ich,idf))

% 9.7.2016, Joerg Hipp
% 30.5.2017, Joerg Hipp, return foi_target instead of foi

% senity checks
if ~isfield(cfg,'fsample'), error('specific cfg.fsample!'), end
n_sens=size(dat,1); n_sample=size(dat,2);
if n_sens>n_sample, fprintf('Warning, number of channles (%i) larger than number of samples (%i), check orientation of data matrix\n',n_sens,n_sample), end

if ~isfield(cfg,'verbose'),            cfg.verbose=1;                        end % if true, plot progress in command line
if ~isfield(cfg,'oct_bw'),             cfg.oct_bw = 0.5;                     end % wavelet frequency resolution -> q=f/sigma_f=
if ~isfield(cfg,'oct_delta'),          cfg.oct_delta = cfg.oct_bw/4;         end % wavelet frequency resolution -> q=f/sigma_f=
if ~isfield(cfg,'foi_start'),          cfg.foi_start = 2;                    end % wavelet frequency resolution -> q=f/sigma_f=
if ~isfield(cfg,'foi_end'),            cfg.foi_end = 32;                     end % wavelet frequency resolution -> q=f/sigma_f=
if ~isfield(cfg,'window_shift'),       cfg.window_shift=0.25;                end % fraction the window is shifted of total window size
if ~isfield(cfg,'kernel_width'),       cfg.kernel_width=5;                   end
if ~isfield(cfg,'allow_fraction_nan'), cfg.allow_fraction_nan=0;             end
if ~isfield(cfg,'norm');               cfg.norm='linear';                    end % 'linear', 'log2'
if ~isfield(cfg,'freq_shift_factor');  cfg.freq_shift_factor=1;              end

%% spectral parameter
foi        = 2.^[log2(cfg.foi_start):cfg.oct_delta:log2(cfg.foi_end)]; % [Hz] frequencies to analyze
foi_target = foi;
foi        = foi*cfg.freq_shift_factor;
foi_tmp    = 2.^[(log2(cfg.foi_start)-cfg.oct_delta/2):cfg.oct_delta:(log2(cfg.foi_end)+cfg.oct_delta/2)];
foi_tmp    = foi_tmp*cfg.freq_shift_factor;
foi_delta  = [foi_tmp(1:end-1);foi_tmp(2:end)]';          % [Hz] range for integration
foi_min    = 2*foi/(2^cfg.oct_bw+1);                      % arithmetic mean
foi_max    = 2*foi/(2^-cfg.oct_bw+1);
sigma_freq = (foi_max-foi_min)/(2*sqrt(2*log(2)));        % std in freq domain
sigma_time = 1./(2*pi*sigma_freq);                        % std in time domain

if cfg.verbose, tic, fprintf('Morelet Wavelet transform ['), end
for ifoi=1:length(foi)
    if cfg.verbose, fprintf('.'), end
    % convolution kernel
    n_win = round(cfg.kernel_width*sigma_time(ifoi)*cfg.fsample+1);
    n_shift = round(n_win*cfg.window_shift);
    t = ((1:n_win)-n_win/2-0.5)/cfg.fsample;
    z = t./sigma_time(ifoi);
    TAPER = exp(-(1/2)*z.^2);
    TAPER = TAPER/sqrt(sum(abs(TAPER).^2));
    iEXP = exp(1i*2*pi*foi(ifoi)*t);
    KERNEL = (TAPER.*iEXP).';
    
    % collect info on nan sections
    idx_up = find(diff([0,isnan(sum(dat))])==1);
    idx_down = find(diff([0,isnan(sum(dat))])==-1);
    nan_width = zeros(1,size(dat,2));
    for icnt=1:length(idx_up)-1
        nan_width(idx_up(icnt):idx_down(icnt)) = idx_down(icnt)-idx_up(icnt);
    end
    if length(idx_up)>length(idx_down)
        nan_width(idx_up(end):end) = length(nan_width)+1-idx_up(end);
    end
    
    % memory allocation
    DAT      = nan(n_sens,length(1:n_shift:n_sample-n_win+1));
    frac_nan = nan(1,size(DAT,2));
    % convolution
    cnt = 0;
    for isection = 1:n_shift:size(dat,2)-n_win+1
        section = double(dat(:,isection:isection+n_win-1));
        nan_width_section = nan_width(isection:isection+n_win-1);
        n_nan = sum(isnan(section(1,:)));
        cnt=cnt+1;
        frac_nan(cnt) = n_nan/size(section,2);
        if n_nan==0
            DAT(:,cnt) = section*KERNEL*sqrt(2/cfg.fsample);
        elseif n_nan<size(section,2)*cfg.allow_fraction_nan & ...
                max(nan_width_section)<size(section,2)*cfg.allow_fraction_nan
            idx_valid = find(~isnan(section(1,:)));
            KERNEL_tmp = KERNEL(idx_valid)/sqrt(sum(abs(TAPER(idx_valid)).^2));
            DAT(:,cnt) = section(:,idx_valid)*KERNEL_tmp*sqrt(2/cfg.fsample);
        else
            DAT(:,cnt) = nan(size(section,1),1);
        end % valid section
    end
    
    % derive measures for frequency-transformed data
    idx_valid = find(~isnan(DAT(1,:)));
    n_valid = length(idx_valid);
    DAT = DAT(:,idx_valid);
    frac_nan = frac_nan(idx_valid);
    
    % power
    pow(:,ifoi)        = mean(abs(DAT).^2,2);                % [uV^2*Hz^-1] assuming the input signal was in uV
    pow_full(:,ifoi)   = mean(abs(DAT(:,frac_nan==0)).^2,2);
    pow_median(:,ifoi) = median(abs(DAT).^2,2);
    pow_var(:,ifoi)    = var(abs(DAT).^2,[],2);
    n(ifoi)            = n_valid;
end % loop foi
if cfg.verbose, fprintf('] %.1f sec\n',toc), end
unit='ÂµVÂ²/Hz'; % assuming the input signal was in uV

cfg.norm

% Normalization
if strcmp(cfg.norm,'log2')
    w=diff(foi_delta')./cfg.oct_delta;
    pow = repmat(w,[n_sens,1]).*pow; % [uV^2/log2(Hz)] assuming the input signal was in uV
    pow_full = repmat(w,[n_sens,1]).*pow_full;
    pow_median = repmat(w,[n_sens,1]).*pow_median;
    pow_var = repmat(w,[n_sens,1]).*pow_var;
    unit='µV²/log2(Hz)'; % assuming the input signal was in uV
end
w_lin2log=diff(foi_delta')./cfg.oct_delta;

