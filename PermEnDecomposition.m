% Code modified by JF from original code by Pedro Mediano
% See here https://github.com/pmediano/EntRate

%     The PermEnDecomposition function takes in two sets of time series
%     data (x1 and x2) and computes the average difference in permutation
%     entropy between the two datasets. The function also decomposes the
%     permutation entropy difference in spectral, phasic, and interaction
%     terms using phase randomization techniques. The function accepts an
%     optional argument nb_samples which specifies the number of surrogate
%     samples to use (default is 1000). The function returns a struct res
%     containing all elements of the decomposition, including the
%     permutation entropy for each dataset, the total permutation entropy
%     difference, and the spectral, phasic, and interaction components of
%     the permutation entropy difference. The code begins by checking the
%     input arguments and initializing some variables. It then sets the
%     number of workers to use for parallelization and defines some
%     function handles depending on whether the code is being run in Octave
%     or Matlab. The code then loops over the different lag values and
%     computes the permutation entropy for each dataset using the
%     CalcPermEn function. Next, the code performs the phase randomization
%     on the data, and computes the permutation entropy and its spectral
%     and phasic components using the PermEn function. The code then saves
%     the results to the res struct and returns it at the end of the
%     function.
% 

function [ res ] = PermEnDecomposition(x1, s1, x2, s2, nb_samples, new_fsample)
%% PermEnDECOMPOSITION Decomposes spectral, phasic, and interaction components
% of the PermEn difference between two datasets.
%
%   R = PermEnDECOMPOSITION(X1, X2), where X1,X2 are 1D cell arrays of 1-by-T time
%   series vectors, computes the average difference in PermEn(X1) - PermEn(X2) and
%   decomposes it in spectral, phasic, and interaction terms using phase
%   randomisation techniques.
%
%   R = PermEnDECOMPOSITION(X1, X2, N) uses N surrogate samples (default: 1000).
%
%   R = PermEnDECOMPOSITION(..., S) appends the results to pre-existing struct S.
%
% Results are returned in a struct R with all elements of the decomposition.
% These are:
%
%   - pe1: average PermEn of dataset X1
%   - pe2: average PermEn of dataset X2
%   - diff: total PermEn difference (X1 - X2)
%   - pai1: phase-amplitud component of PermEn in X1
%   - pai2: phase-amplitud component of PermEn in X2
%   - phase: phasic component of PermEn difference (X1 - X2)
%   - spec: spectral component of PermEn difference (X1 - X2)
%
% Reference:
%   Mediano, P. A., Rosas, F. E., Barrett, A. B., & Bor, D. (2020). Decomposing
%   spectral and phasic differences in non-linear features between datasets.
%   arXiv:2009.10015.
%   https://arxiv.org/pdf/2009.10015.pdf
%
% Adapted from code by Pedro Mediano, Oct 2019
% Modified by Joel Frohlich, August 2021

dbstop if error

assert(strcmp(s1,s2),'Sleep and wake are from different subjects')

%% Load channel info

load Angelman_lay lay

%% Argument checks and parameter initialisation
if nargin < 3 || isempty(nb_samples)
  nb_samples = 200;
end
if nargin < 4  || isempty(new_fsample)
    new_fsample = 125; 
end

nchan = size(x1{1},1);
assert(nchan==size(x2{1},1),'Number of channels doesn''t match between conditions')

x1_valid = iscell(x1) && all(cellfun(@length, x1) == length(x1{1}));
x2_valid = iscell(x2) && all(cellfun(@length, x2) == length(x2{1}));
if ~(x1_valid && x2_valid)
  error('Inputs must be 1D cell arrays of equally-sized trials');
end
nb_windows_1 = length(x1);
nb_windows_2 = length(x2);

PermEn = @(u,v,w) CalcPermEn(u,v,w);



isoctave = exist('OCTAVE_VERSION', 'builtin');
if isoctave
    % We're in Octave -- disable parallelisation by default
    nb_workers = 0;
else
    % We're in Matlab -- enable parfor only if not run in a worker thread already
    isOnWorker = ~isempty(getCurrentTask());
    if isOnWorker
        nb_workers = 0;
    else
        nb_workers = 32;
    end
end

% Set function handle for symmetric IFFT depending on Octave/Matlab version
if isoctave
    symm_ifft = @(z) real(ifft(z));
else
    symm_ifft = @(z) ifft(z, 'symmetric');
end

%% loop through lag values

lags = [8 16 32 64 128]; % lags for permutation entropy in seconds

pe1 = nan(length(lags), nchan);
pe2 = nan(length(lags), nchan);

pe1_f1  = nan([nb_samples, length(lags), nchan]);
pe1_f12 = nan([nb_samples, length(lags), nchan]);
pe2_f2  = nan([nb_samples, length(lags), nchan]);
pe2_f12 = nan([nb_samples, length(lags), nchan]);

progcnt = 0; % progress counter

for ich = 1:nchan
    for ilag = 1:length(lags)
        %% Pre-compute spectra, phases and true PermEn
        psd1  = cellfun(@(v)   abs(fft(v)), x1, 'UniformOutput', false);
        phi1  = cellfun(@(v) angle(fft(v)), x1, 'UniformOutput', false);
        psd2  = cellfun(@(v)   abs(fft(v)), x2, 'UniformOutput', false);
        phi2  = cellfun(@(v) angle(fft(v)), x2, 'UniformOutput', false);
        % phi_all = [phi1; phi2];
        phi_all = [phi1 phi2]; % JF edit - 07/02/21
        assert(all(size(x1)==size(x2)))
        tau = repmat({lags(ilag)},size(x1,1),size(x1,2));
        chan = repmat(lay.label(ich),size(x1,1),size(x1,2));
        pe1(ilag,ich) = mean(cellfun(PermEn, x1, tau, chan));
        pe2(ilag,ich) = mean(cellfun(PermEn, x2, tau, chan));


        %% Initialise arrays and compute surrogate PermEn values
        %{
            Notation:
              x1     -- actual time series of condition 1
              x1_f1  -- time series with spectra from condition 1 and phases swapped within condition 1
              x1_f12 -- time series with spectra from condition 1 and phases pooled from both conditions
        %}


        parfor (isamp=1:nb_samples, nb_workers)
        %for isamp = 1:nb_samples  
            % Swapped within condition
            phi1_shuf = phi1(randperm(nb_windows_1));
            phi2_shuf = phi2(randperm(nb_windows_2));
            x1_f1   = cellfun(@(s, p) symm_ifft(s.*exp(1j.*p)), psd1, phi1_shuf, 'UniformOutput', false);
            tau = repmat({lags(ilag)},size(x1_f1,1),size(x1_f1,2));
            chan = repmat(lay.label(ich),size(x1,1),size(x1,2));
            pe1_f1(isamp,ilag,ich) = mean(cellfun(PermEn, x1_f1, tau, chan));

            x2_f2   = cellfun(@(s, p) symm_ifft(s.*exp(1j.*p)), psd2, phi2_shuf, 'UniformOutput', false);
            tau = repmat({lags(ilag)},size(x2_f2,1),size(x2_f2,2));
            chan = repmat(lay.label(ich),size(x1,1),size(x1,2));
            pe2_f2(isamp,ilag,ich) = mean(cellfun(PermEn, x2_f2, tau, chan));

            % Swapped across conditions
            shuf_idx = randperm(nb_windows_1 + nb_windows_2);
            phi1_shuf = phi_all(shuf_idx(1:nb_windows_1));
            phi2_shuf = phi_all(shuf_idx(nb_windows_1 + 1:end));

            x1_f12  = cellfun(@(s, p) symm_ifft(s.*exp(1j.*p)), psd1, phi1_shuf, 'UniformOutput', false);
            tau = repmat({lags(ilag)},size(x1_f12,1),size(x1_f12,2));
            chan = repmat(lay.label(ich),size(x1,1),size(x1,2));
            pe1_f12(isamp,ilag,ich) = mean(cellfun(PermEn, x1_f12, tau, chan));

            x2_f12  = cellfun(@(s, p) symm_ifft(s.*exp(1j.*p)), psd2, phi2_shuf, 'UniformOutput', false);
            tau = repmat({lags(ilag)},size(x2_f12,1),size(x2_f12,2));
            chan = repmat(lay.label(ich),size(x1,1),size(x1,2));
            pe2_f12(isamp,ilag,ich) = mean(cellfun(PermEn, x2_f12, tau, chan));
            if isamp == 50
                fprintf('     Getting there ...\n')
            elseif isamp == 100
                fprintf('     Halfway there ...\n')
            elseif isamp == 150
                fprintf('     Almost there!\n')
            end
        end
        progcnt = progcnt + 1;
        %if mod(progcnt,2) == 0 % sets how often progress is reported
            fprintf('\n%2.1f%% complete for %s ...\n',...
                (progcnt/(length(lags)*nchan))*100, s1)
        %end
    end
end
%% Put results together and return

% Uncomment below to average before saving in the results structure
pe1_f1  = squeeze(mean(pe1_f1));
pe1_f12 = squeeze(mean(pe1_f12));
pe2_f2  = squeeze(mean(pe2_f2));
pe2_f12 = squeeze(mean(pe2_f12));


res.pe1 = pe1;
res.pe2 = pe2;
res.diff = pe1 - pe2;
res.pai1 = pe1 - pe1_f1;
res.pai2 = pe2 - pe2_f2;
phase1 = pe1_f1 - pe1_f12;
phase2 = pe2_f2 - pe2_f12;
res.phase  = phase1 - phase2;
res.spec = pe1_f12 - pe2_f12; % the phases are pooled, so difference must be due to spectra
res.interact = res.pai1 - res.pai2; 
res.lags = lags;


end

