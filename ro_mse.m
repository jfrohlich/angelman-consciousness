
function [mse,scale,n_valid,r] = ro_mse(dat,cfg)

% Computes multi scale entropy (if scales is set to 1 -> sample entropy)
%
% Modified from original function (ro_mse) 09/20/18 by Joel Frohlich for
% use with parallel computing toolbox and updated 02/21/19 to allow the
% dynamic tolerance (different r for each scale)
%
% Notes: 100 time points minimum are recommended at each time scale, 900
% preferred (Grandy et al., 2016 citing Pincus and Goldberger, 1994). For
% 20 time scales sampled at 200 Hz, this gives 10 sec as shortest possible
% window, 90 sec even better (if feasible with computing resources).
%
% Input:
% - dat ... [channels x samples], works with NaN sections
% - cfg ... struct, see script itself for default parameters
% - n_valid .. # of valid patterns that are compared (n_valid/fsampl -> net time used)
% Output:
% - mse ... [channels x scales] multiscale entropy
% - scales ... [scales] vector with scales (same as cfg.scales)
%
% Implementation follows:
% (1) cfg.type='Costa': Costa, M., Goldberger, A.L., and Peng, C.-K. (2005). Multiscale entropy analysis of biological signals. Phys Rev E Stat Nonlin Soft Matter Phys 71, 21906.
% (2) cfg.type='Xie':   Xie, H.-B., He, W.-X., and Liu, H. (2008). Measuring time series regularity using nonlinear similarity-based sample entropy. Physics Letters A 372, 7140?7146.
% See Grandy et al., 2016 for justification of concatenating shorter snipets (as may be the case for data with intermitted nan sesctions). Grandy, T.H., Garrett, D.D., Schmiedek, F., and Werkle-Bergner, M. (2016). On the estimation of brain signal entropy from sparse neuroimaging data. Scientific Reports 6, 23073.
%
% % Example:
% rng('default'), rng(1)
% x=randn(5,2^10); x(:,100:200)=nan; x(:,1010:1300)=nan;
% cfg = [];
% cfg.verbose=1;
% cfg.r = 0.15*nanstd(x,[],2);
% cfg.m = 2;
% cfg.scales=1:10;
% cfg.type='Costa';
% mse=ro_mse(x,cfg)
%
% % Example 2 (produce figure as in Xie et al. (Figure 1)
% rng('default'), rng(5)
% x=randn(1,100);
% r=0.01:0.01:1;
% clear cfg
% cfg.m=2;
% cfg.scales=1;
% for icnt=1:length(r)
%     cfg.r = r(icnt);
%     cfg.type='Costa';
%     mse_costa(icnt)=ro_mse(x,cfg);
%     cfg.type='Xie';
%     mse_xie(icnt)=ro_mse(x,cfg);
% end
% figure('Color','w')
% plot(r,mse_xie,'ro')
% hold on
% plot(r,mse_costa,'k*')
%
% Example 3
% load ./papers/multiscale_entropy/packsource.mat
% fsample=100;
% datasource(:,sampletype100hz~=1)=nan;
% x = datasource;
% cfg = [];
% cfg.verbose=1;
% cfg.r = 0.15*nanstd(x,[],2);
% cfg.m = 2;
% cfg.scales=1:20;
% cfg.type='Costa';
% [mse,scales,n_valid]=ro_mse(x,cfg)
% n_valid/fsample
%
%
% Relevant papers on mse and autism:
% - Bosl, W., Tierney, A., Tager-Flusberg, H., and Nelson, C. (2011). EEG complexity as a biomarker for autism spectrum disorder risk. BMC Medicine 9, 18.
%   - EEG, method from Xie, m=2, r=0.15*std, filter: 0.1 to 100, fsample=250, 20 scales, 20 sec segemnt / subject
% - Ghanbari, Y., Bloy, L., Edgar, J.C., Blaskey, L., Verma, R., and Roberts, T.P.L. (2015). Joint Analysis of Band-Specific Functional Connectivity and Signal Complexity in Autism. J Autism Dev Disord 45, 444?460.
%   - MEG, method from Costa et al, m=2, r=0.2*std, filter: 0.03-150 Hz, fsample=500, >40s mean ~80 sec of data, data where bandpass filtered! (5 ranges and broad band), 60 scales
% - Bosl, W.J., Loddenkemper, T., and Nelson, C.A. (2017). Nonlinear EEG biomarker profiles for autism and absence epilepsy. Neuropsychiatric Electrophysiology 3, 1.
%   - EEG, method from Xie, m=2, r=0.15*std, filter: 0.1 to 100, fsample=200, 20 scales, 30 sec segemnt / subject
% - Catarino, A., Churches, O., Baron-Cohen, S., Andrade, A., King, H. (2011). Atypical EEG complexity in autism spectrum conditions: a multiscale entropy analysis
%   - EEG, method from Costa, m=2, r= 0.15*std. 40 scale factos, N= 40 000 data points at Fs = 1000Hz
% 2.4.2017, Joerg Hipp

[n_chan,~]=size(dat);

if ~isfield(cfg,'verbose'),            cfg.verbose=0;                      end % if true, plot progress in command line
if ~isfield(cfg,'scales'),             cfg.scales = [1:20];                end % vector with scales to derive e.g. 1:10, if set to 1 -> sample entropy
if ~isfield(cfg,'m'),                  cfg.m=2;                            end % embedded dimonesion, default is 2 as in Costa et al., 2005
if ~isfield(cfg,'tolerance'),          cfg.tolerance = 0.2;                end % proportion of SD to use as state space radius
if ~isfield(cfg,'r'),                  cfg.r=cfg.tolerance*nanstd(dat,[],2);end % vector default is 0.2*std(dat,[],2); the shorter the dataset is, the larger the minimum tolerance r is needed (Xie et al, 2008)
if ~isfield(cfg,'type'),               cfg.type='Costa';                   end % 'Costa' (default), or 'Xie'
if ~isfield(cfg,'gpu'),                cfg.gpu= false;                     end % Optimize for GPU?
if ~isfield(cfg,'dynr')                cfg.dynr= false;                    end % dynamic tolerance (calcualte separetly for each scale)?

if cfg.verbose, fprintf('Compute MultiScale Entropy (%i scales, m=%i, r=%.2e)\n',length(cfg.scales),cfg.m,mean(cfg.r)), end
mse = nan(n_chan,length(cfg.scales)); % preallocation

% make gpu arrays
if cfg.gpu
    mse = gpuArray(mse);
end

for ichan=1:n_chan
    if cfg.verbose, tic, fprintf('Channel %i of %i ..',ichan,n_chan); end
    for iscale=1:length(cfg.scales)
        x=dat(ichan,:);
        idx_up = find(diff(isnan([1,x]))==1);
        idx_down = find(diff(isnan([1,x]))==-1);
        if length(idx_down)<length(idx_up), idx_down(end+1)=length(x); end
        for icnt=length(idx_up):-1:1
            if idx_up(icnt)~=idx_down(icnt)
                x(idx_up(icnt)+1:idx_down(icnt)-1)=[]; % shrink NaN sections to just one NaN
            end
        end
        if iscale>1
            if cfg.dynr && ichan == 1
                cfg.r = [cfg.r nan(size(cfg.r,1),1)];
            end
            kernel= ones(1,cfg.scales(iscale))/cfg.scales(iscale);
            % crude moving average filter that implements the coarse
            % graining procedure
            x=conv(x,kernel,'valid'); % check! alternatively: x=conv(dat',kernel,'valid')';
            x=x(:,1:cfg.scales(iscale):end);
        end
        Norig=length(x);       
        xMat = zeros(cfg.m+1,Norig-cfg.m);
        for i = 1:cfg.m+1
            xMat(i,:) = x(i:Norig-cfg.m+i-1);
        end
        xMat(:,isnan(mean(xMat))) = []; % remove NaNs
        n_valid=size(xMat,2);
        if cfg.gpu
            xMat = gpuArray(xMat); % transfer data to GPU
        end
        if cfg.dynr
            cfg.r(ichan,iscale) = gather(cfg.tolerance*nanstd(x));
            %fprintf('     Tolerance set to r = %1.2f\n',cfg.r(ichan,iscale))
            
        end
        switch cfg.type
            case 'Costa'
                dist_m = abs(pdist(xMat(1:cfg.m,:)','chebychev'));
                n_m = length(find(dist_m < cfg.r(ichan,iscale) & dist_m > 0));
                dist_mp1 = abs(pdist(xMat(1:cfg.m+1,:)','chebychev'));
                n_mp1 = length(find(dist_mp1 < cfg.r(ichan,iscale) & dist_mp1 > 0));
                mse(ichan,iscale)=log(n_m/n_mp1);
            case 'Xie'
                xMat_m=xMat(1:cfg.m,:); xMat_m=xMat_m-repmat(mean(xMat_m),[cfg.m,1]);
                xMat_mp1=xMat; xMat_mp1=xMat_mp1-repmat(mean(xMat_mp1),[cfg.m+1,1]);
                n_m=0;
                dist_m = abs(pdist(xMat_m(1:cfg.m,:)','chebychev'));
                n_m=n_m+sum(1./(1+exp((dist_m-0.5)/cfg.r(ichan,iscale))));
                n_mp1=0;
                dist_mp1 = abs(pdist(xMat_mp1(1:cfg.m+1,:)','chebychev'));
                n_mp1=n_mp1+sum(1./(1+exp((dist_mp1-0.5)/cfg.r(ichan,iscale))));
                mse(ichan,iscale)=log(n_m/n_mp1);
        end
    end
    if cfg.verbose, fprintf('. done (%.2f sec)\n',toc), end
end

scale = cfg.scales;

if cfg.gpu
    mse = gather(mse);
end

r = cfg.r;

end
