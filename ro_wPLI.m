function[wpli,dwpli,foi] = ro_wPLI(data,fs,foi)

if nargin < 3
    % same frequencies from wavelet analysis, 8 bins per octave
    foi = 2.^[-1:0.125:5];
end
tic
Pxy = nan(size(data,3),size(data,1),size(data,1),length(foi));

for itrl = 1:size(data,3)
    for ich = 1:size(data,1)
        for jch = 1:size(data,1)
            if ich > jch % only do one half triangle
                % if there are no nans (artifacts)
                if ~any(isnan(data(ich,:,itrl))) && ~any(isnan(data(jch,:,itrl)))
                    % cross spectral density
                    [tmp,fout] = cpsd(data(ich,:,itrl),data(jch,:,itrl),fs*2,fs,foi,fs);
                    assert(all(round(fout,3) == round(foi,3)))
                    Pxy(itrl,ich,jch,:) = tmp;
                end
            end
        end
    end
    fprintf('     %2.1f%% complete\n',itrl/size(data,3)*100)
end

% reshape as trial x channel combinations x frequency
% Pxy = reshape(Pxy,size(Pxy,1),1,size(Pxy,2));

% compute wPLI using Fieldtrip
wpli = nan;
%wpli = ft_connectivity_wpli(Pxy,'debias',false,'dojack',false,'feedback','none');
dwpli = ft_connectivity_wpli(Pxy,'debias',true,'dojack',false,'feedback','none');

end

