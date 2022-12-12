function surr_y = surrogate_data(x,fs)

% Wrapper function to generate surrogates for normalizing complexity measures 

if size(x,1) > 1
	error('not a scalar time series')
end

% Do surrogate data for each block of valid data without NaNs

valid = x(~isnan(x)); % remove all NaNs! 

if isempty(valid) % No valid data
    surr_y = nan(size(x));
else
    y = valid;
    % Calls the function that Daniel gave me
    surr_y = surrogate(y, 1, 'FT', 0, fs);
end

end
