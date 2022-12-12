function[out] = myround(num,option)

% This is a wrapper version of the native MATLAB round function that allows
% user to specify whether we should round towards nearest OOD integer
% (option = 'odd') or nearest EVEN integer (option = 'even'). Setting option to
% 'closest' just uses the native MATLAB round() function. Takes row vectors
% or scalars as input.

if nargin < 2
    option = 'closest';
end

assert(size(num,1) == 1);

out = nan(size(num)); % allocation

switch option
    case 'even'
        for i = 1:length(num)
            down = floor(num(i));
            if mod(down,2) == 0 % if rounding down gives up an even number
                out(i) = down;
            else
                out(i) = down+1; % round up instead! 
            end
        end
    case 'odd'
        for i = 1:length(num)
            down = floor(num(i));
            if mod(down,2) == 1 % if rounding down gives up an odd number
                out(i) = down;
            else
                out(i) = down+1; % round up instead! 
            end
        end
    case 'either'
        for i = 1:length(num)
            out(i) = round(num(i));
        end
    otherwise
        error('Option not detected')
end

end
        