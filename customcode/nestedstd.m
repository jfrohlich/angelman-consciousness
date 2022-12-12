function[sd] = nestedstd(x,index,dim)

% This function will compute the standard deviation along dim for nested data (e.g., multiple
% observations for some subjects) given the indicies for unique subjects, after averaging within subjects. 

% output args: sd is the stardard dev after averaging within subjects

if nargin < 3
    dim = 1;
end

% Convert cell array if necessary
if iscell(index)
    index2 = nan(1,length(index));
    for i = 1:length(index)
        index2(i) = str2num(index{i});
    end
    index = index2;
end


U = unique(index);
mu1 = [];
    
for i = 1:length(U)
    where = index == U(i);
    switch dim
        case 1
            data = x(where,:,:,:); % handled matrices up to 4d
            datamean = nanmean(data,dim); % average of that's subjects data
            mu1 = cat(dim,mu1,datamean);
        case 2
            data = x(:,where,:,:); % handled matrices up to 4d
            datamean = nanmean(data,dim); % average of that's subjects data
            mu1 = cat(dim,mu1,datamean);
        case 3
            data = x(:,:,where,:); % handled matrices up to 4d
            datamean = nanmean(data,dim); % average of that's subjects data
            mu1 = cat(dim,mu1,datamean);
        case 4
            data = x(:,:,:,where); % handled matrices up to 4d
            datamean = nanmean(data,dim); % average of that's subjects data
            mu1 = cat(dim,mu1,datamean);
        otherwise
            error('Dimension > 4 not allowed')
    end
end

try
    assert(size(mu1,dim)==length(U))
catch
    keyboard
end

sd = nanstd(mu1,[],dim);
    
    