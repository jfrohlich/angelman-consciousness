function[cohensd] = cohen_d(x1,x2)

% Assumes paired samples data 

n1 = length(x1);
n2 = length(x2);

assert(n1==n2,'unpaired samples')

pooled_sd = sqrt( ((n1-1)*std(x1)^2 + (n2-1)*std(x2)^2) / (n1+n2-1) );
cohensd = (mean(x1)-mean(x2)) / pooled_sd;

end