function[] = mypcolor(data)

% wrapper function for pcolor that doesn't chop off the last row and column
[m,n] = size(data);
data2 = nan(m+1,n+1);
data2(2:m+1,1:n) = data;
pcolor(flipud(data2))

end