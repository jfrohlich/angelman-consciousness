function [signal_pe nb symbols] = PE(signal,kernel,tau)
% results = pe(signal,kernel)
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% (c) Jacobo Sitt
% ---------------------------------------------------------------------

% make sure first dimension is used
if length(size(signal)) == 1 && size(signal,1) == 1, signal = signal'; end

% identify symbols from all possible permutations
symbols = perms(1:kernel);
code_symbols = symbols(:,1)*9+symbols(:,2)*3+symbols(:,3);

%%% with this transformation the code symbols are:
% 321 = 9 x 3 + 3 x 2 + 1 = 34 , pos 1
% 312 = 9 x 3 + 3 x 1 + 2 = 32 , pos 2
% 231 = 9 x 2 + 3 x 3 + 1 = 28 , pos 3
% 213 = 9 x 2 + 3 x 1 + 3 = 24 , pos 4
% 132 = 9 x 1 + 3 x 3 + 2 = 20 , pos 5
% 123 = 9 x 1 + 3 x 2 + 3 = 18 , pos 6

%%% initialize the ordered signal (shape = kernel x channels x size(signal,1)-tau*(kernel-1)
order_signal = zeros(kernel,size(signal,2),size(signal,1)-tau*(kernel-1));

for k = (size(signal,1)-tau*(kernel-1)):-1:1
   
[unused order_signal(:,:,k)] = sort([signal(k,:) ;signal(k+tau,:); signal(k+2*tau,:)],1,  'descend');
  
end

%%% code symbols with the same transformation as before
order_signal2 = squeeze(order_signal(1,:,:))*9+squeeze(order_signal(2,:,:))*3+squeeze(order_signal(3,:,:));

%%% search for each symbols
signal_pe = zeros(size(order_signal2)); % initialize the simbolic timeseries
n = 0;
for s = code_symbols'
    n =n +1;
    ii = find(order_signal2==s); % find symbol position in data
    signal_pe(ii) = n; % generate the signal in the symbolic space
end

%for i = 1:size(signal_pe,1)
for i = 1:min(size(signal_pe)) %-JF edit 09/01/21, this fixes the bug that crashes if data is just one channel
    for s = 1:6
        if min(size(signal_pe)) > 1
            nb(i,s) = length(find(signal_pe(i,:)==s));
        else
            nb(i,s) = length(find(signal_pe==s));
        end
    end
end




end
