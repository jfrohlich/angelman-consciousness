function [SMI_out, wSMI_out, PE] = SMI_and_wSMI(signal,sym_prob,kernel)

% ---------------------------------------------------------------------
% Data should be on the symbolic space
% - signal should be a matrix of (channels x samples)
% - sym_prob is the probability of the symbols in each channel 
%  should be a matrix in the form of (channels x number of symbols)
% output:
% Symbolic mutual information, wSMI
% ---------------------------------------------------------------------
% (c) Jacobo Sitt
% ---------------------------------------------------------------------

channels = size(signal,1);
samples  = size(signal,2);

n_symbols = factorial(kernel);

%%%% intialize output variables
wSMI_out   = zeros(channels*(channels-1)/2,1);
SMI_out    = zeros(channels*(channels-1)/2,1);

joint_prob = zeros(n_symbols^2,channels*(channels-1)/2);


%%%% penalization symbols
% 1 <-> 1 & 1 <-> 6
% 2 <-> 2 & 2 <-> 5
% 3 <-> 3 & 3 <-> 4
% 4 <-> 4 & 4 <-> 3
% 5 <-> 5 & 5 <-> 2
% 6 <-> 6 & 6 <-> 1

forb = [6 5 4 3 2 1];

%%% permutation entropy calculation
PE            = sym_prob.*log(sym_prob); % bug here??? JF 13.05.21
PE(isnan(PE)) = 0;
PE            = -sum(PE,2);

index = 0;
for ch1 = 1:(channels-1) %%% from channel 
    
    for ch2 = (ch1+1):channels %%% to channel 
 
        index = index+1;
       
        %%% create a meta_signal using both signals
        meta_signal = signal(ch1,:) + (signal(ch2,:)-1)*n_symbols; 
        
        %%% compute joint probabilities (shape: n_symbolsˆ2 x channel pairs)
        for xxx = 1:n_symbols^2
            joint_prob(xxx,index) = sum(meta_signal == xxx)/samples;
        end
        
        
        out(index) = 0;
        
        for symbols_ch1 = 1:n_symbols
            for symbols_ch2 = 1:n_symbols
                
                meta_symbol = symbols_ch1 + (symbols_ch2-1)*n_symbols;
                
                %%% establish weight
                if symbols_ch1 == symbols_ch2 || forb(symbols_ch1) == symbols_ch2
                    w = 0;
                else
                    w = 1;
                end
                
                
                if joint_prob(meta_symbol,index) > 0 
                    
                    if sym_prob(ch1,symbols_ch1) == 0 || sym_prob(ch2,symbols_ch2) == 0
                        error('error: joint prob>0 but single prob == 0')
                    end
                    
                    %%%% Mutual information formula 
                    aux = joint_prob(meta_symbol,index)*...
                        log(joint_prob(meta_symbol,index)/sym_prob(ch1,symbols_ch1)/sym_prob(ch2,symbols_ch2));
                    
                    wSMI_out(index) = wSMI_out(index) + w*aux;
                    SMI_out(index)  =  SMI_out(index) +   aux;
                    
                
                end
                
 
                
            end
        end
        
    end
    
end


%%%% normalization btw 0 and 1
   wSMI_out = wSMI_out /  log(n_symbols) ;
   SMI_out = SMI_out /  log(n_symbols) ;


