% Convert GPU array back to regular array 

pth = './MSE_TD/';
files = dir(sprintf('%s*.mat',pth));

for ifile = 1:length(files)
    fprintf('%s\n',files(ifile).name)
    MSEout = load(sprintf('%s%s',pth,files(ifile).name),'MSEout');
    
    % recursively find the real structure
    while isfield(MSEout,'MSEout')
        MSEout = MSEout.MSEout;
    end
    
    if isa(MSEout.mse,'gpuArray')
        % get data out of GPU array
        MSEout.mse = gather(MSEout.mse);
    end
    
    save(sprintf('%s%s',pth,files(ifile).name),'MSEout');
end

%%

pth = './MSE_Feb/Xie/sleep';
files = dir(sprintf('%s*.mat',pth));

for ifile = 1:length(files)
    fprintf('%s\n',files(ifile).name)
    MSEout = load(sprintf('%s%s',pth,files(ifile).name),'MSEout');
    
    % recursively find the real structure
    while isfield(MSEout,'MSEout')
        MSEout = MSEout.MSEout;
    end
    
    if isa(MSEout.mse,'gpuArray')
        % get data out of GPU array
        MSEout.mse = gather(MSEout.mse);
    end
    
    save(sprintf('%s%s',pth,files(ifile).name),'MSEout');
end
%%

pth = './MSE_Feb/Xie/wake';
files = dir(sprintf('%s*.mat',pth));

for ifile = 1:length(files)
    fprintf('%s\n',files(ifile).name)
    MSEout = load(sprintf('%s%s',pth,files(ifile).name),'MSEout');
    
    % recursively find the real structure
    while isfield(MSEout,'MSEout')
        MSEout = MSEout.MSEout;
    end
    
    if isa(MSEout.mse,'gpuArray')
        % get data out of GPU array
        MSEout.mse = gather(MSEout.mse);
    end
    
    save(sprintf('%s%s',pth,files(ifile).name),'MSEout');
end