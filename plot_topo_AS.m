function [] = plot_topo_AS(dat,cfg)

if nargin<2; cfg = []; end

if ~isfield(cfg,'comment'); cfg.comment = 'no'; end

load angelman_lay.mat lay

cfg.parameter = 'avg';
cfg.xlim = [0,0];
cfg.layout = lay;
cfg.gridscale = 400;
cfg.comment = 'no';
cfg.markersize = 8;

if length(dat) == 18 % if banana montage
    newLocs = nan(18,2);
    A = {'FP1','F3','C3','P3','FP2','F4','C4','P4','FP1','F7','T3','T5','FP2',...
    'F8','T4','T6','FZ','CZ'};
    B = {'F3','C3','P3','O1','F4','C4','P4','O2','F7','T3','T5','O1','F8','T4', ...
    'T6','O2','CZ','PZ'};
    for ich = 1:18
        I = find(strcmp(upper(lay.label),A{ich}));
        J = find(strcmp(upper(lay.label),B{ich}));
        newLocs(ich,:) = mean([lay.pos(I,:); lay.pos(J,:)]);
    end
    lay.pos = newLocs;
    lay.width = lay.width(1:18);
    lay.height = lay.height(1:18);
    lay.label = strcat(A,'-',B); % cell array of new labels 
    display('18 channels detected, assuming bipolar banana montage')
end

cfg.parameter = 'avg';
cfg.zlim = 'maxabs';%[0,0];
cfg.layout = lay;

dummy.avg = dat;
dummy.var = zeros(size(dat));
dummy.dof = ones(size(dat));
dummy.time = 0;
dummy.dimord = 'chan_time';
dummy.label = lay.label;


ft_topoplotER_JF(cfg,dummy);
