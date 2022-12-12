%clear vars
load angelman_lay.mat lay

POS = lay.cfg.elec.chanpos; % 3D Cartesian coordinates
ED = nan(size(POS,1),size(POS,1)); % preallocation for Euclidian distances
for ich = 1:size(POS,1)
    for jch = 1:size(POS,1)
        ED(ich,jch) = sqrt((POS(ich,1)-POS(jch,1))^2 + (POS(ich,2)-POS(jch,2))^2 + ...
            (POS(ich,3)-POS(jch,3))^2);
    end
end

myfigure
grid on
hold on
stem3(POS(contains(lay.label,'F'),1),POS(contains(lay.label,'F'),2),POS(contains(lay.label,'F'),3),'r','fill')
stem3(POS(contains(lay.label,'C'),1),POS(contains(lay.label,'C'),2),POS(contains(lay.label,'C'),3),'k','fill')
stem3(POS(contains(lay.label,'T'),1),POS(contains(lay.label,'T'),2),POS(contains(lay.label,'T'),3),'g','fill')
stem3(POS(contains(lay.label,'O'),1),POS(contains(lay.label,'O'),2),POS(contains(lay.label,'O'),3),'m','fill')
stem3(POS(contains(lay.label,'P'),1),POS(contains(lay.label,'P'),2),POS(contains(lay.label,'P'),3),'y','fill')

for ich = 1:length(lay.label)
    text(POS(ich,1),POS(ich,2),POS(ich,3),lay.label{ich})
end
%%

figure('color','w')
mypcolor(ED), axis square, mycolorbar
xticks([1:19])
xticklabels(lay.label)
yticks([1:19])
yticklabels(lay.label)

figure
histogram(ED(:),50)
xlabel('Euclidian distance (mm)')
ylabel('Count')
makefighandsome
print('-dpng','./Figures/ChannelDistances.png')
print('-dsvg','./Figures/ChannelDistances.svg')

%%

source = ED == 0;
neighbors = ED <= 80 & ED > 0;
shortrange = ED > 80 & ED < 130;
longrange = ED >= 130;

myfigure2
mypcolor(neighbors+shortrange*2+longrange*3)
axis square
colormap gray
xticks([1.5:19.5])
xticklabels(lay.label)
yticks([1.5:19.5])
yticklabels(flipud(lay.label))
title('Channel pairs','fontsize',18)
makefighandsome, axis square
%set(yAX,'FontSize', 1)
%set(xAX,'FontSize', 10)
axis([1 20 1 20])

print('-dpng','./Figures/electrode_distances.png')
print('-dsvg','./Figures/electrode_distances.svg')