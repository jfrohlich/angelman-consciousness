box off
axis square
set(gca,'linewidth',5)
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 30)
set(xAX,'color','k')
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 30)
set(yAX,'color','k')
set(gca, 'TickDir', 'out')
set(gcf,'color','w')
set(gca,'Layer','top')
%xAx.Color = 'k';