box off
set(gca,'linewidth',3)
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 20)
set(xAX,'color','k')
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 20)
set(yAX,'color','k')
set(gca, 'TickDir', 'out')
set(gcf,'color','w')
set(gca,'Layer','top')
axis normal
