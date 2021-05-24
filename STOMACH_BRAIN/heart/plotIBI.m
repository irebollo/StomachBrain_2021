
function plotIBI(ecg_n,cr,p,v,thresh,IBI_s)
figure(47894);clf
wholescreen = get(0,'ScreenSize');
pp = wholescreen;
pp(2) = wholescreen(4)/2;
pp(4) = wholescreen(4)/2;
set(gcf,'position',pp)

subplot(2,1,1)
plot(ecg_n,'b');
hold on
plot(cr(1001:end-1000),'r')
scatter(p,v)

axis tight
hline(thresh,'linestyle','--','color','k')
xlabel('samples')
ylabel('a.u.')
legend('ECG','corr','R-peak','location','eastoutside')
zoom
title('zoom in to check signal')

hhist = subplot(2,2,3);
hist(IBI_s,30)
xlabel('IBI (s)')
set(gca,'FontSize',14)
set(hhist,'tag','hhist')