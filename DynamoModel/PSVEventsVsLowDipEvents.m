clear
close all
clc
Colors = brewermap(8,'Dark2');


load('./DataFiles/myData.mat');
aThStart = 0.5;
aThStop = 0.035;

%% PSV events
[stPSV,spPSV]=FindEventsActivity(a,aThStart,aThStop);
NPSV = length(stPSV);
fprintf('Number of PSV events: %g\n',NPSV)
% find reversals
RevIndPSV = FindRevs(stPSV,spPSV,g10);
fprintf('Number of reversals: %g\n',sum(RevIndPSV))
% reversal rate
revRate = sum(RevIndPSV)/(t(end)*1e-3);
fprintf('Rev. rate: %g\n',revRate)
disp(' ')


%% Low dipole events
[stD,spD]=FindLowDipEvents(D,0.1,0.8);
Nld = length(stD);
fprintf('Number of low dipole events: %g\n',Nld)
% find reversals
RevIndD = FindRevs(stD,spD,g10);
fprintf('Number of reversals: %g\n',sum(RevIndD))
% reversal rate
revRate = sum(RevIndD)/(t(end)*1e-3 );
fprintf('Rev. rate: %g\n',revRate)

figure, hold on
plot(t*1e-3,loga,'k')
plot([t(1) t(end)]*1e-3,median(loga)*[1 1],'Color',Colors(2,:),'LineWidth',2)
plot([t(1) t(end)]*1e-3,mean(loga)*[1 1],'Color',Colors(3,:),'LineWidth',2)
xlim([0,t(end)*1e-3])
set(gcf,'Color','w')
set(gca,'FontSize',16)
xlabel('Time (Myr)')
ylabel('log(PSV index)')
box off


figure(102)
h1=histogram(loga,'Normalization','pdf');
figure(100),hold on
stairs(h1.BinEdges(2:end),h1.Values,'Color',Colors(1,:),'LineWidth',2);
set(gcf,'Color','w')
set(gca,'FontSize',16)
xlabel('log(PSV index)')
ylabel('pdf')
box off
title('PSV Index')
close(102)

figure(102)
h1=histogram(loga(loga<log10(0.5)),'Normalization','pdf');
figure(100),hold on
stairs(h1.BinEdges(2:end),h1.Values,'Color',Colors(2,:),'LineWidth',2);
set(gcf,'Color','w')
set(gca,'FontSize',16)
xlabel('log(PSV index)')
ylabel('pdf')
box off
legend('All','Normal state')
close(102)

disp(' ')

fprintf('Median (log) of PSV index: %g\n',10^median(loga))
fprintf('Mean (log) of PSV index: %g\n',10^mean(loga))

disp(' ')

fprintf('Median (log) of PSV index in normal state: %g\n',10^median(loga(loga<log10(0.5))))
fprintf('Mean (log) of PSV index in normal state: %g\n',10^mean(loga(loga<log10(0.5))))

disp(' ')
fprintf('Conclusion: 0.035 is a good candidate for a typical PSV index\n')