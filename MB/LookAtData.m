clear
close all
clc
Colors = brewermap(8,'Set1');
WhichColor = [3 2 5 1];

aCat = [.5 3 15];
aThStart = 0.5; 
aThStop  = 0.035;
T = load('GGFMBdipolecoeffs.txt');
% time g10 g11 h11

t   =-T(:,1);
g10 = T(:,2);
g11 = T(:,3);
h11 = T(:,4);
 
D = 2.586*sqrt(g10.^2 + g11.^2 + h11.^2)*1e-3;
coLat = -atan2(sqrt(g11.^2 +h11.^2),g10);
Lrad = pi/2-coLat-pi;
L = Lrad*180/pi;
 
Dbar = 80;
a =  ((pi/2)-abs(Lrad))./(pi*D)*Dbar;
loga = log10(a);

%% get Level
lev = zeros(length(a),1);
for kk=1:length(a)
    if a(kk)<aCat(1)
        % nothing
    elseif a(kk)>=aCat(1) && a(kk)<aCat(2)
        lev(kk) = 1;
    elseif a(kk)>=aCat(2) && a(kk)<aCat(3)
        lev(kk) = 2;
    elseif a(kk)>=aCat(3)
        lev(kk) = 3;
    else
        error('Something is wrong')
    end
end

%% start and end of the event
stt = find(a>=aThStart,1,'last');
stp = find(a(1:stt)<=aThStop,1,'last');

%% Event duration
EvDur = (stp-stt)*(t(2)-t(1));
fprintf('Event duration: %g ka\n',EvDur)


% figure(102)
% h1=histogram(loga,10,'Normalization','pdf');
% figure(100),hold on
% stairs(h1.BinEdges(2:end),h1.Values,'Color',Colors(1,:),'LineWidth',2);
% set(gcf,'Color','w')
% set(gca,'FontSize',16)
% xlabel('log(PSV index)')
% ylabel('pdf')
% box off
% title('PSV Index')
% close(102)
% 
% figure(102)
% h1=histogram(loga(loga<log10(0.5)),10,'Normalization','pdf');
% figure(100),hold on
% stairs(h1.BinEdges(2:end),h1.Values,'Color',Colors(2,:),'LineWidth',2);
% set(gcf,'Color','w')
% set(gca,'FontSize',16)
% xlabel('log(PSV index)')
% ylabel('pdf')
% box off

disp(' ')

fprintf('Median (log) of PSV index: %g\n',10^median(loga))
fprintf('Mean (log) of PSV index: %g\n',10^mean(loga))

disp(' ')

fprintf('Median (log) of PSV index in normal state: %g\n',10^median(loga(loga<log10(0.5))))
fprintf('Mean (log) of PSV index in normal state: %g\n',10^mean(loga(loga<log10(0.5))))

% disp(' ')
% fprintf('Conclusion: 0.08 is a good candidate for a typical PSV index\n')


figure
subplot(411), hold on
plot(t,D,'k')
plot(t(stt:-1:stp),D(stt:-1:stp),'Color',Colors(WhichColor(3),:),'LineWidth',3) 
% plot([t(1) t(end)],[0 0],'k')
ylabel({'Dipole moment'; '(ZAm^2)'})
xlabel('Time (ka)')
set(gca,'Fontsize',16)
box off

subplot(412), hold on
plot(t,L,'k')
plot(t(stt:-1:stp),L(stt:-1:stp),'Color',Colors(WhichColor(3),:),'LineWidth',3) 
ylim([-95 95])
yticks([-90 -45 0 45 90])
ylabel({'Pole latitude';'(degrees)'})
xlabel('Time (ka)')
set(gca,'Fontsize',16)
box off

subplot(413), hold on
plot(t,loga,'k')
plot(t(stt:-1:stp),loga(stt:-1:stp),'Color',Colors(WhichColor(3),:),'LineWidth',3) 
plot([t(1) t(end)],log10(aThStart)*[1 1],'Color',Colors(4,:))
plot([t(1) t(end)],log10(aThStop)*[1 1],'Color',Colors(1,:))
ylabel('log_{10}(P_{i_D})')
xlabel('Time (ka)')
ylim([-3 2])
set(gca,'Fontsize',16)

subplot(414), hold on
% stais(t,lev,'k')
stairs(t(stt+1:-1:stp+1),lev(stt+1:-1:stp+1),'Color',Colors(WhichColor(3),:),'LineWidth',2) 
ylabel('Level')
xlabel({'Time'; '(ka)'})
ylim([0 3])
yticks(0:3)
xlim([-900 -700])
set(gca,'Fontsize',16)
box off
set(gcf,'Color','w')

f = gcf;
f.Position = [100 100 400 800];


% stt = find(a>=aThStart,1,'last');
% stp = find(a(1:stt)<=aThStop,1,'last')
% 
% figure, hold on
% plot(loga,'k')
% plot(loga(1:stt),'r')
% 
% plot(stt,loga(stt),'.','MarkerSize',20) 
% plot(stp,loga(stp),'.','MarkerSize',20) 
% plot([1 length(a)],log10(aThStart)*[1 1],'Color',Colors(4,:))
% plot([1 length(a)],log10(aThStop)*[1 1],'Color',Colors(1,:))
% ylabel('Act. ind. (log)')
% xlabel('Time (kyr)')
% ylim([-3 2])
% set(gca,'Fontsize',16)

lev0Inds = find(lev==0);
lev1Inds = find(lev==1);
lev2Inds = find(lev==2);
lev3Inds = find(lev==3);




figure, hold on
DGrid = logspace(0,2,1000);
LGrid = 0:.5:90;
PlotPSVIndexContours(aCat,LGrid,DGrid,inf,2)
PlotPSVIndexNegContours(aCat,LGrid,DGrid,inf,2)

plot(D,L,'k')
plot(D(lev0Inds),L(lev0Inds),'.','Color',Colors(WhichColor(1),:),'LineWidth',2,'MarkerSize',20)
plot(D(lev1Inds),L(lev1Inds),'.','Color',Colors(WhichColor(2),:),'LineWidth',2,'MarkerSize',20)
plot(D(lev2Inds),L(lev2Inds),'.','Color',Colors(WhichColor(3),:),'LineWidth',2,'MarkerSize',20)
plot([1e0 1e2],[0 0],'k')
plot(D(stt),L(stt),'o','Color',Colors(7,:),'LineWidth',4,'MarkerSize',10)
plot(D(stp),L(stp),'x','Color',Colors(8,:),'LineWidth',4,'MarkerSize',20)
set(gcf,'Color','w')
set(gca,'XScale','log')
set(gca,'FontSize',16)
f = gcf;
f.Position = [100 100 400 800];
xlabel('Dipole moment')
ylabel({'Pole latitude';'(degrees)'})




