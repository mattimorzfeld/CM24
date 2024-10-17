clear
close all
clc
load('myLSMOD2.mat')

Colors = brewermap(8,'Set1');
WhichColor = [3 2 5 1];
t = t-1.990;

% figure
% subplot(211)
% plot(t,L,'k')
% ylim([-90 90])
% xlim([t(1) t(end)])
% ylabel('Pole latitude')
% xlabel('Time (kyr)')
% set(gca,'FontSize',16)
% box off
% 
% subplot(212)
% plot(t,D,'k')
% xlim([t(1) t(end)])
% xlabel('Time (kyr)')
% ylabel('Dipole moment')
% set(gcf,'Color','w')
% set(gca,'FontSize',16)
% box off


figure
subplot(411)
plot(t,D,'k')
xlim([t(1) t(end)])
xlabel('Time (ka)')
ylabel({'Dipole moment';'(ZAm^2)'})
set(gcf,'Color','w')
set(gca,'FontSize',16)
box off

subplot(412)
plot(t,L,'k')
ylim([-90 90])
xlim([t(1) t(end)])
ylabel({'Pole latitude';'(degrees)'})
xlabel('Time (ka)')
set(gca,'FontSize',16)
box off

subplot(413)
plot(t,log10(a),'k')
xlim([t(1) t(end)])
xlabel('Time (ka)')
ylabel('log_{10}(P_{i_D})')
set(gcf,'Color','w')
set(gca,'FontSize',16)
box off

%% Search for events
aCat = [.5 3 15];
aThStart = 0.5;
aThStop = 0.035;

[start,stop]=FindEventsActivity(a,aThStart,aThStop);
Duration = t(stop)-t(start);
NumEvents = length(start);

NumCats = length(aCat)+1;
LCats = cell(NumCats,1);
DCats = cell(NumCats,1);
logaCats = cell(NumCats,1);
aCats = cell(NumCats,1);
IndsCats = cell(NumCats,1);
TinCat = zeros(NumCats,1);
AllInds = 1:length(D);

for kk=1:NumCats
    if kk==1
        inds = a<aCat(1);
    elseif kk==NumCats
        inds = a>aCat(end);
    else
        inds = ((a<aCat(kk)) & (a>aCat(kk-1)));
    end
    LCats{kk} = L(inds);
    DCats{kk} = D(inds);
    logaCats{kk} = loga(inds);
    aCats{kk} = a(inds);
    IndsCats{kk} = AllInds(inds);
    TinCat(kk) = length(LCats{kk})/length(D);
end

%% Label by categories
cat = zeros(length(D),1);
for kk=1:NumCats
    cat(IndsCats{kk}) = kk;
end
WhatCat = zeros(length(NumEvents),1);
CatDurations = cell(NumCats-1,1);
for kk=1:NumEvents
    stt = start(kk);
    stp = stop(kk);
    WhatCat(kk) = max(cat(stt:stp));
    CatDurations{WhatCat(kk)-1} = [CatDurations{WhatCat(kk)-1}, t(stp)-t(stt)];
end

subplot(411),hold on
plot(t(start:stop),D(start:stop),'Color',Colors(WhichColor(3),:),'LineWidth',2)
% ylim([-1.6 1.6])

subplot(412),hold on
plot(t(start:stop),L(start:stop),'Color',Colors(WhichColor(3),:),'LineWidth',2)
ylim([-95 95])
yticks([-90 -45 0 45 90])

subplot(413),hold on
plot(t(start:stop),loga(start:stop),'Color',Colors(WhichColor(3),:),'LineWidth',2)
plot([t(1) t(end)],log10(aThStart)*[1 1],'Color',Colors(4,:))
plot([t(1) t(end)],log10(aThStop)*[1 1],'Color',Colors(1,:))

subplot(414), hold on
stairs(t,cat-1,'k-')
stairs(t(start-1:stop),cat(start-1:stop)-1,'Color',Colors(WhichColor(3),:),'LineWidth',2)
xlim([t(1) t(end)])
xlabel('Time (ka)')
ylabel('Level')
set(gcf,'Color','w')
set(gca,'FontSize',16)
ylim([0 3])
yticks(0:3)
box off

f = gcf;
f.Position = [100 100 400 800];



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
plot(D(start),L(start),'o','Color',Colors(7,:),'LineWidth',4,'MarkerSize',10)
plot(D(stop),L(stop),'x','Color',Colors(8,:),'LineWidth',4,'MarkerSize',20)
set(gcf,'Color','w')
set(gca,'XScale','log')
set(gca,'FontSize',16)
f = gcf;
f.Position = [100 100 400 800];
xlabel('Dipole moment')
ylabel({'Pole latitude';'(degrees)'})



 t(start)
 t(stop) 
 t(start)-t(stop)