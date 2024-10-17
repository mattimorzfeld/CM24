clear
close all
clc

Colors = brewermap(8,'Set1');
WhichColor = [3 2 5 1];

load('./DataFiles/myData.mat');
absL = abs(L);
aCat = [.5 3 15];
aThStart = 0.5;
aThStop = 0.035;

%% Find all events
[start,stop]=FindEventsActivity(a,aThStart,aThStop);
NumEvents = length(start);
fprintf('Number of PSV events: %g\n',NumEvents)

%% Find category of each event
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

%% Assing events catergories to indeces
Cat1Inds = find(WhatCat-1==1);
Cat2Inds = find(WhatCat-1==2);
Cat3Inds = find(WhatCat-1==3);

%% find the cat-1 event that reverses
for oo=1:length(Cat1Inds)
    ind = Cat1Inds(oo);
    st_tmp = start(ind);
    sp_tmp  = stop(ind);
    if Ds(st_tmp)*Ds(sp_tmp)<0
        fprintf('Cat-1 event number %g reverses\n',oo)
    end
end

%% select a Cat-1 event
ind = Cat1Inds(608);
st = start(ind);
sp =  stop(ind);
stLong = st-1e2;
spLong = sp+1e2;

figure
subplot(411), hold on
plot(t(stLong:spLong)-t(st),D(stLong:spLong),'k','LineWidth',1)
plot(t(st:sp)-t(st),D(st:sp),'Color',Colors(WhichColor(2),:),'LineWidth',3)
ylabel({'Dipole moment'})
xlabel('Time (ka)')
xlim([t(stLong)-t(st) t(spLong)-t(st)])
ylim([0 1.6])
set(gca,'Fontsize',16)
% title('Cat-1')
box off

subplot(412), hold on
plot(t(stLong:spLong)-t(st),L(stLong:spLong),'k','LineWidth',1)
plot(t(st:sp)-t(st),L(st:sp),'Color',Colors(WhichColor(2),:),'LineWidth',3)

ylabel({'Pole latitude'; '(degrees)'})
xlabel('Time (ka)')
xlim([t(stLong)-t(st) t(spLong)-t(st)])
ylim([-95 95])
yticks([-90 -45 0 45 90])
set(gca,'Fontsize',16)
box off

subplot(413), hold on
plot(t(stLong:spLong)-t(st),loga(stLong:spLong),'k','LineWidth',1)
plot(t(st:sp)-t(st),loga(st:sp),'Color',Colors(WhichColor(2),:),'LineWidth',3)
plot([t(stLong)-t(st) t(spLong)-t(st)],log10(aThStart)*[1 1],'Color',Colors(4,:))
plot([t(stLong)-t(st) t(spLong)-t(st)],log10(aThStop)*[1 1],'Color',Colors(1,:))
ylabel('log_{10}(P_{i_D})')
xlabel('Time (ka)')
xlim([t(stLong)-t(st) t(spLong)-t(st)])
ylim([-5 4])
set(gca,'Fontsize',16)
box off
set(gcf,'Color','w')

subplot(414), hold on
stairs(t(stLong:spLong)-t(st),cat(stLong:spLong)-1,'k','LineWidth',1)
stairs(t(st-1:sp)-t(st),cat(st-1:sp)-1,'Color',Colors(WhichColor(2),:),'LineWidth',3)
ylabel('Level')
xlabel('Time (ka)')
xlim([t(stLong)-t(st) t(spLong)-t(st)])
ylim([0 3])
yticks(0:3)
set(gca,'Fontsize',16)
box off
set(gcf,'Color','w')

f = gcf;
f.Position = [100 100 400 800];

sc = 20;
transp = 1;
figure(21), hold on
for kk=1:length(aCat)+1
    s1=scatter(DCats{kk},LCats{kk},sc,'filled','MarkerEdgeColor',Colors(WhichColor(kk),:),'MarkerFaceColor',Colors(WhichColor(kk),:));
    s1.MarkerFaceAlpha = transp;
    s1.MarkerEdgeAlpha = transp;
end

DGrid = logspace(-4,0,1000);
LGrid = 0:.5:90;
figure(21)
PlotPSVIndexContours(aCat,LGrid,DGrid,inf,2)
PlotPSVIndexNegContours(aCat,LGrid,DGrid,inf,2)
set(gcf,'Color','w')
set(gca,'FontSize',16)
xlabel('Dipole moment')
ylabel({'Pole latitude'; '(degrees)'})
set(gca,'XScale','log')


figure(21)
hold on, plot(D(st:sp),L(st:sp),'-','Color',Colors(6,:),'LineWidth',2)
hold on, plot(D(st),L(st),'o','Color',Colors(7,:),'LineWidth',4,'MarkerSize',10)
hold on, plot(D(sp),L(sp),'x','Color',Colors(8,:),'LineWidth',4,'MarkerSize',20)
plot([1e-4 4.5],[0 0],'k')
xlim([1e-3 2.5])

f = gcf;
f.Position = [100 100 400 800];

%% select Cat-2 event
ind = Cat2Inds(24); %41
st = start(ind);
sp =  stop(ind);
stLong = st-1e2;
spLong = sp+1e2;

figure
subplot(411), hold on
plot(t(stLong:spLong)-t(st),D(stLong:spLong),'k','LineWidth',1)
plot(t(st:sp)-t(st),D(st:sp),'Color',Colors(WhichColor(3),:),'LineWidth',3)
ylabel({'Dipole moment'})
xlabel('Time (ka)')
xlim([t(stLong)-t(st) t(spLong)-t(st)])
ylim([0 1.6])
set(gca,'Fontsize',16)
box off

subplot(412), hold on
plot(t(stLong:spLong)-t(st),L(stLong:spLong),'k','LineWidth',1)
plot(t(st:sp)-t(st),L(st:sp),'Color',Colors(WhichColor(3),:),'LineWidth',3)

ylabel({'Pole latitude'; '(degrees)'})
xlabel('Time (ka)')
xlim([t(stLong)-t(st) t(spLong)-t(st)])
ylim([-95 95])
yticks([-90 -45 0 45 90])
set(gca,'Fontsize',16)
box off

subplot(413), hold on
plot(t(stLong:spLong)-t(st),loga(stLong:spLong),'k','LineWidth',1)
plot(t(st:sp)-t(st),loga(st:sp),'Color',Colors(WhichColor(3),:),'LineWidth',3)
plot([t(stLong)-t(st) t(spLong)-t(st)],log10(aThStart)*[1 1],'Color',Colors(4,:))
plot([t(stLong)-t(st) t(spLong)-t(st)],log10(aThStop)*[1 1],'Color',Colors(1,:))
ylabel('log_{10}(P_{i_D})')
xlabel('Time (ka)')
xlim([t(stLong)-t(st) t(spLong)-t(st)])
ylim([-5 4])
set(gca,'Fontsize',16)
box off
set(gcf,'Color','w')

subplot(414), hold on
stairs(t(stLong:spLong)-t(st),cat(stLong:spLong)-1,'k','LineWidth',1)
stairs(t(st-1:sp)-t(st),cat(st-1:sp)-1,'Color',Colors(WhichColor(3),:),'LineWidth',3)
ylabel('Level')
xlabel('Time (ka)')
xlim([t(stLong)-t(st) t(spLong)-t(st)])
ylim([0 3])
yticks(0:3)
set(gca,'Fontsize',16)
box off
set(gcf,'Color','w')

f = gcf;
f.Position = [100 100 400 800];

sc = 20;
transp = 1;
figure(22), hold on
for kk=1:length(aCat)+1
    s1=scatter(DCats{kk},LCats{kk},sc,'filled','MarkerEdgeColor',Colors(WhichColor(kk),:),'MarkerFaceColor',Colors(WhichColor(kk),:));
    s1.MarkerFaceAlpha = transp;
    s1.MarkerEdgeAlpha = transp;
end

DGrid = logspace(-4,0,1000);
LGrid = 0:.5:90;
figure(22)
PlotPSVIndexContours(aCat,LGrid,DGrid,inf,2)
PlotPSVIndexNegContours(aCat,LGrid,DGrid,inf,2)
set(gcf,'Color','w')
set(gca,'FontSize',16)
xlabel('Dipole moment')
ylabel({'Pole latitude'; '(degrees)'})
set(gca,'XScale','log')


figure(22)
hold on, plot(D(st:sp),L(st:sp),'-','Color',Colors(6,:),'LineWidth',2)
hold on, plot(D(st),L(st),'o','Color',Colors(7,:),'LineWidth',4,'MarkerSize',10)
hold on, plot(D(sp),L(sp),'x','Color',Colors(8,:),'LineWidth',4,'MarkerSize',20)
plot([1e-4 4.5],[0 0],'k')
xlim([1e-3 2.5])

f = gcf;
f.Position = [100 100 400 800];


%% select Cat-3 event
ind = Cat3Inds(30); % 27 is extreme!
st = start(ind);
sp =  stop(ind);
stLong = st-1e2;
spLong = sp+1e2;

figure
subplot(411), hold on
plot(t(stLong:spLong)-t(st),D(stLong:spLong),'k','LineWidth',1)
plot(t(st:sp)-t(st),D(st:sp),'Color',Colors(WhichColor(4),:),'LineWidth',3)
ylabel({'Dipole moment'})
xlabel('Time (ka)')
xlim([t(stLong)-t(st) t(spLong)-t(st)])
ylim([0 1.6])
set(gca,'Fontsize',16)
box off

subplot(412), hold on
plot(t(stLong:spLong)-t(st),L(stLong:spLong),'k','LineWidth',1)
plot(t(st:sp)-t(st),L(st:sp),'Color',Colors(WhichColor(4),:),'LineWidth',3)

ylabel({'Pole latitude'; '(degrees)'})
xlabel('Time (ka)')
xlim([t(stLong)-t(st) t(spLong)-t(st)])
ylim([-95 95])
yticks([-90 -45 0 45 90])
set(gca,'Fontsize',16)
box off

subplot(413), hold on
plot(t(stLong:spLong)-t(st),loga(stLong:spLong),'k','LineWidth',1)
plot(t(st:sp)-t(st),loga(st:sp),'Color',Colors(WhichColor(4),:),'LineWidth',3)
plot([t(stLong)-t(st) t(spLong)-t(st)],log10(aThStart)*[1 1],'Color',Colors(4,:))
plot([t(stLong)-t(st) t(spLong)-t(st)],log10(aThStop)*[1 1],'Color',Colors(1,:))
ylabel('log_{10}(P_{i_D})')
xlabel('Time (ka)')
xlim([t(stLong)-t(st) t(spLong)-t(st)])
ylim([-5 4])
set(gca,'Fontsize',16)
box off
set(gcf,'Color','w')

subplot(414), hold on
stairs(t(stLong:spLong)-t(st),cat(stLong:spLong)-1,'k','LineWidth',1)
stairs(t(st-1:sp)-t(st),cat(st-1:sp)-1,'Color',Colors(WhichColor(4),:),'LineWidth',3)
ylabel('Level')
xlabel('Time (ka)')
xlim([t(stLong)-t(st) t(spLong)-t(st)])
ylim([0 3])
yticks(0:3)
set(gca,'Fontsize',16)
box off
set(gcf,'Color','w')

f = gcf;
f.Position = [100 100 400 800];

sc = 20;
transp = 1;
figure(23), hold on
for kk=1:length(aCat)+1
    s1=scatter(DCats{kk},LCats{kk},sc,'filled','MarkerEdgeColor',Colors(WhichColor(kk),:),'MarkerFaceColor',Colors(WhichColor(kk),:));
    s1.MarkerFaceAlpha = transp;
    s1.MarkerEdgeAlpha = transp;
end

DGrid = logspace(-4,0,1000);
LGrid = 0:.5:90;
figure(23)
PlotPSVIndexContours(aCat,LGrid,DGrid,inf,2)
PlotPSVIndexNegContours(aCat,LGrid,DGrid,inf,2)
set(gcf,'Color','w')
set(gca,'FontSize',16)
xlabel('Dipole moment')
ylabel({'Pole latitude'; '(degrees)'})
set(gca,'XScale','log')


figure(23)
hold on, plot(D(st:sp),L(st:sp),'-','Color',Colors(6,:),'LineWidth',2)
hold on, plot(D(st),L(st),'o','Color',Colors(7,:),'LineWidth',4,'MarkerSize',10)
hold on, plot(D(sp),L(sp),'x','Color',Colors(8,:),'LineWidth',4,'MarkerSize',20)
plot([1e-4 4.5],[0 0],'k')
xlim([1e-3 2.5])

f = gcf;
f.Position = [100 100 400 800];


