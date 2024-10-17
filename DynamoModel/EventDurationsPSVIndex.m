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

[start,stop]=FindEventsActivity(a,aThStart,aThStop);
Duration = t(stop)-t(start);

figure
histogram(Duration,'Normalization','pdf','FaceColor','k','FaceAlpha',0.2)
set(gcf,'Color','w')
set(gca,'FontSize',16)
box off
xlabel('Event duration (kyr)')
ylabel('pdf')

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
    LCats{kk} = absL(inds);
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
CatRevs = zeros(NumCats-1,1);
NumCatEvents = zeros(NumCats-1,1);
RevVsEx = cell(NumCats-1,1);
for kk=1:NumEvents
    stt = start(kk);
    stp = stop(kk);
    WhatCat(kk) = max(cat(stt:stp));
    NumCatEvents(WhatCat(kk)-1) = NumCatEvents(WhatCat(kk)-1)+1;
    CatDurations{WhatCat(kk)-1} = [CatDurations{WhatCat(kk)-1}, t(stp)-t(stt)];
    if g10(stt)*g10(stp)<0 % reversal 
        CatRevs(WhatCat(kk)-1)=CatRevs(WhatCat(kk)-1)+1;
        RevVsEx{WhatCat(kk)-1} = [RevVsEx{WhatCat(kk)-1} 1];
    else
        RevVsEx{WhatCat(kk)-1} = [RevVsEx{WhatCat(kk)-1} 0];
    end
end

fprintf('Number of events: %g\n',sum(NumCatEvents))
fprintf('Cat 1: %g (%g percent) \n',NumCatEvents(1),NumCatEvents(1)/sum(NumCatEvents)*100)
fprintf('Cat 2: %g (%g percent) \n',NumCatEvents(2),NumCatEvents(2)/sum(NumCatEvents)*100)
fprintf('Cat 3: %g (%g percent) \n',NumCatEvents(3),NumCatEvents(3)/sum(NumCatEvents)*100)

disp(' ')
disp(' ')

fprintf('Number of reversals per category: %g \n',sum(CatRevs))
fprintf('Cat 1: %g (%g percent) \n',CatRevs(1),CatRevs(1)/NumCatEvents(1)*100)
fprintf('Cat 2: %g (%g percent) \n',CatRevs(2),CatRevs(2)/NumCatEvents(2)*100)
fprintf('Cat 3: %g (%g percent) \n',CatRevs(3),CatRevs(3)/NumCatEvents(3)*100)




figure
b=bar(1:NumCats-1,NumCatEvents/sum(NumCatEvents)*100,'k','FaceAlpha',0.3,'LineWidth',2);
set(gcf,'Color','w')
set(gca,'FontSize',16)
box off
xlabel('Category')
ylabel('Count (per cent)')


figure(102)
hold on,h1 = histogram(CatDurations{1},15,'Normalization','pdf');
hold on,h2 = histogram(CatDurations{2},10,'Normalization','pdf');
hold on,h3 = histogram(CatDurations{3},10,'Normalization','pdf');
xlabel('Event duration (kyr)'),ylabel('pdf')

figure,hold on
stairs(h1.BinEdges(2:end),h1.Values,'Color',Colors(WhichColor(2),:),'LineWidth',3);
stairs(h2.BinEdges(2:end),h2.Values,'Color',Colors(WhichColor(3),:),'LineWidth',3);
stairs(h3.BinEdges(2:end),h3.Values,'Color',Colors(WhichColor(4),:),'LineWidth',3);
legend('Cat-1','Cat-2','Cat-3','Location','NorthEast')
set(gcf,'Color','w')
set(gca,'FontSize',16)
xlabel('Event duration (kyr)'),ylabel('pdf')
close(102)


%% Kernel density estimates
[f1,xi1] = ksdensity(CatDurations{1},'Support','positive');
[f2,xi2] = ksdensity(CatDurations{2},'Support','positive');
[f3,xi3] = ksdensity(CatDurations{3},'Support','positive');
% no categories
[f,xi]=ksdensity(Duration,'Support','positive');

figure(103),hold on
% plot(xi,f,'Color',[.7 .7 .7],'LineWidth',3);
plot(xi1,f1,'Color',Colors(WhichColor(2),:),'LineWidth',3);
plot(xi2,f2,'Color',Colors(WhichColor(3),:),'LineWidth',3);
plot(xi3,f3,'Color',Colors(WhichColor(4),:),'LineWidth',3);
legend('Cat-1','Cat-2','Cat-3','Location','NorthEast')
set(gcf,'Color','w')
set(gca,'FontSize',16)
xlabel('Event duration (kyr)'),ylabel('pdf')
xlim([0 50])

%% Mean, median and std of durations
% disp(' '),disp(' ')
% fprintf('Cat-1: Mean event duration:%g kr\n',mean(CatDurations{1}))
% fprintf('Cat-1: Median event duration:%g kr\n',median(CatDurations{1}))
% fprintf('Cat-1: Standard deviation of event duration:%g kr\n',std(CatDurations{1}))
% 
% disp(' '),disp(' ')
% fprintf('Cat-2: Mean event duration:%g kr\n',mean(CatDurations{2}))
% fprintf('Cat-2: Median event duration:%g kr\n',median(CatDurations{2}))
% fprintf('Cat-2: Standard deviation of event duration:%g kr\n',std(CatDurations{2}))
% 
% disp(' '),disp(' ')
% fprintf('Cat-3: Mean event duration:%g kr\n',mean(CatDurations{3}))
% fprintf('Cat-3: Median event duration:%g kr\n',median(CatDurations{3}))
% fprintf('Cat-3: Standard deviation of event duration:%g kr\n',std(CatDurations{3}))


%% Separate reversals from no-reversals
EvDursCat2 = CatDurations{2};
EvDursCat2e = EvDursCat2(RevVsEx{2}==0);
EvDursCat2r = EvDursCat2(RevVsEx{2}==1);
EvDursCat3 = CatDurations{3};
EvDursCat3e = EvDursCat3(RevVsEx{3}==0);
EvDursCat3r = EvDursCat3(RevVsEx{3}==1);


disp(' '),disp(' ')
fprintf('Event duration, Cat-1: Mean, median, sig: %g, %g, %g kr\n',...
    mean(CatDurations{1}),median(CatDurations{1}),std(CatDurations{1}))

disp(' ')
fprintf('Event duration, Cat-2: Mean, median, sig: %g, %g, %g kr\n',...
    mean(EvDursCat2),median(EvDursCat2),std(EvDursCat2))
fprintf('Event duration, Cat-2, no revs.: Mean, median, sig: %g, %g, %g kr\n',...
    mean(EvDursCat2e),median(EvDursCat2e),std(EvDursCat2e))
fprintf('Event duration, Cat-2, revs.: Mean, median, sig: %g, %g, %g kr\n',...
    mean(EvDursCat2r),median(EvDursCat2r),std(EvDursCat2r))

disp(' ')
fprintf('Event duration, Cat-3: Mean, median, sig: %g, %g, %g kr\n',...
    mean(EvDursCat3),median(EvDursCat3),std(EvDursCat3))
fprintf('Event duration, Cat-3, no revs.: Mean, median, sig: %g, %g, %g kr\n',...
    mean(EvDursCat3e),median(EvDursCat3e),std(EvDursCat3e))
fprintf('Event duration, Cat-3, revs.: Mean, median, sig: %g, %g, %g kr\n',...
    mean(EvDursCat3r),median(EvDursCat3r),std(EvDursCat3r))

%% KDEs
% Cat-2
[f2,xi2] = ksdensity(EvDursCat2e,'Support','positive');
[f3,xi3] = ksdensity(EvDursCat2r,'Support','positive');

figure(104),hold on
plot(xi2,f2,'--','Color',Colors(WhichColor(3),:),'LineWidth',3);
plot(xi3,f3,':','Color',Colors(WhichColor(3),:),'LineWidth',3);
set(gcf,'Color','w')
set(gca,'FontSize',16)
xlabel('Event duration (kyr)'),ylabel('pdf')
xlim([0 50])

% Cat-3
[f2,xi2] = ksdensity(EvDursCat3e,'Support','positive');
[f3,xi3] = ksdensity(EvDursCat3r,'Support','positive');

plot(xi2,f2,'--','Color',Colors(WhichColor(4),:),'LineWidth',3);
plot(xi3,f3,':','Color',Colors(WhichColor(4),:),'LineWidth',3);
legend('Cat-2, No rev.','Cat-2, Rev.','Cat-3, No rev.','Cat-3, Rev.','Location','NorthEast')
set(gcf,'Color','w')
set(gca,'FontSize',16)
xlabel('Event duration (kyr)'),ylabel('pdf')
xlim([0 50])


