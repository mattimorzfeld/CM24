clear
close all
clc

Colors = brewermap(9,'Set2');
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

%% Time to reach level 2
TtRMaxLev2 = zeros(length(Cat2Inds),1);
TtRMaxLev2e = zeros(length(Cat2Inds),1);
TtRMaxLev2r = zeros(length(Cat2Inds),1);
TtDecfLev2 = zeros(length(Cat2Inds),1);
TtDecfLev2e = zeros(length(Cat2Inds),1);
TtDecfLev2r = zeros(length(Cat2Inds),1);
for oo=1:length(Cat2Inds)
    ind = Cat2Inds(oo);
    st_tmp = start(ind);
    sp_tmp  = stop(ind);
    levs = cat(st_tmp:sp_tmp)-1;
    ReachLev2Ind = find(levs==2,1,'first');
    DecayLev2Ind = find(levs==2,1,'last');
    TtRMaxLev2(oo) = t(st_tmp+ReachLev2Ind)-t(st_tmp);
    TtDecfLev2(oo) = t(sp_tmp)-t(st_tmp+DecayLev2Ind);
    if g10(st_tmp)*g10(sp_tmp)<0 % reversal
        TtRMaxLev2r(oo) = t(st_tmp+ReachLev2Ind)-t(st_tmp);
        TtDecfLev2r(oo) = t(sp_tmp)-t(st_tmp+DecayLev2Ind);
    else 
        TtRMaxLev2e(oo) = t(st_tmp+ReachLev2Ind)-t(st_tmp);
        TtDecfLev2e(oo) = t(sp_tmp)-t(st_tmp+DecayLev2Ind);
    end
end
TtRMaxLev2e = TtRMaxLev2e(TtRMaxLev2e~=0);
TtDecfLev2e = TtDecfLev2e(TtDecfLev2e~=0);
TtRMaxLev2r = TtRMaxLev2r(TtRMaxLev2r~=0);
TtDecfLev2r = TtDecfLev2r(TtDecfLev2r~=0);

disp(' ')
fprintf('Avg. time (2 std) to reach level 2: %g (%g)\n',mean(TtRMaxLev2),std(TtRMaxLev2))
fprintf('Avg. time (2 std) to reach level 2 (no rev): %g (%g)\n',mean(TtRMaxLev2e),std(TtRMaxLev2e))
fprintf('Avg. time (2 std) to reach level 2 (rev): %g (%g)\n',mean(TtRMaxLev2r),std(TtRMaxLev2r))
disp(' ')
fprintf('Avg. time (2 std) to decay from level 2: %g (%g)\n',mean(TtDecfLev2),std(TtDecfLev2))
fprintf('Avg. time (2 std) to decay from level 2 (no rev): %g (%g)\n',mean(TtDecfLev2e),std(TtDecfLev2e))
fprintf('Avg. time (2 std) to decay from level 2 (rev): %g (%g)\n',mean(TtDecfLev2r),std(TtDecfLev2r))


% histograms
figure, 
subplot(221),hold on
histogram(TtRMaxLev2,12,'Normalization','pdf')
histogram(TtRMaxLev2e,12,'Normalization','pdf')
histogram(TtRMaxLev2r,12,'Normalization','pdf')
legend('All','No Rev.','Rev.')
xlabel('Time (kyr)')
ylabel('pdf')
set(gcf,'Color','w')
set(gca,'FontSize',16)
title('Time to reach lev. 2')

subplot(223),hold on
histogram(TtDecfLev2,15,'Normalization','pdf')
histogram(TtDecfLev2e,15,'Normalization','pdf')
histogram(TtDecfLev2r,15,'Normalization','pdf')
legend('All','No Rev.','Rev.')
xlabel('Time (kyr)')
ylabel('pdf')
set(gcf,'Color','w')
set(gca,'FontSize',16)
title('Time to decay from lev. 2')


%% Time to reach level 3
TtRMaxLev3 = zeros(length(Cat3Inds),1);
TtRMaxLev3e = zeros(length(Cat3Inds),1);
TtRMaxLev3r = zeros(length(Cat3Inds),1);
TtDecfLev3 = zeros(length(Cat3Inds),1);
TtDecfLev3e = zeros(length(Cat3Inds),1);
TtDecfLev3r = zeros(length(Cat3Inds),1);
for oo=1:length(Cat3Inds)
    ind = Cat3Inds(oo);
    st_tmp = start(ind);
    sp_tmp  = stop(ind);
    levs = cat(st_tmp:sp_tmp)-1;
    ReachLev3Ind = find(levs==3,1,'first');
    DecayLev3Ind = find(levs==3,1,'last');
    TtRMaxLev3(oo) = t(st_tmp+ReachLev3Ind)-t(st_tmp);
    TtDecfLev3(oo) = t(sp_tmp)-t(st_tmp+DecayLev3Ind);
    if g10(st_tmp)*g10(sp_tmp)<0 % reversal
        TtRMaxLev3r(oo) = t(st_tmp+ReachLev3Ind)-t(st_tmp);
        TtDecfLev3r(oo) = t(sp_tmp)-t(st_tmp+DecayLev3Ind);
    else 
        TtRMaxLev3e(oo) = t(st_tmp+ReachLev3Ind)-t(st_tmp);
        TtDecfLev3e(oo) = t(sp_tmp)-t(st_tmp+DecayLev3Ind);
    end
end
TtRMaxLev3e = TtRMaxLev3e(TtRMaxLev3e~=0);
TtDecfLev3e = TtDecfLev3e(TtDecfLev3e~=0);

TtRMaxLev3r = TtRMaxLev3r(TtRMaxLev3r~=0);
TtDecfLev3r = TtDecfLev3r(TtDecfLev3r~=0);
disp(' ')
fprintf('Avg. time (2 std) to reach level 3: %g (%g)\n',mean(TtRMaxLev3),std(TtRMaxLev3))
fprintf('Avg. time (2 std) to reach level 3 (no rev): %g (%g)\n',mean(TtRMaxLev3e),std(TtRMaxLev3e))
fprintf('Avg. time (2 std) to reach level 3 (rev): %g (%g)\n',mean(TtRMaxLev3r),std(TtRMaxLev3r))
disp(' ')
fprintf('Avg. time (2 std) to decay from level 3: %g (%g)\n',mean(TtDecfLev3),std(TtDecfLev3))
fprintf('Avg. time (2 std) to decay from level 3 (no rev): %g (%g)\n',mean(TtDecfLev3e),std(TtDecfLev3e))
fprintf('Avg. time (2 std) to decay from level 3 (rev): %g (%g)\n',mean(TtDecfLev3r),std(TtDecfLev3r))


% histograms 
subplot(222),hold on
histogram(TtRMaxLev3,12,'Normalization','pdf')
histogram(TtRMaxLev3e,12,'Normalization','pdf')
histogram(TtRMaxLev3r,12,'Normalization','pdf')
legend('All','No Rev.','Rev.')
xlabel('Time (kyr)')
ylabel('pdf')
set(gcf,'Color','w')
set(gca,'FontSize',16)
title('Time to reach lev. 3')

subplot(224),hold on
histogram(TtDecfLev3,15,'Normalization','pdf')
histogram(TtDecfLev3e,15,'Normalization','pdf')
histogram(TtDecfLev3r,15,'Normalization','pdf')
legend('All','No Rev.','Rev.')
xlabel('Time (kyr)')
ylabel('pdf')
set(gcf,'Color','w')
set(gca,'FontSize',16)
title('Time to decay from lev. 3')
