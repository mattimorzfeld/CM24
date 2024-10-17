clear
close all
clc
Colors = brewermap(8,'Dark2');


load('./DataFiles/myData.mat');
aThStart = 0.5;
aThStop = 0.035;
aCat = [.5 3 15];;

dt = t(2)-t(1);
absL = abs(L);

[start,stop]=FindEventsActivity(a,aThStart,aThStop);
Duration = t(stop)-t(start);
IndDuration = stop-start;
maxIndDuration = max(IndDuration);

NumCats = length(aCat)+1;
IndsCats = cell(NumCats,1);
AllInds = 1:length(D);

for kk=1:NumCats
    if kk==1
        inds = a<aCat(1);
    elseif kk==NumCats
        inds = a>aCat(end);
    else
        inds = ((a<aCat(kk)) & (a>aCat(kk-1)));
    end
    IndsCats{kk} = AllInds(inds);
end
%% Label time series by categories
cat = zeros(length(D),1);
for kk=1:NumCats
    cat(IndsCats{kk}) = kk;
end
%% Find Category of each event
WhatCat = zeros(length(start),1);
for kk=1:length(start)
    stt = start(kk);
    stp = stop(kk);
    WhatCat(kk) = max(cat(stt:stp));
end

HowMuchBefore = 250;
HowMuchAfter  = 500;
WhereToStart = find(start>HowMuchBefore,1);
WhereToStop = find(start+HowMuchAfter>length(D),1);
if isempty(WhereToStop)
    WhereToStop=length(start);
else
    WhereToStop = WhereToStop(1);
end
NumEvents = WhereToStop-WhereToStart;
WhatCat = WhatCat(WhereToStart:WhereToStop);

DPlot = zeros(HowMuchBefore+HowMuchAfter+1,NumEvents);
LPlot = zeros(HowMuchBefore+HowMuchAfter+1,NumEvents);
logaPlot = zeros(HowMuchBefore+HowMuchAfter+1,NumEvents);
tPlot = (1:(HowMuchBefore+HowMuchAfter+1))*dt;
WhichCat = zeros(NumEvents,1);

counter = 0;
for kk=WhereToStart:WhereToStop
    counter = counter+1;
    inds = start(kk)-HowMuchBefore:start(kk)+HowMuchAfter;
    DPlot(:,counter)=D(inds);
    LPlot(:,counter)=absL(inds);
    logaPlot(:,counter)=loga(inds);
    %% exclude previous events
    if start(kk)-HowMuchBefore<start(kk-1) 
        EmptyInds = 1:find(inds==stop(kk-1));
        DPlot(EmptyInds,counter)=nan;
        LPlot(EmptyInds,counter)=nan;
        logaPlot(EmptyInds,counter)=nan;
    end
    if start(kk)-HowMuchBefore<stop(kk-1) 
        EmptyInds = 1:find(inds==stop(kk-1));
        DPlot(EmptyInds,counter)=nan;
        LPlot(EmptyInds,counter)=nan;
        logaPlot(EmptyInds,counter)=nan;  
    end
    %% exclude next events
    if kk<WhereToStop
        if start(kk)+HowMuchAfter>start(kk+1)
            EmptyInds = find(inds==start(kk+1));
            DPlot(end:-1:EmptyInds,counter)=nan;
            LPlot(end:-1:EmptyInds,counter)=nan;
            logaPlot(end:-1:EmptyInds,counter)=nan;
        end
    end

end


for kk=min(WhatCat):max(WhatCat)
    inds = find(WhatCat==kk);

    %% histograms
    bins = linspace(0,1.8,50);
    PlotDensity(tPlot,bins,DPlot(:,inds))
    xlabel('Time (kyr)')
%     ylabel('Normalized |g_1^0|')
    set(gca,'FontSize',32)
    box off
    set(gcf,'Color','w')
    hold on, plot(HowMuchBefore*dt*[1 1],[0,max(D)],'w-','LineWidth',2)
    hold on, plot(tPlot, ones(size(tPlot)),'k--','LineWidth',2)
    drawnow
    FileName = strcat('./PNGs/Dip_Cat',num2str(kk-1),'.png');
    saveas(gcf,FileName)
    
    bins = linspace(0,95,50);
    PlotDensity(tPlot,bins,LPlot(:,inds))
    xlabel('Time (kyr)')
%     ylabel('Pole latitude')
    set(gca,'FontSize',32)
    box off
    set(gcf,'Color','w')
    hold on, plot(HowMuchBefore*dt*[1 1],[0,90],'w-','LineWidth',2)
    hold on, plot(tPlot, mean(absL)*ones(size(tPlot)),'k--','LineWidth',2)
    drawnow
    FileName = strcat('./PNGs/Lat_Cat',num2str(kk-1),'.png');
    saveas(gcf,FileName)
    
    bins = logspace(-2,3,50);
    bins = log10(bins);
    % bins = linspace(min(loga),max(loga),100);
    PlotDensity(tPlot,bins,logaPlot(:,inds))
    xlabel('Time (kyr)')
%     ylabel('Log_{10} of PSV index')
    set(gca,'FontSize',32)
    box off
    set(gcf,'Color','w')
    hold on,
    plot(HowMuchBefore*dt*[1 1],[-2,3],'w-','LineWidth',2)
    plot(tPlot,log10(0.5)*ones(size(tPlot)),'k--','LineWidth',2)
    plot(tPlot,log10(0.05)*ones(size(tPlot)),'k--','LineWidth',2)
    drawnow
    FileName = strcat('./PNGs/loga_Cat',num2str(kk-1),'.png');
    saveas(gcf,FileName)
end