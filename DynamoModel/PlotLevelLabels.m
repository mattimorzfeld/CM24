clear
close all
clc

load('./DataFiles/myData.mat');

absL = abs(L);

Colors = brewermap(8,'Set1');
WhichColor = [3 2 5 1];

aCat = [.5 3 15];

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
TinCat

figure(1), hold on
for kk=1:length(aCat)+1
    scatter(DCats{kk},LCats{kk},120,'filled','MarkerFaceColor',Colors(WhichColor(kk),:),'MarkerFaceAlpha',1)
end

DGrid = logspace(-4,0,1000);
LGrid = 0:.5:90;
figure(1)
PlotPSVIndexContours(aCat,LGrid,DGrid,1,3)

set(gcf,'Color','w')
set(gca,'FontSize',16)
xlabel('Dipole moment')
ylabel('Pole latitude (degree)')
set(gca,'XScale','log')
xlim([1e-3 3])
% FileName = strcat('./PNGs/LabeledData_Transparent.png');
% saveas(gcf,FileName)

%% Label by categories
cat = zeros(length(D),1);
for kk=1:NumCats
    cat(IndsCats{kk}) = kk;
end

figure
stairs(t*1e-3,cat-1,'k','LineWidth',1)
set(gcf,'Color','w')
set(gca,'FontSize',16)
xlabel('Time (Myr)')
ylabel('Level')
ylim([0 NumCats])
box off

%% Transition matrix: ~40 years
dt = t(2)-t(1); % about 43 years
T = zeros(NumCats);
for kk=1:length(D)-1
    i = cat(kk);
    j = cat(kk+1);
    T(i,j) = T(i,j)+1;
end
% normalize
st = sum(T');
Tn = diag(1./st)*T*100

figure
imagesc(Tn)
colormap([1 1 1; brewermap([],'Blues')])
set(gcf,'Color','w')
set(gca,'FontSize',16)
xlabel('Level')
ylabel('Level')
yticklabels({'normal','1','2','3','4'})
yticks([1:5])
xticklabels({'normal','1','2','3','4'})
xticks([1:5])
colorbar

%% Figure 
% TnPlot = zeros(4);
TnPlot =-triu(Tn)+tril(Tn);
imagesc(TnPlot)
colormap(brewermap(100,'RdBu'))
set(gcf,'Color','w')
set(gca,'FontSize',16)
xlabel('Level')
ylabel('Level')
yticklabels({'normal','1','2','3'})
yticks([1:4])
xticklabels({'normal','1','2','3'})
xticks([1:4])
clim([-2 2])

% %% Transition matrix: 400 years
% dt = t(2)-t(1); % about 43 years
% Dt =  1000;     % about DT*dt years
% T = zeros(NumCats);
% for kk=1:length(D)-Dt
%     tmpCats = cat(kk+1:kk+Dt);
%     CatsEncountered = zeros(NumCats,1);
%     for jj=1:NumCats
%         tmp = (tmpCats==jj);
%         CatsEncountered(jj) = sum(tmp);
%     end
%     i = cat(kk);
%     for jj=1:NumCats
%         T(i,jj) = T(i,jj)+CatsEncountered(jj);
%     end
% end
% % normalize
% st = sum(T');
% Tn = diag(1./st)*T*100
% 
% figure
% imagesc(Tn)
% colormap([1 1 1; brewermap([],'Blues')])
% set(gcf,'Color','w')
% set(gca,'FontSize',16)
% xlabel('Category')
% ylabel('Category')
% yticks([0:5])
% xticks([0:5])
% colorbar