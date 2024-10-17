function [fhand,p5,p95] = PlotDensity(x,bins,X)
Colors = brewermap(8,'Dark2');
n = length(x);
xAxis = (bins(2:end)+bins(1:end-1))/2;
nBins = length(bins)-1;

Density = zeros(n,nBins);
p5 = zeros(n,1);
p95 = zeros(n,1);
m = zeros(n,1);
med = zeros(n,1);

for kk=1:n
   figure(99)
   tmp = X(kk,:);
   inds = find(~isnan(tmp));
   a = histogram(tmp,bins,'Normalization','pdf');
   Density(kk,:) = a.Values;
   p5(kk) = prctile(tmp,5);
   p95(kk) = prctile(tmp,95);
   m(kk) = mean(tmp(inds));
   med(kk) = median(tmp(inds));
   close 99
end

figure
fhand = pcolor(x,xAxis,log10(Density)');
% mycolormap = parula(256); mycolormap(1,:) = 0; colormap(mycolormap)
mycolormap =  inferno(100);
mycolormap(1,:) = 0; colormap(mycolormap)
set(fhand,'EdgeColor','none')
% colorbar
hold on


stairs(x,p5,'-w','linewidth',2)
stairs(x,p95,'-w','linewidth',2)
stairs(x,m  ,'Color',Colors(1,:),'linewidth',3)
stairs(x,med,'Color',Colors(3,:),'linewidth',3)
m