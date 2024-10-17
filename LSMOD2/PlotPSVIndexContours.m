function PlotPSVIndexContours(levels,LGrid,DGrid,CN,LW)


[mD,mL] = meshgrid(DGrid,LGrid);
ma = zeros(size(mD));
for kk=1:size(mD,1)
    for jj=1:size(mD,2)
        Lrad = LGrid(kk) * (pi/180);
        a =  ((pi/2)-abs(Lrad))./(pi*DGrid(jj))*80;

        ma(kk,jj) = a;
    end
end
hold on


if isfinite(CN)
    [M,c] = contour(mD,mL,ma,levels,'k','ShowText','off',...
        'LineWidth',LW,'LabelSpacing',100);
else
    [M,c] = contour(mD,mL,ma,levels,'Color',[.7 .7 .7],'ShowText','off',...
        'LineWidth',LW,'LabelSpacing',100);    
end