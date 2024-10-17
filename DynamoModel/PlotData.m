clear
close all
clc
Colors = brewermap(8,'Dark2');

% load('./DataFiles/ghl1.mat');
% t = time;
% g10 = ghl1(:,1);
% g11 = ghl1(:,2);
% h11 = ghl1(:,3);
% D = sqrt(g10.^2 + g11.^2 + h11.^2);
% D = D/mean(D);
% Ds = D.*sign(g10);
% 
% coLat = -atan2(sqrt(g11.^2 +h11.^2),g10);
% Lrad = pi/2-coLat-pi;
% L = Lrad*180/pi;
% 
% a =  ((pi/2)-abs(Lrad))./(pi*D);
% loga = log10(a);

load('./DataFiles/myData.mat');
ChunkLength = floor(length(t)/3);
t = t*1e-3;

for kk=1:3
    inds = (kk-1)*ChunkLength+1:kk*ChunkLength;
    figure
    plot(t(inds),Ds(inds),'k')
    xlim([t(inds(1)) t(inds(end))])
    ylim([-2 2])
    yticks([-1 0 1])
    set(gcf,'Color','w')
    set(gca,'FontSize',40)
    box off
    f = gcf;
    f.Position = [100 100 800 400];
end

for kk=1:3
    inds = (kk-1)*ChunkLength+1:kk*ChunkLength;
    figure
    plot(t(inds),L(inds),'k')
    xlim([t(inds(1)) t(inds(end))])
    ylim([-90 90])
    yticks([-90 -45 0 45 90])
    set(gcf,'Color','w')
    set(gca,'FontSize',40)
    box off
    f = gcf;
    f.Position = [100 100 800 400];
end

for kk=1:3
    inds = (kk-1)*ChunkLength+1:kk*ChunkLength;
    figure
    plot(t(inds),loga(inds),'k')
    xlim([t(inds(1)) t(inds(end))])
    ylim([-4 2])
    yticks([-3 -1 1])
    set(gcf,'Color','w')
    set(gca,'FontSize',40)
    box off
    f = gcf;
    f.Position = [100 100 800 400];
end


ChunkLength = floor(length(t)/10);
kk=2;
inds = (kk-1)*ChunkLength+1:kk*ChunkLength;
figure
plot(t(inds),Ds(inds),'k')
xlim([t(inds(1)) t(inds(end))])
ylim([-2 2])
yticks([-1 0 1])
set(gcf,'Color','w')
set(gca,'FontSize',40)
box off
f = gcf;
f.Position = [100 100 1200 400];
ylabel({'Dipole moment','(signed)'})


figure
plot(t(inds),L(inds),'k')
xlim([t(inds(1)) t(inds(end))])
ylim([-90 90])
yticks([-90 -45 0 45 90])
set(gcf,'Color','w')
set(gca,'FontSize',40)
box off
f = gcf;
f.Position = [100 100 1200 400];
ylabel('Pole latitude')

figure
plot(t(inds),loga(inds),'k')
xlim([t(inds(1)) t(inds(end))])
ylim([-4 2])
yticks([-3 -1 1])
set(gcf,'Color','w')
set(gca,'FontSize',40)
box off
f = gcf;
f.Position = [100 100 1200 400];
ylabel({'PSV index','(log)'})
xlabel('Time (Myr)')
