clear,clc, close all
fold = 'C:\Users\Zoe\Google Drive\paper 3\diffusion study\100\36C\';
name = [fold(51:53) ' ' fold(55:end-1)];
% name = [fold(20:22) ' ' fold(24:26)];
MyFileInfo = dir([fold 'centered\*cropped centered 2000nm.mat']);
nfiles = size(MyFileInfo,1);
col = colormap(hsv(nfiles));
SD = nan(nfiles,30000);
Daway = nan(nfiles,30000);
Sigma = nan(nfiles,1);
index = 0:0.05:1;
deltat = 0.00186;
exposuret = 0.0017;
Dall = nan(nfiles,length(index)-1);
Dallerr = nan(nfiles,length(index)-1,1);
%
for j = 1:nfiles
    text = MyFileInfo(j).name(1:8);
    load([fold 'centered\' MyFileInfo(j).name])
    Sigma(j) = mean(a(:,10))*10^-3;
    load([fold 'agg study\ratio used Shelby.mat'])
    agg = agg_used(j,1);
    load([fold 'agg study\' text ' ' num2str(agg*100) ' 400.mat'])
    a = trackArtData;
    get_rid_of_all_long_traj 
    % code above get rid of trajectories longer than 32 steps
    X = a(:,5); % x in micron
    Y = a(:,6); % y
    traj_number = a(:,1); % cnt for traj
    for i = 2:size(a,1)
        if traj_number(i) == traj_number(i-1)
            meanX = (X(i)+X(i-1))/2;
            meanY = (Y(i)+Y(i-1))/2;
            SD(j,i) = (X(i)-X(i-1))^2+(Y(i)-Y(i-1))^2;
            Daway(j,i) = sqrt(meanX^2+meanY^2);
        else
            SD(j,i) = nan;
            Daway(j,i) = nan;
        end
    end
end

for j = 1:nfiles
    daway = Daway(j,:);
    Singlestep = sqrt(SD(j,:));
    daway(isnan(daway)) = [];
    Singlestep(isnan(Singlestep)) = [];
    [m n] = histc(daway, index);
    sigma = Sigma(j);
    for i  = 1:length(index)-1
        st = Singlestep(n == i);
        st(st == 0) = [];
        if length(st) > 5
            bins = linspace(0,max(st),30);
            [hd md] = histc(st,bins);
            hd = hd/sum(hd);
            x = bins+bins(2)/2;    
            ftype = fittype(['a*x.*exp(-x.^2/(4*' num2str(deltat) '*d))/(2*' num2str(deltat) '*d)']);
            fres = fit(x',(hd/bins(2))',ftype,'start',[1 1],'lower',[0.8 0.155],'upper',[1.2 5])
            fittingpara = confint(fres);
            d = (fittingpara(1,2)+fittingpara(2,2))/2;
            Dallerr(j,i) = abs(fittingpara(1,2)-fittingpara(2,2))/2;
            Dall(j,i) = (d-sigma^2/(2*deltat))/(1-exposuret/(3*deltat));
        else
            Dallerr(j,i) = nan;
            Dall(j,i) = nan;
        end
    end
    figure(1)
    plot(index(1:end-1)+(index(2)-index(1))/2,Dall(j,:),'o-.','color',col(j,:),'LineWidth',3)
    hold on
end

for i  = 1:length(index)-1    
    for j = 1: nfiles
        if Dallerr(j,i) > 1
            Dall(j,i) = nan;
            Dallerr(j,i) = nan;
        end
    end
end
Dtot = nanmean(Dall);
Dtoterr = nanstd(Dall);
plot(index(1:end-1)+(index(2)-index(1))/2,Dtot,'sk','MarkerSize',20,'LineWidth',3)
xlabel('Distance away from NP center (nm)','fontsize',25)
ylabel('D_x_y (\mum^2/s)','fontsize',25)
set(gca,'fontsize',25,'linewidth',4)
set(gcf,'position', [12   45.5  1138.5  572.5])
saveas(gcf,[fold name 'D with correction 50 1000.jpg'],'jpg')
save([fold 'Dmb average NP.mat'])