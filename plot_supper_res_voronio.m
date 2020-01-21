clear,clc
fold = 'C:\Users\Zoe\Google Drive\paper 2\diffusion study\100\';

load([fold num2str(loc) '.mat'])
[v,c] = voronoin([data(:,3) data(:,4)]); % lx and ly are the center of each steps
area = nan(1,length(c));
for i = 1:length(c) % this loop finds all area in each polygon
    ind = c{i}';
    if ind ~= 1 % all(c{i}~=1)
        area(i) = polyarea(v(ind,1),v(ind,2));
    end
end
arealogn = log10(area);
bins1 = linspace(-1,5,254);
[h1, m1] = histc(arealogn,bins1);
col = colormap(flipud(hot(255)));
figure(1)
% set(gcf,'position', [16.00        477.00       1863.00        489.00])
% subplot(1,3,1)
for i = 1:length(c)
    ind = c{i}';
    patch(v(ind,1),v(ind,2),col(m1(i)+1,:),'EdgeColor','none'); % use color i.
end
xlim([min(data(:,3))-150 max(data(:,3))-150])
ylim([min(data(:,4))-150 max(data(:,4))-150])
axis tight
box off
set(gca,'XTickLabel',[],'YTickLabel',[]);
