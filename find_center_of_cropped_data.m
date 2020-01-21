clear,clc,close all
text1 = '2 '; % drift corrected number
text3 = '1'; % NP number
text2 = '1000'; %cropped size
text = [text1 'cropped ' text2 ' nm' text3];
fold = 'C:\Zoe\phase study\20191218\E\';
load([fold text '.mat']);
data = cropped;
figure(5),hist(data(:,2),50)
x = data(:,3);
y = data(:,4);
figure(1),plot(x,y,'.');axis equal

mX = [x y];
bins = 50;
vXEdge = linspace(min(x),max(x),bins);
vYEdge = linspace(min(y),max(y),bins);
z = hist2d (mX, vXEdge, vYEdge);
z = transpose(z);
[X Y] = meshgrid(vXEdge,vYEdge);
X = X(1:end-1,1:end-1);
Y = Y(1:end-1,1:end-1);
figure(2), surf(X,Y,z)
shading flat
view([0 90])
axis equal tight
% saveas(gcf,[fold text3 ' super res.jpg'],'jpg')
%% 2D gaussian fit - good for fit a bell curve localization
param = [50 7 7 25 25 0]; % 4 and 5 th element are x0 and y0
fresult = fit2DGauss_funct(z,param)
figure(3), plot(x,y,'x')
hold on; plot3(min(x)+fresult.x0*(max(x)-min(x))/(bins),min(y)+fresult.y0*(max(y)-min(y))/(bins),25,'ko','markerfacecolor','k');
hold off;
axis equal tight
center = [min(x)+fresult.x0*(max(x)-min(x))/(bins) min(y)+fresult.y0*(max(y)-min(y))/(bins)];
%%  center cropped data
load([fold text1 'cropped 2000 nm' text3 '.mat']);
a = cropped;
% a = cropped(cropped(:,2)>1000,:);
a(:,3) = a(:,3) - center(1); % center is in nm but a is in micron size.
a(:,4) = a(:,4) - center(2);
save([fold text1 text3 'cropped centered 2000nm.mat'],'a','center')