% this is for images from andor (128x128) images are used for diffusion
% study and already know the used NP for linking (NPused)
% Run one cell at a time!
clear,clc

text1 = '0809';
text2 = 'B';
npused = [5];
fold = ['C:\Users\Zoe\Google Drive\paper 3\diffraction image intensity study\POPC TF\' text1 '\' text2 '\'];
name = [text1 ' ' text2];
file2 = [fold 'POPC TF_01.tif'];
fileNP = 'NP405.csv';
NP405_tot = csvread([fold fileNP],1,0);
NP405 = NP405_tot(npused,:);
A = imread(file2);
% A = im2double(A);
A = cast(A,'double');
res = 80;
SizeA = 0:res:10160;
clear i j X Y
figure(4)
imagesc(SizeA,SizeA,A) % 0:511 pixels, with 80nm per pixel
hold on
plot(NP405(:,3),NP405(:,4),'ro','Linewidth',3)
hold off
axis equal tight
saveas(gcf,[fold 'NPloc.jpg'],'jpg')
B = nan(15,15,size(NP405,1));

% crop out each NP405 loc on membrane and find I_c and I_s
clear i j X Y
X_NP = NP405(:,3);
Y_NP = NP405(:,4);
X_NP_mem = nan(size(X_NP));
Y_NP_mem = nan(size(Y_NP));
smallboxwidth = 16; % make this an even number, final size = this size plus 1
cropped = nan(smallboxwidth+1,smallboxwidth+1);
% % % % If need to break the for loop, use the code below load existing data and
% % % % change i to where it breaked
% load([fold text ' cropped diff limited membrane2.mat'])
for i =1:size(NP405,1)
    X_NP_mem(i) = find(abs(SizeA - X_NP(i)) == min(abs(SizeA - X_NP(i))));
    Y_NP_mem(i) = find(abs(SizeA - Y_NP(i)) == min(abs(SizeA - Y_NP(i))));
    % X_NP_mem and Y_NP_mem are the element NP at membrane location in 
    % diffraction limited ima0ge maxtrix unit in px
    [m, n] = size(A);
    ymintemp = max([(Y_NP_mem(i)-smallboxwidth/2) 1]);
    xmintemp = max([(X_NP_mem(i)-smallboxwidth/2) 1]);
    ymaxtemp = min([(Y_NP_mem(i)+smallboxwidth/2) m]);
    xmaxtemp = min([(X_NP_mem(i)+smallboxwidth/2) n]);
    
    xrange = length(xmintemp:xmaxtemp);
    yrange = length(ymintemp:ymaxtemp);
    cropped(1:yrange,1:xrange) = A(ymintemp:ymaxtemp,xmintemp:xmaxtemp);
    z = cropped(:,:) - min(min(cropped(:,:)));
    z = z./max(max(z))*10;
    z(isnan(z) == 1) = 0;
%     param = [5 5 5 7 7 0];
    param = [5 2 2 (smallboxwidth)/2+1 (smallboxwidth)/2+1 0]; % 4 and 5 th element are x0 and y0
    fresult = fit2DGauss_diff_funct(z,param)
%     center = [ceil(fresult.x0+X_NP_mem(i)-(smallboxwidth/2+1)) ceil(fresult.y0+Y_NP_mem(i)-(smallboxwidth/2+1))]
    center = [round(fresult.x0+X_NP_mem(i)-(smallboxwidth/2+1)) round(fresult.y0+Y_NP_mem(i)-(smallboxwidth/2+1))]
    figure(5),
    imagesc(X_NP_mem(i)-smallboxwidth/2:X_NP_mem(i)+smallboxwidth/2,...
                Y_NP_mem(i)-smallboxwidth/2:Y_NP_mem(i)+smallboxwidth/2,cropped);
    hold on;
%     text(fresult.x0+X_NP_mem(i)-(smallboxwidth/2+1),fresult.y0+Y_NP_mem(i)-(smallboxwidth/2+1),'center')
    plot(center(1),center(2),'ko','LineWidth',4)
    hold off;
    axis tight

    aaaa = input('Is this a good location? (1 = yes, 0 = no): ');
    if aaaa == 1
        B(:,:,i) = A((center(2)-7):(center(2)+7),(center(1)-7):(center(1)+7));
        imagesc(X_NP_mem(i)-smallboxwidth/2:X_NP_mem(i)+smallboxwidth/2,...
            Y_NP_mem(i)-smallboxwidth/2:Y_NP_mem(i)+smallboxwidth/2,B(:,:,i));
        hold on;
        plot(center(1),center(2),'ko','LineWidth',4)
        hold off;
        drawnow
%         saveas(gcf,[fold 'images\' num2str(i) '_2Dfit.jpg'],'jpg')
        saveas(gcf,[fold num2str(i) '_2Dfit.jpg'],'jpg')
    else
%         disp(['i = ' num2str(i)])
        bbbb = input('Do you want to manually pick the maximums (1 = yes, 0 = escape the spot): ');
        if bbbb == 1
            X = input('what is X for maximum? (where X from imagesc plot): ');
            Y = input('what is Y for maximum? (where Y from imagesc plot): ');
            center = [X Y];
            B(:,:,i) = A((center(2)-7):(center(2)+7),(center(1)-7):(center(1)+7));
            % draw the new loc
            figure(6)
            imagesc(X_NP_mem(i)-smallboxwidth/2:X_NP_mem(i)+smallboxwidth/2,...
                Y_NP_mem(i)-smallboxwidth/2:Y_NP_mem(i)+smallboxwidth/2,B(:,:,i));
%             hold on
%             plot(center(1),center(2),'ko','LineWidth',4)
%             hold off
            drawnow
%             saveas(gcf,[fold 'images\' num2str(i) 'handpicked.jpg'],'jpg')
            saveas(gcf,[fold num2str(i) 'handpicked.jpg'],'jpg')
        
        else
            continue
        end
    end
    clear cropped
%     save([fold 'images\cropped diff limited membrane.mat'])
end

% make mask to find I_1, I_3, I_5, I_7
a = nan(15,15,size(B,3));
b = a; b(7:9,7:9,:) = 1;
c = a; c([5 11],7:9,:) = 1;c(7:9,[5 11],:) = 1; c([6 10],[6 10],:) = 1;
d = a; d([3 13],6:10,:) = 1;d(6:10,[3 13],:) = 1; d([4 12],[5 11],:) = 1; d([5 11],[4 12],:) = 1;
e = a; e([1 15],6:10,:) = 1;e(6:10,[1 15],:) = 1;e([4 5 11 12],[2 14],:) = 1;e([2 14],[4 5 11 12],:) = 1;e([3 13],[3 13],:)=1;
I1 = B.*b;
I2 = B.*c;
I3 = B.*d;
I4 = B.*e;
avgI1 = mean(I1(~isnan(I1)));
avgI2 = mean(I2(~isnan(I2)));
avgI3 = mean(I3(~isnan(I3)));
avgI4 = mean(I4(~isnan(I4)));
stdI1 = std(I1(~isnan(I1)));
stdI2 = std(I2(~isnan(I2)));
stdI3 = std(I3(~isnan(I3)));
stdI4 = std(I4(~isnan(I4)));
% ave_of_atimesb = sum(sum(a.*b))/sum(sum(b))
I1(isnan(I1)) = [];
I1 = reshape(I1,9,size(I1,2)/9);
I1 = mean(I1);

I3(isnan(I3)) = [];
I3 = reshape(I3,28,size(I3,2)/28);
I3 = mean(I3);

figure(7)
errorbar([0 3:2:7],[avgI1 avgI2 avgI3 avgI4],[stdI1 stdI2 stdI3 stdI4],'ko','LineWidth',3)
xlabel('distance away from center (px)','fontsize',30)
ylabel('<I>_{all particles} (photon)','fontsize',30)
title(['I vs location'],'fontsize',30)
set(gca,'fontsize',30,'linewidth',4)
xlim([-0.5 7.5])
saveas(gcf,[fold 'images\I vs distance.jpg'],'jpg')
save([fold 'images\cropped diff limited membrane.mat']);
save(['C:\Users\Zoe\Google Drive\paper 3\diffraction image intensity study\POPC TF\' name 'cropped diff limited membrane.mat'],'B');