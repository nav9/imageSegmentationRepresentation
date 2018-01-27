clc;clear all;close all;

%%Estimation
spr=10;spc=3;spi=1;%subplot vars
paramVarThresh = 800;%parameter error variation threshold
figure(1);
I = imread('E:\Module5_StatisticalEstimation\assignment\imageSegmentationRepresentation\bd\69022.jpg');
T = imread('E:\Module5_StatisticalEstimation\assignment\imageSegmentationRepresentation\bd\GT69022.png');
I = double(rgb2gray(I)); T = double(T);
IVAR = I;
subplot(spr, spc, spi);spi=spi+1;imshow(uint8(I));title('original');
subplot(spr, spc, spi);spi=spi+1;imshow(uint8(T));title('ground truth');
ir = size(I,1); ic = size(I,2);%image size
psv = 3; psh = 3;%patch size vertical and patch size horizontal
if psv > ir || psh > ic, disp('PATCH TOO LARGE');end

highestError = -10000;
for i = round(psv/2):ir-round(psv/2)
    err = [];
    for j = round(psh/2):ic-round(psh/2)
        if psv > 2, phv = floor(psv/2); else phv=0;end;if psh > 2, phh = floor(psh/2); else phh=0;end
        IP = I(i-phv:i+phv, j-phh:j+phh);
        A = []; B = [];
        for k=i-phv:i+phv
            for kk=j-phh:j+phh
                A = [A; k kk 1];
                B = [B; I(k,kk)];                
            end
        end
        X = linsolve(A,B);
        dif = X(3) - I(i,j);
        err = [err dif]; %#ok<*AGROW>
        %if highestError < X(3), highestError = X(3);fprintf('highestError = %d at [%d,%d]\n', highestError,i,j);end
        if dif > paramVarThresh, IVAR(i,j) = 255; else IVAR(i,j) = 0;end
%         subplot(spr, spc, spi);spi=spi+1;imshow(uint8(IP));
    end
    if spi < (spr*spc)-2, subplot(spr, spc, spi);spi=spi+1;plot(linspace(1,length(err),length(err)),err);end
end

figure(2);
imshow(uint8(IVAR));title(strcat('Shows parameter variation points. Error threshold at ', num2str(paramVarThresh)));

figure(3);
spr=2;spc=2;spi=1;%subplot vars
H = imread('E:\Module5_StatisticalEstimation\assignment\imageSegmentationRepresentation\bd\69022.jpg');
H = rgb2hsv(H);
subplot(spr, spc, spi);spi=spi+1;imshow(H(:,:,1));
subplot(spr, spc, spi);spi=spi+1;imshow(H(:,:,2));
subplot(spr, spc, spi);spi=spi+1;imshow(H(:,:,3));