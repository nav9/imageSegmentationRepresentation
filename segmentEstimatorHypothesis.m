%%Image Segmentation via Estimation
clc;clear all;close all;
%% Initializations
spr=4; spc=1; spi=1;%subplot vars
paramVarThresh = 0.95;%parameter error variation threshold
figNum = 1;
figure(figNum);figNum = figNum + 1;
alpha = 0.6799;%alpha = 0.5359;%alpha = 0.5199;

for roc = 1:1
zForAlpha = norminv(alpha);    
fprintf('ROC %d with alpha=%f\n', roc, alpha);
OI = imread('bd/87015.jpg');T = imread('bd/GT87015.png');nam='Snake';
% OI = imread('bd/100098.jpg');T = imread('bd/GT100098.png');nam='Bear';
% OI = imread('bd/69022.jpg');T = imread('bd/GT69022.png');nam='Kangaroo';
% OI = imread('bd/Plane.jpg');T = imread('bd/GTPlane.png');nam='Plane';
I = double(rgb2gray(OI)); T = double(T);
HI = rgb2hsv(OI);
D = I;%stores difference in estimated intensity vs actual intensity
Y = D;%stores estimates
I_FINAL = I; X1 = I; X2 = I; X3 = I;
ir = size(I,1); ic = size(I,2);%image size
psv = 99;%patch size vertical (use only odd numbers)
psh = 99;%patch size horizontal (use only odd numbers)
psvw = (psv-1)/2; pshw = (psh-1)/2;%patch width on either side
if psvw < 1, psvw=1;end; if pshw < 1, pshw=1;end;%bounds check
%create a patch of i's and j's for use in A
isr = linspace(1, ir, ir)'; jsc = linspace(1, ic, ic); iis = [];jjs = [];
for i = 1:ic;iis = [iis isr];end; for i = 1:ir;jjs = [jjs; jsc];end; %#ok<*AGROW>

% subplot(spr, spc, spi);spi=spi+1;imshow(uint8(I));title('original');
% subplot(spr, spc, spi);spi=spi+1;imshow(uint8(T));title('ground truth');
% subplot(spr, spc, spi);spi=spi+1;plot3(iis(:), jjs(:), I(:));title('3D intensity');
if psv > ir || psh > ic || psv < 1 || psh < 1, fprintf('PATCH OUT OF BOUNDS\n');return;end
if ~mod(psv, 2) || ~mod(psh, 2), fprintf('USE ONLY ODD NUMBERS FOR PATCH\n');return;end

popMean = mean(I(:)); popVar = var(I(:)); popStd = std(I(:));
foreBackRatio = hist(T(:), [255 0]);
foreBackRatio = foreBackRatio(1) / foreBackRatio(2);%prior probability of foreground pix
histFG=zeros(1, 256); histBG=histFG; histTP=histFG; histTN=histFG; histFP=histFG; histFN=histFG;
%% Processing
highestError = -10000;
% zStore = [];xa=0;xb=0;xc=0;%deleteme
% anim = [];
tic
for i = 1:ir
    err = zeros(1, ic);
    priorProb = foreBackRatio;
    for j = 1:ic
        pvs=i-psvw; pve=i+psvw; phs=j-pshw; phe=j+pshw;%patch st & en pos
        if pvs<1,pvs=1;end; if phs<1,phs=1;end;%bounds check
        if pve>ir,pve=ir;end; if phe>ic,phe=ic;end;%bounds check
        
        IP = I(pvs:pve, phs:phe);
        Ai = iis(pvs:pve, phs:phe); Ai = Ai(:);
        Aj = jjs(pvs:pve, phs:phe); Aj = Aj(:);
        A = [Ai Aj ones(length(Ai), 1)];
        B = I(pvs:pve, phs:phe); B = B(:);
        %--- Least squares
%         W = (IP-mean(mean(IP)).*ones(size(IP,1),size(IP,2))).^2;%variance
%         W = 1./W; W = diag(W(:));%1/variance = weight          
        
        X = (A'*A)\(A'*B);%X = inv(A'*A)*(A'*B);%X = linsolve(A,B);
        X1(i,j)=X(1);X2(i,j)=X(2);X3(i,j)=X(3);
        Y(i,j) = i*X1(i,j) + j*X2(i,j) + X3(i,j);%linear        
        D(i,j) = Y(i,j) - I(i,j);        
    end
    %---Hypothesis: Take mean after eliminating outliers and compare
    if i==1, dV = D(i,:); else dV = [D(i,:) D(i-1,:)];end
    maV = max(dV); miV = min(dV); haV = maV - miV; indV = find(dV > miV+haV);
    dV(indV) = [];%eliminate half the outlier error values
    popMean = mean(dV(:));
    for j=1:ic
        zScore = (D(i,j) - popMean) / popStd;
        if zForAlpha > zScore, 
            I_FINAL(i,j) = 0; histBG(I(i,j)+1) = histBG(I(i,j)+1) + 1;
        else
            I_FINAL(i, j) = 255; histFG(I(i,j)+1) = histFG(I(i,j)+1) + 1;
        end
        %---Sensitivity Specificity
        if I_FINAL(i,j)==0 && T(i,j)==0, histTN(I(i,j)+1)=histTN(I(i,j)+1)+1;end
        if I_FINAL(i,j)==255 && T(i,j)==255, histTP(I(i,j)+1)=histTP(I(i,j)+1)+1;end
        if I_FINAL(i,j)==0 && T(i,j)==255, histFN(I(i,j)+1)=histFN(I(i,j)+1)+1;end
        if I_FINAL(i,j)==255 && T(i,j)==0, histFP(I(i,j)+1)=histFP(I(i,j)+1)+1;end
    end
%     if i==159,plot(linspace(1,size(D(i,:),2),size(D(i,:),2)),D(i,:),'k');hold on;end%residual    
end

% plot(linspace(1,ic,ic), D(159, :),'b');

% %% R square
% totSumSq = sum((Y(:) - mean(Y(:)).*ones(size(Y(:)))).^2);
% SumSqRes = sum(D(:).^2);
% rSq=1-(SumSqRes/totSumSq);
% fprintf('R squared = %f\n', rSq);

% %% Thresholding
% tm = mean(D(:));
% I_FINAL(D > tm*paramVarThresh) = 255;
% I_FINAL(D < tm*paramVarThresh) = 0;
fprintf('Time taken: %fs\n', toc);

% plot(linspace(1,length(zStore),length(zStore))', zStore, 'b');hold on;
% plot(anim, 0, 'r');hold off;
% plot(linspace(1,length(histFG),length(histFG)),histBG,linspace(1,length(histBG),length(histBG)), histFG);
% legend('foreground', 'background'); title(strcat(nam,'. Non linear regression result')); xlabel('intensities'); ylabel('number of pixels');

% %%Errors
% figure(figNum);figNum = figNum + 1; mesh(D);title(strcat(nam,': Errors'));

%%Sensitivity Specificity
figure(figNum);figNum = figNum + 1;
plot(linspace(1,length(histTP),length(histTP)),histTP,  linspace(1,length(histTN),length(histTN)), histTN,  linspace(1,length(histFP),length(histFP)), histFP,  linspace(1,length(histFN),length(histFN)), histFN);
legend('TP', 'TN', 'FP', 'FN'); title(strcat(nam,'. Sensitivity Specificity')); xlabel('intensities'); ylabel('number of pixels');
%% Dice Score
dice = 2*nnz(I_FINAL&T)/(nnz(I_FINAL)+nnz(T));
fprintf('Dice score difference between ground truth and segmented image:\n %f%%\n', dice*100); 
        
figure(figNum);figNum = figNum + 1;
imshow(uint8(I_FINAL));title(strcat('Shows parameter variation points. Error threshold at ', num2str(paramVarThresh)));

% figure(figNum);figNum = figNum + 1;
% imhist(D);
sum(histFP)
sum(histTP)
%%Plot ROC curve
figure(figNum);figNum = figNum + 1;
plot(sum(histFP), sum(histTP), '*r');hold on;title('ROC curve');xlabel('1-specificity=FP=Pfa');ylabel('sensitivity=TP=Pd');
alpha=alpha+0.04;
end % roc
