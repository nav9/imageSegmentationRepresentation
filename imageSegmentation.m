clc;clear all;close all;

% %---Segmentation via thresholding
% spr=1;spc=1;spi=1;%subplot vars
% figure(1);
% for i=1:1
%     filename = strcat('E:\Module4_ImageProcessing\Assignment\peds\frame_',num2str(i),'.jpg');
%     I = imread(filename);
%     tic
%     %-----preprocess
%     I = histeq(I);II=I;
%     I = imopen(I, [1 1 1;1 1 1;1 1 1]);
%     I = imsharpen(I);
%     %-----intensity slice
%     I(I<25) = 0;
%     I(I>180) = 255;
%     I(I>165 & I<181) = 255;
%     I(I>40 & I<130) = 80;  
%     I(I>130 & I<150 | I == 130) = 80; 
%     I(I>151 & I<166) = 255;  
%     fprintf('Time taken: %fs\n', toc);
%     subplot(spr, spc, spi);spi=spi+1;imshow(I);title(num2str(i));
%     if i==1,
%         G = imread(strcat('E:\Module4_ImageProcessing\Assignment\peds\groundTruth.bmp'));
%         G = rgb2gray(G);
%         ge = entropy(G); ie = entropy(I);
%         fprintf('Score difference between ground truth and segmented image: %f', sum(sum(ie-ge)));
%     end    
% end

%---quadtree
spr=1;spc=1;spi=1;%subplot vars
meanT1 = 106;%threshold divider for sigma
sigBright = 28; sigDark = 35;%acceptable sigma for homogenity in block
finBright = 255; finDark = 0;%final values to set
markingThresh = 122;%threshold below which values are marked 0
figure(1);
lev2=[];
for img = 1:1
    I = imread(strcat('E:\Module4_ImageProcessing\Assignment\peds\frame_',num2str(img),'.jpg'));
    tic
    I = histeq(I);
    II = -1.*ones(size(I));
    imgR = size(I,1); imgC = size(I,2);
    lev = [1 1 imgR imgC];%quadtree stack. Topmost level
    lev2 = [1 1 imgR imgC];%Just to keep record. [r,c,rowHt,colWidth]
    
    while size(lev,1) > 0,%quadtree creation and traversal       
        sq = lev(size(lev,1), :);%read last row of lev
        r=sq(1);c=sq(2);h=sq(3);w=sq(4);
        if II(r,c)>=0, 
            lev(size(lev,1), :) = [];%remove from stack
            continue;
        end%already marked
        SI = double(I(r:r+h-1, c:c+w-1));%image section
        me = mean(SI(:));
        st = std(SI(:));
        if me < meanT1, threshold = sigDark; else threshold = sigBright;end
        if st > threshold && h > 1 && w > 1,%if can split into 4
            hr = floor(h/2);
            hc = floor(w/2);
            nl = [r c hr hc;
                  r+hr c h-hr hc;
                  r c+hc hr w-hc;
                  r+hr c+hc h-hr w-hc];
            lev = [lev; nl]; %#ok<AGROW>   
            lev2 = [lev2; nl]; %#ok<AGROW> 
        else          
            if me < markingThresh, II(r:r+h, c:c+w) = finDark; 
            else II(r:r+h, c:c+w) = finBright;end
            %II(r:r+h, c:c+w) = me;%mark mean
            lev(size(lev,1), :) = [];%remove from stack
        end
    end
    fprintf('Time taken: %fs\n', toc);
    subplot(spr, spc, spi);spi=spi+1;imshow(uint8(II));title(num2str(img));
    if img==1,
        G = imread(strcat('E:\Module4_ImageProcessing\Assignment\peds\groundTruth.png'));
        G = rgb2gray(G);
        ge = entropy(G); ie = entropy(II);
        fprintf('Score difference between ground truth and segmented image:\n %f', sum(sum(ie-ge)));
        dice = 2*nnz(II&G)/(nnz(II)+nnz(G));
        fprintf('Dice score difference between ground truth and segmented image:\n %f\n', dice);        
    end
end
%---display quadtree extents
figure(2);
I=imread(strcat('E:\Module4_ImageProcessing\Assignment\peds\frame_1.jpg'));
II = cat(3,I,I,I);
for i=1:size(lev2,1)
    r=lev2(i,1)+lev2(i,3);
    c=lev2(i,2)+lev2(i,4);
    II(lev2(i,1):r,c,2)=255;
    II(r,lev2(i,2):c,2)=255;
end
imshow(II);title('quadtree iterations');

