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
