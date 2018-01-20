clc;clear all;close all;

% %---object description
% spr=1;spc=1;spi=1;%subplot vars
% figure(1);
% for i=1:1
%     fprintf('%d --------------\n',i);
%     I = imread(strcat('E:\Module4_ImageProcessing\Assignment\peds\frame_',num2str(i),'.jpg'));
%     tic
%     I = imsharpen(histeq(I));
%     I = double(I);
%     I = edge(I, 'canny', 0.7, 0.9);%threshold,sigma
% %     I = bwmorph(I, 'bridge');
% %     I = edge(I, 'canny', 0.7, 0.1);
%     I = bwmorph(I, 'spur');I = bwmorph(I, 'spur');%using twice removes spurs better
%     I = bwmorph(I, 'thin');
%     %---trace edges
%     I = -1.*I;%mark all 1's as -1 so that markers below can function well
%     marker=1;
%     I(1,:)=0;I(size(I,1),:)=0;I(:,1)=0;I(:,size(I,2))=0;
%     for r=2:size(I,1)-1%-1 prevents edge tracer window from breaching boundary
%         for c=2:size(I,2)-1
%             if I(r,c)==-1,%found an edge point
%                 %initialize tracing edge
%                 foundMore=1;nr=r;nc=c;
%                 while foundMore,
%                     I(nr,nc) = marker;
%                     if I(nr+1,nc)==-1,nr=nr+1;%7
%                     elseif I(nr,nc+1)==-1,nc=nc+1;%1
%                     elseif I(nr-1,nc+1)==-1,nr=nr-1;nc=nc+1;%2
%                     elseif I(nr-1,nc)==-1,nr=nr-1;%3
%                     elseif I(nr+1,nc+1)==-1,nr=nr+1;nc=nc+1;%8
%                     elseif I(nr,nc-1)==-1,nc=nc-1;%5
%                     elseif I(nr+1,nc-1)==-1,nr=nr+1;nc=nc-1;%6
%                     elseif I(nr-1,nc-1)==-1,nr=nr-1;nc=nc-1;%4
%                     else foundMore=0;marker = marker + 1;
%                     end
%                 end
%             end
%         end
%     end
%     fprintf('%d objects traced\n', marker);
%     %---remove redundancies
%     joinSize=4;%joined=[0];
%     for r=2:size(I,1)-joinSize%-5 prevents edge tracer window from breaching boundary
%         for c=2:size(I,2)-joinSize
%             %if ~ismember(I(r,c), joined),%found an edge point
%             if I(r,c)>0,
%                 v = [I(r,c)]; %#ok<*NBRAK>
%                 for k=1:joinSize
%                     v = [v I(r+k, c)]; %#ok<*AGROW>
%                     v = [v I(r+k, c+k)];
%                     v = [v I(r, c+k)];
%                 end
%                 v = unique(v);v(v < 1) = [];
%                 if length(v) > 1,
%                     for k=2:length(v)
%                         I(I==v(k))=v(1);
%                     end
%                 end
%             end
%         end
%     end
%     uni = unique(I); uni(uni<1)=[];%remove zero
%     
%     fprintf('%d objects after cleaning\n', length(uni));
%     I = uint8(I);
%     maxPedWd = 17;minPedWd=5;maxPedHt=7;
%     %----- Representation
%     for k=1:length(uni)
%         %[r, c] = ind2sub(size(I), find(I, uni(k))); %#ok<*SAGROW>
%         obj=[];%get each pedestrian's coordinate extents
%         peri = [];%minCol,maxCol,row
%         for r=1:size(I,1)
%             maC=0;miC=size(I,2);
%             for c=1:size(I,2)
%                 if I(r,c)==uni(k),
%                     obj=[obj; [r, c]];
%                     if c > maC, maC = c;end
%                     if c < miC, miC = c;end                      
%                 end            
%             end
%             if maC>0,peri = [peri; miC maC r];end
%         end
%         %remove small and wide objects
%         objWidth = abs(max(obj(:,2))-min(obj(:,2)));
%         objHeight = abs(max(obj(:,1))-min(obj(:,1)));
%         if objWidth > maxPedWd || objWidth < minPedWd || objHeight < maxPedHt, 
%             for r=1:size(obj,1),I(obj(r,1),obj(r,2))=0;end;continue;
%         end        
%         %extract only boundary values
%         finalB = [];
%         for kk=1:size(peri,1)
%             for c=peri(kk,1):peri(kk,2)
%                 if c~=peri(kk,1) && c~=peri(kk,2), I(peri(kk,3),c) = 0;end
%             end
%         end
%         %fprintf('%fs\n',toc);
%         %----get the chain of dots
%         tic
%         found=0;%find the first dot
%         for r=1:size(I,1)
%             for c=1:size(I,2),if I(r,c)==uni(k),found=1;break;end;end
%             if found, break;end
%         end
%         traceDir = 1;%1=cw, 2=ccw, 3=stop
%         dir = 0;%0=st,1=e,2=ne,3=n,4=nw...8=se
%         cw=[r,c,dir]; ccw=[]; T=zeros(size(I));
%         while traceDir < 3,
%             T(r,c)=1;
%             if I(r+1,c)==uni(k) && T(r+1,c)==0,r=r+1;dir =7;%7
%             elseif I(r,c+1)==uni(k) && T(r,c+1)==0,c=c+1;dir=1;%1
%             elseif I(r-1,c+1)==uni(k) && T(r-1,c+1)==0,r=r-1;c=c+1;dir=2;%2
%             elseif I(r-1,c)==uni(k) && T(r-1,c)==0,r=r-1;dir=3;%3                
%             elseif I(r+1,c+1)==uni(k) && T(r+1,c+1)==0,r=r+1;c=c+1;dir=8;%8
%             elseif I(r,c-1)==uni(k) && T(r,c-1)==0,c=c-1;dir=5;%5
%             elseif I(r+1,c-1)==uni(k) && T(r+1,c-1)==0,r=r+1;c=c-1;dir=6;%6
%             elseif I(r-1,c-1)==uni(k) && T(r-1,c-1)==0,r=r-1;c=c-1;dir=4;%4
%             else if traceDir==1,r=cw(1,1);c=cw(1,2);end; traceDir = traceDir + 1;
%             end
%             if traceDir==1, cw = [cw; [r,c,dir]];end
%             if traceDir==2, ccw = [ccw; [r,c,dir]];end
%         end
%         cw = [cw; ccw(size(ccw,1):-1:1,:)];%chaincode
% %          fprintf('%d: ',k);for fp=1:size(cw,1),fprintf('%d,',cw(fp,3));end;fprintf('\n');
% %    
%         miR = min(cw(:,1)); maR = max(cw(:,1)); miC = min(cw(:,2)); maC = max(cw(:,2));
%         cen = [miR+(maR-miR)/2 miC+(maC-miC)/2];
% %          fprintf('centroid = %d,%d\n', round(cen(1)), round(cen(2)));
%         
%         %---signature
%          sig = [];sigR=[];%fprintf('Sig: %d: ',k);
%         for kk=1:size(cw,1)
%             sig = [sig; (180/pi)*atan(abs(cen(1)-cw(kk,1))/abs(cen(2)-cw(kk,2)))];
%             sigR = [sigR; sqrt((cen(1)-cw(kk,1))^2+(cen(2)-cw(kk,2))^2)];
%         end
%         fprintf('%fs\n',toc);
% %         for fp=1:size(sig,1),fprintf('%d,',round(sig(fp)));end;fprintf('\n');
% %         for fp=1:size(sig,1),fprintf('%d,',round(sigR(fp)));end;fprintf('\n');
%     end
%     
%     %I = 3.*uint8(I)-24;
%     %I = cat(3, I,ones(size(I)),ones(size(I)));
%     I=9.*I;
%     subplot(spr, spc, spi);spi=spi+1;imshow(I);title(num2str(i));
% end

