%---pixel adjacency and connectivity
%TODO: check for possibility of bug where one or more pixels are outside
%image boundary for 4 or 8 adjacency

function pixelAdjacency() 
    clc; close all; clear all;
    whichAdj = 2;%1=4 adjacency, 2=8 adjacency, 3=m adjacency
    
%     I = imread('moon.tif'); 
%     II = I(200:240, 30:300); %II = im2bw(I(200:220, 30:300)); 
    II = [1 1 1 0 0 0 0 0;
          0 0 0 1 0 0 0 0;
          0 0 0 0 0 0 1 1;];
    A = -1*zeros(size(II));
    label = 1;%starting label. Should be greater than 0 because of the if conditions in adjacency functions
    [ir, ic] = size(II);
    s = [1:200];%set
    
    %---first pass
    for i=1:ir,
        for j=1:ic,
            if isempty(s(s==II(i,j))), A(i,j) = 0;%p is not a member of set
            else
                if whichAdj==1, [foundNeighbour, A] = adj4(i,j,s,II,A);end
                if whichAdj==2, [foundNeighbour, A] = adj8(i,j,s,II,A);end
                if ~foundNeighbour, label = label + 1; A(i,j) = label;end%couldn't find adjacent neighbour of same set so assigning new label
            end
        end
    end
    
    %---second pass
    for i=1:ir,
        for j=1:ic,
            if A(i,j)==0 || isempty(s(s==II(i,j))), continue;end%nothing to consolidate || p is not member of set
            if whichAdj==1, [A] = adj4Consolidate(II, A, s, i, j);end
            if whichAdj==2, [A] = adj8Consolidate(II, A, s, i, j);end
        end
    end
    
    %---plot
%     sc = 1; sr = 2; si = 1;
%     subplot(sc,sr,si);si=si+1; imshow(II); title('original');
%     subplot(sc,sr,si);si=si+1; imshow(uint8(A)); title('adjacency adjusted');
    fprintf('%d labels in 1st pass \n%d labels after 2nd pass\n', label-1, length(unique(A(:)))-1);
    disp('Adjacency adjusted matrix');disp(A);
end

function [found, A] = adj4(r, c, s, I, A)%check for 4 adjacency
    found = 0;    
    [ir, ic] = size(I);    
    if r-1 >= 1, q = I(r-1, c); if ~isempty(s(s==q)), if A(r-1, c) > 0, A(r, c) = A(r-1, c); found = 1;return;end;end;end%if within limits, if q belongs to set, if existing label > 0, assign label of neighbour
    if c-1 >= 1, q = I(r, c-1); if ~isempty(s(s==q)), if A(r, c-1) > 0, A(r, c) = A(r, c-1); found = 1;return;end;end;end%if within limits, if q belongs to set, if existing label > 0, assign label of neighbour
    if c+1 <= ic, q = I(r, c+1); if ~isempty(s(s==q)), if A(r, c+1) > 0, A(r, c) = A(r, c+1); found = 1;return;end;end;end%if within limits, if q belongs to set, if existing label > 0, assign label of neighbour
    if r+1 <= ir, q = I(r+1, c); if ~isempty(s(s==q)), if A(r+1, c) > 0, A(r, c) = A(r+1, c); found = 1;return;end;end;end%if within limits, if q belongs to set, if existing label > 0, assign label of neighbour    
end

function [found, A] = adj8(r, c, s, I, A)%check for 8 adjacency
    [ir, ic] = size(I);
    [found, A] = adj4(r,c,s,I,A);
    if ~found, return;end
    %---check diagonals
    if r-1 >= 1 && c-1 >= 1, q = I(r-1, c-1); if ~isempty(s(s==q)), if A(r-1, c-1) > 0, A(r, c) = A(r-1, c-1); found = 1;return;end;end;end%if within limits, if q belongs to set, if existing label > 0, assign label of neighbour
    if r+1 <= ir && c-1 >= 1, q = I(r+1, c-1); if ~isempty(s(s==q)), if A(r+1, c-1) > 0, A(r, c) = A(r+1, c-1); found = 1;return;end;end;end%if within limits, if q belongs to set, if existing label > 0, assign label of neighbour
    if r-1 >= 1 && c+1 <= ic, q = I(r-1, c+1); if ~isempty(s(s==q)), if A(r-1, c+1) > 0, A(r, c) = A(r-1, c+1); found = 1;return;end;end;end%if within limits, if q belongs to set, if existing label > 0, assign label of neighbour
    if r+1 <= ir && c+1 <= ic, q = I(r+1, c+1); if ~isempty(s(s==q)), if A(r+1, c+1) > 0, A(r, c) = A(r+1, c+1); found = 1;return;end;end;end%if within limits, if q belongs to set, if existing label > 0, assign label of neighbour    
end

function [A] = adj4Consolidate(I, A, s, r, c)    
    [ir, ic] = size(I);
    if r-1 >= 1, q = I(r-1, c); if ~isempty(s(s==q)), if A(r-1, c) > 0, A(r-1, c) = A(r,c); end;end;end%if within limits, if q belongs to set, if existing label > 0, assign label to neighbour
    if c-1 >= 1, q = I(r, c-1); if ~isempty(s(s==q)), if A(r, c-1) > 0, A(r, c-1) = A(r, c); end;end;end%if within limits, if q belongs to set, if existing label > 0, assign label to neighbour
    if c+1 <= ic, q = I(r, c+1); if ~isempty(s(s==q)), if A(r, c+1) > 0, A(r, c+1) = A(r, c); end;end;end%if within limits, if q belongs to set, if existing label > 0, assign label to neighbour
    if r+1 <= ir, q = I(r+1, c); if ~isempty(s(s==q)), if A(r+1, c) > 0, A(r+1, c) = A(r, c); end;end;end%if within limits, if q belongs to set, if existing label > 0, assign label to neighbour    
end

function [A] = adj8Consolidate(I, A, s, r, c)
    if A(r,c)==0 || isempty(s(s==I(r,c))), return;end%nothing to consolidate || p is not member of set
    [A] = adj4Consolidate(I, A, s, r, c);
    [ir, ic] = size(I);
    if r-1 >= 1 && c-1 >= 1, q = I(r-1, c-1); if ~isempty(s(s==q)), if A(r-1, c-1) > 0, A(r-1, c-1) = A(r,c); end;end;end%if within limits, if q belongs to set, if existing label > 0, assign label to neighbour
    if r+1 <= ir && c-1 >= 1, q = I(r+1, c-1); if ~isempty(s(s==q)), if A(r+1, c-1) > 0, A(r+1, c-1) = A(r, c); end;end;end%if within limits, if q belongs to set, if existing label > 0, assign label to neighbour
    if r-1 >= 1 && c+1 <= ic, q = I(r-1, c+1); if ~isempty(s(s==q)), if A(r-1, c+1) > 0, A(r-1, c+1) = A(r, c); end;end;end%if within limits, if q belongs to set, if existing label > 0, assign label to neighbour
    if r+1 <= ir && c+1 <= ic, q = I(r+1, c+1); if ~isempty(s(s==q)), if A(r+1, c+1) > 0, A(r+1, c+1) = A(r, c); end;end;end%if within limits, if q belongs to set, if existing label > 0, assign label to neighbour    
end

% Two pixels p and q with values from V are m-adjacent if : 
% q is in N4(p) or
% q is in ND(p) and the set N4(p) intersection N4(q) has no pixel whose values are from V (no intersection)
