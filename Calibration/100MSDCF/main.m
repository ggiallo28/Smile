close all; clear all; clc;
% A = [];
% for i = 2: 7
%     orig = imread(['_DSC007',num2str(i),'.JPG']);
%     I = imcrop(orig);
%     R = double(I(:,:,1));
%     G = double(I(:,:,2));
%     B = double(I(:,:,3));
%     lumMat = (0.2126 * R) + (0.7152 * G) + (0.0722 * B); % 
%     K = zeros(size(I));
%     K(:,:,1) = R./(255-lumMat);
%     K(:,:,2) = G./(255-lumMat);
%     K(:,:,3) = B./(255-lumMat);
%     K = im2uint8(K);
%     A = [A;K];
% end
figure,imshow(imread('checkerboard.jpg'));
checker_vector = reshape([[0,0,0;255,0,255];[0,0,0;0,255,255];[0,0,0;255,255,0];[255,255,255;255,0,0];[255,255,255;0,255,0];[255,255,255;0,0,255]],[2,6,3]);
checker_center = [0.5*size(checker_vector,2),0.5*size(checker_vector,2)+1];
orig = imread(['_DSC007',num2str(2),'.JPG']);
figure, O = imcrop(orig);
% [edge_magnitude, edge_orientation] = coloredges(I);
% img1=rgb2ycbcr(I);
% dx1=edge(img1(:,:,1),'canny');
% med = medfilt2(edge_magnitude,[5,5]);
% res = imclose(dx1,strel('disk',5));
% res = median(median(med)).*im2double(res1);
% med = res + med;
% T = adaptthresh(edge_magnitude,0.3,'ForegroundPolarity','bright');
% BW = imbinarize(edge_magnitude,T);
% BW = medfilt2(BW,[9,9]);
% imshow([BW;med;T]);
cymk = rgb2cmyk(O);
R = im2double(O(:,:,1));
G = im2double(O(:,:,2));
B = im2double(O(:,:,3));
Cyano = im2double(cymk(:,:,1));
Magenta = im2double(cymk(:,:,2));
T = adaptthresh(Magenta,0,'ForegroundPolarity','bright');
BW = imbinarize(Magenta,T);
figure
imshow(BW)
[Gx, Gy] = imgradientxy(R,'central');
figure,
imshowpair(Gx,Gy,'falsecolor');
Yellow = im2double(cymk(:,:,3));
Key = im2double(cymk(:,:,4));
lumMat = (0.2126 * R) + (0.7152 * G) + (0.0722 * B); % 
K = zeros(size(O));
K(:,:,1) = R./(1-lumMat);
K(:,:,2) = G./(1-lumMat);
K(:,:,3) = B./(1-lumMat);
% K(:,:,1) = R./(R+G+B);
% K(:,:,2) = G./(R+G+B);
% K(:,:,3) = B./(R+G+B);
I = im2uint8(K);
% nColors = 7;
% kmeans: http://stackoverflow.com/questions/25691735/matlab-color-based-segmentation
%% Segmentazione
obj_red = createMask(I,'red');
obj_green = createMask(I,'green');
obj_blue = createMask(I,'blue');
obj_yellow = createMask(I,'yellow');
obj_chess = [obj_red;obj_green;obj_blue;obj_yellow];
%% Plot Results
inv_BW = uint8(obj_red.black_white|obj_green.black_white|obj_blue.black_white|obj_yellow.black_white);
col_BW = uint8(obj_red.color_mask|obj_green.color_mask|obj_blue.color_mask|obj_yellow.color_mask);
BW = 255*(inv_BW-col_BW);
fuse = obj_red.masked_rgb+obj_green.masked_rgb+obj_blue.masked_rgb+obj_yellow.masked_rgb+repmat(BW,[1 1 3]);
imshow(fuse); hold on;
for l=1:size(obj_chess,1)
    obj = obj_chess(l);
    for i = 1:size(obj.chess,2)
        for k=1:size(obj.chess(i).center_x,1)
            for j=1:size(obj.chess(i).center_x,2)
                scatter(obj.chess(i).center_x(k,j),obj.chess(i).center_y(k,j)) % Riferimento Y verso il basso
            end
        end
    end
end
hold off;
%% Separa i riflessi
enlarged = imdilate(fuse,strel('disk',20));
CC = bwconncomp(enlarged);
maskCenter = false(size(fuse)); maskCenter(CC.PixelIdxList{2}) = 1;  maskCenter = maskCenter(:,:,1) | maskCenter(:,:,2)| maskCenter(:,:,3);
maskLeft = false(size(fuse)); maskLeft(CC.PixelIdxList{1}) = 1;  maskLeft = maskLeft(:,:,1) | maskLeft(:,:,2) | maskLeft(:,:,3);
maskRight = false(size(fuse)); maskRight(CC.PixelIdxList{3}) = 1;  maskRight = maskRight(:,:,1) | maskRight(:,:,2) | maskRight(:,:,3);
positions = [{'Left'} {'Center'} {'Right'}];
types = [{'Primary'} {'Secondary'} {'Real'}];
order = []; k =1;
% Order Vector of chess
for l=1:size(obj_chess,1)
    for i = 1:size(obj_chess(l).chess,2)
        order = [order,[l;i;obj_chess(l).chess(i).centroid(1)]]; k = k+1;
    end
end
% Find Objects Position
for l=1:size(obj_chess,1)
    left = 0;
    right = 0;
    arr = [];
    for i = 1:size(obj_chess(l).chess,2)
        v(1) = sum(sum(obj_chess(l).chess(i).mask&maskLeft));
        v(2) = sum(sum(obj_chess(l).chess(i).mask&maskCenter));
        v(3) = sum(sum(obj_chess(l).chess(i).mask&maskRight));
        idx = find(v==max(v));
        obj_chess(l).chess(i).position = positions(idx);
        if idx==2
            obj_chess(l).chess(i).type =types(3);
        end
    end
end
% Adds label
label = cell(3,k-1); k =1;
for l=1:size(obj_chess,1)
    for i = 1:size(obj_chess(l).chess,2)
        label(1,k) = cellstr(obj_chess(l).name);
        label(2,k) = cellstr(obj_chess(l).chess(i).background);
        label(3,k) = cellstr(obj_chess(l).chess(i).position);
        k = k+1;
    end
end
% Ordering of chess
minimum = order(3,1); % Inizializzazione
for l=1:size(order,2)-1
    for i=l:size(order,2)
        if(order(3,l)>order(3,i))
            tmp = order(:,l);
            tmp2 = label(:,l);
            order(:,l) = order(:,i);
            label(:,l) = label(:,i);
            order(:,i) = tmp;
            label(:,i) = tmp2;
        end
    end
end
% Find Reflection type using colors
for l=1:size(obj_chess,1)
    left = 0;
    right = 0;
    arr = [];
    for i = 1:size(obj_chess(l).chess,2)
        v(1) = sum(sum(obj_chess(l).chess(i).mask&maskLeft));
        v(2) = sum(sum(obj_chess(l).chess(i).mask&maskCenter));
        v(3) = sum(sum(obj_chess(l).chess(i).mask&maskRight));
        idx = find(v==max(v));
        if(idx==1)
           left = left +1; 
        end
        if(idx==3)
           right = right +1; 
        end
        arr = [arr,obj_chess(l).chess(i).position];
    end
    if(left == 2 && right == 1)
       idx = find(strcmp(arr,positions(1)));
       obj_chess(l).chess(idx(1)).type = types(2);
       obj_chess(l).chess(idx(2)).type = types(1);
       idx = find(strcmp(arr,positions(3)));
       obj_chess(l).chess(idx(1)).type = types(1);
    end
    if(right == 2 && left == 1)
       idx = find(strcmp(arr,positions(3)));
       obj_chess(l).chess(idx(1)).type = types(1);
       obj_chess(l).chess(idx(2)).type = types(2);  
       idx = find(strcmp(arr,positions(1)));
       obj_chess(l).chess(idx(1)).type = types(1);
    end
    if(left == 1 && right == 1)
        for i = 1:size(obj_chess(l).chess,2)
            if(~strcmp(obj_chess(l).chess(i).position,positions(2)))
                curr_color = name2code(obj_chess(l).name);
                curr_background = name2code(obj_chess(l).chess(i).background);
                idc = find(sum(checker_vector(1,:,:),3)==sum(curr_background));
                if(find(sum(checker_vector(2,idc(1):idc(end),:),3)==sum(curr_color(1,:))))
                    curr_color = curr_color(1,:);
                else
                	curr_color = curr_color(2,:);
                end
                idx = find(order(3,:) == obj_chess(l).chess(i).centroid(1));
                labs = find(strcmp(label(3,:),positions(2)));
                if(idx < labs(1) && idx~=1)
                    color = name2code(obj_chess(order(1,idx-1)).name);
                    background = name2code(obj_chess(order(1,idx-1)).chess(order(2,idx-1)).background);
                    idc = find(sum(checker_vector(1,:,:),3)==sum(background));
                    if(find(sum(checker_vector(2,idc(1):idc(end),:),3)==sum(color(1,:))))
                        color = color(1,:);
                    else
                        color = color(2,:);
                    end
                    check = checker_vector(2,:,:);
                    ii = find(check(:,:,1) == color(1) & check(:,:,2) == color(2) & check(:,:,3) == color(3));
                    jj = find(check(:,:,1) == curr_color(1) & check(:,:,2) == curr_color(2) & check(:,:,3) == curr_color(3));                   
                    if(jj>ii)
                         obj_chess(l).chess(i).type = types(2);
                    else
                         obj_chess(l).chess(i).type = types(1);
                    end                 
                elseif(idx > labs(end) && idx~=size(order,2))
                    color = name2code(obj_chess(order(1,idx+1)).name);
                    background = name2code(obj_chess(order(1,idx+1)).chess(order(2,idx-1)).background);
                    idc = find(sum(checker_vector(1,:,:),3)==sum(background));
                    if(find(sum(checker_vector(2,idc(1):idc(end),:),3)==sum(color(1,:))))
                        color = color(1,:);
                    else
                        color = color(2,:);
                    end
                    check = checker_vector(2,:,:);
                    ii = find(check(:,:,1) == color(1) & check(:,:,2) == color(2) & check(:,:,3) == color(3));
                    jj = find(check(:,:,1) == curr_color(1) & check(:,:,2) == curr_color(2) & check(:,:,3) == curr_color(3));
                    if(jj<ii)
                         obj_chess(l).chess(i).type = types(2);
                    else
                         obj_chess(l).chess(i).type = types(1);
                    end
                end
            end 
        end
    end
end

imshow(maskCenter|maskLeft|maskRight)
figure, imshow(enlarged);
%% Sposta punti
maskC  = false(size(fuse,1),size(fuse,2));
maskL1 = false(size(fuse,1),size(fuse,2));
maskL2 = false(size(fuse,1),size(fuse,2));
maskR1 = false(size(fuse,1),size(fuse,2));
maskR2 = false(size(fuse,1),size(fuse,2));

for i=1:size(obj_chess,1)
    for j=1:size(obj_chess(i).chess,2) %Inew = I.*repmat(M,[1,1,3]);
        if(strcmp(obj_chess(i).chess(j).position,'Left') && strcmp(obj_chess(i).chess(j).type,'Secondary'))
            maskL2 = maskL2 | obj_chess(i).chess(j).mask;
        end
        if(strcmp(obj_chess(i).chess(j).position,'Right') && strcmp(obj_chess(i).chess(j).type,'Primary'))
            maskR1 = maskR1 | obj_chess(i).chess(j).mask;
        end
        if(strcmp(obj_chess(i).chess(j).position,'Left') && strcmp(obj_chess(i).chess(j).type,'Primary'))
            maskL1 = maskL1 | obj_chess(i).chess(j).mask;
        end
        if(strcmp(obj_chess(i).chess(j).position,'Right') && strcmp(obj_chess(i).chess(j).type,'Secondary'))
            maskR2 = maskR2 | obj_chess(i).chess(j).mask;
        end
        if(strcmp(obj_chess(i).chess(j).position,'Center'))
            maskC = maskC | obj_chess(i).chess(j).mask;
        end
    end
end

for k=1:5
    switch(k)
        case 1
            mask = maskC;
        case 2
            mask = maskL1;
        case 3
            mask = maskL2;
        case 4
            mask = maskR1;
        case 5
            mask = maskR2;
    end
    idx = find(mask == 1);
    [idx,idy]=ind2sub(size(mask),idx);
    j = boundary(idx,idy,0.2); % Parametro
    mask = poly2mask(idy(j),idx(j), size(mask,1), size(mask,2)); 
    ImgMat = im2double(I); ImgMat = (ImgMat(:,:,1)+ImgMat(:,:,2)+ImgMat(:,:,3))./3; ImgMat = ImgMat.*mask;
    filtered = medfilt2(ImgMat,[10,10]); 
    yfilter = fspecial('sobel').*3;
	xfilter = yfilter';
    Jx = imfilter(filtered,xfilter);
    Jy = imfilter(filtered,yfilter);
    J = sqrt(Jx.*Jx + Jy.*Jy);
    bw = hysthresh(J, 0.1, 0.8);
    figure, imshow(imclose(J,strel('square',5)));
end

    edge_magnitude = edge_magnitude + edge(mask);
    edge_magnitude = bwareaopen(edge_magnitude>0.9,10);
    figure, [rj, cj, re, ce] = findendsjunctions(edge_magnitude, 1);

for k=1:5
    switch(k)
        case 1
            mask = maskC;
        case 2
            mask = maskL1;
        case 3
            mask = maskL2;
        case 4
            mask = maskR1;
        case 5
            mask = maskR2;
    end
    idx = find(mask == 1);
    [idx,idy]=ind2sub(size(mask),idx);
    j = boundary(idx,idy,0.2); % Parametro
    mask = poly2mask(idy(j),idx(j), size(mask,1), size(mask,2)); 
%     ImgMat = im2double(I); ImgMat(:,:,1) = im2double(ImgMat(:,:,1)).*mask;
%     ImgMat(:,:,2) = im2double(ImgMat(:,:,2)).*mask;
%     ImgMat(:,:,3) = im2double(ImgMat(:,:,3)).*mask;
%     noOfColors = 5;%Number of colors to be present in the output
%     s_img = size(ImgMat);
%     inImg = ImgMat;
% 
%     %K-Means
%     r = inImg(:,:,1);
%     g = inImg(:,:,2);
%     b = inImg(:,:,3);
%     inputImg = zeros((s_img(1) * s_img(2)), 3);
%     inputImg(:,1) = r(:);
%     inputImg(:,2) = g(:);
%     inputImg(:,3) = b(:);
%     inputImg = im2double(inputImg);
%     disp('K-Means Processing Started');
%     [idx, C] = kmeans(inputImg, noOfColors, 'EmptyAction', 'singleton');
%     clusters_color = reshape(C,1,noOfColors,3);  
%     figure, imshow(clusters_color);
%     
%     disp('K-Means Processing Completed');
%     palette = C;
%     
% 
%     %Color Mapping
%     idx = uint8(idx);
%     outImg = im2double(zeros(s_img(1),s_img(2),3));
%     temp = reshape(idx, [s_img(1) s_img(2)]);
%     for i = 1 : 1 : s_img(1)
%         for j = 1 : 1 : s_img(2)
%             outImg(i,j,:) = palette(temp(i,j),:);
%         end
%     end
%     imshow(outImg); 
%   
    masked = im2double(I); 
    masked(:,:,1) = masked(:,:,1).*mask;
    masked(:,:,2) = masked(:,:,2).*mask;
    masked(:,:,3) = masked(:,:,3).*mask;
    
    filtered1 = Kuwahara(masked(:,:,1),4*4+1);
    filtered1 = filtered1.*imopen(bwareaopen(filtered1>0.2,1000),strel('square',5));
    filtered2 = Kuwahara(masked(:,:,2),4*4+1);
    filtered2 = filtered2.*imopen(bwareaopen(filtered2>0.2,1000),strel('square',5));
    filtered3 = Kuwahara(masked(:,:,3),4*4+1);
    filtered3 = filtered3.*imopen(bwareaopen(filtered3>0.2,1000),strel('square',5));
    filtered = im2double(cat(3,filtered1,filtered2,filtered3));

    for i=1:size(filtered,1)
        for j=1:size(filtered,2)
            c = reshape(filtered(i,j,:),1,3);
            if(sum(c) == 0)
                continue;
            end
            colors = getColors(i,j,obj_chess,100);
            v = zeros(1,size(colors,2));
            for l=1:size(v,2)
                v(l) = norm(colors(:,l)-c');
            end
            idx = find(v==min(v));
            c = colors(:,max(idx));
            if(max(c)>0)
                c = c./max(c);
            end
            filtered(i,j,:) = reshape(c,1,1,3);
        end
    end
    
    imshow(filtered);
    [edge_magnitude, edge_orientation, Jx, Jy, Jxy] = coloredges(filtered);
    edge_magnitude = bwareaopen(edge_magnitude>0.4,10);
    edge_magnitude = filledgegaps(edge_magnitude,10)+edge(mask);
    [rj, cj, re, ce] = findendsjunctions(edge_magnitude, 1);
end


imshow(edge_magnitude);

Dx = (im2double(Jx)-im2double(Jy))>0.4;
Dx = bwareaopen(Dx,20);
Dx = filledgegapsel(Dx,3,20);
imshow(Dx)

[rj, cj, re, ce] = findendsjunctions(Dx, 1);

Dy = (im2double(Jy)-im2double(Jx))>0.4;
Dy = bwareaopen(Dy,20);

imshowpair(Dx,Dy,'falsecolor');



Dxj = filledgegaps(Dx,10);
Dxj = filledgegapsel(Dxj,3,50);
imshowpair(Dxj,Dx,'falsecolor');

Jx = imopen(Jx>0.4,strel('square',2));
Jx = bwareaopen(Jx>0.4,20);
Jx = medfilt2(Jx,[7,2]);

Jx = bwareaopen(Jx,150);
Jy = bwareaopen(Jy>0.4,100);
Jy = medfilt2(Jy,[2,7]);
Jy = filledgegapsel(Jy,25,3);
Jxy = Jx & Jy;
figure, imshowpair(imdilate(Jxy,strel('disk',5)),filtered,'falsecolor');
        





imshow(RGB)
axis off
title('RGB Segmented Image')

for i=1:size(obj_chess,1)
    for j=1:size(obj_chess(i).chess,2) %Inew = I.*repmat(M,[1,1,3]);
        mask = bwconvhull(obj_chess(i).chess(j).mask);
        mask = imdilate(mask,strel('square',25));
        filtered1 = Kuwahara(im2double(I(:,:,1)).*mask,4*4+1);
        filtered2 = Kuwahara(im2double(I(:,:,2)).*mask,4*4+1);
        filtered3 = Kuwahara(im2double(I(:,:,3)).*mask,4*5+1);
        filtered = cat(3,filtered1,filtered2,filtered3);
        [edge_magnitude, edge_orientation, Jx, Jy, Jxy] = coloredges(filtered);
        Jx = bwareaopen(Jx>0.4,100);
        Jx = medfilt2(Jx,[7,2]);
        Jx = filledgegapsel(Jx,1,25);
        Jx = bwareaopen(Jx,150);
        Jy = bwareaopen(Jy>0.4,100);
        Jy = medfilt2(Jy,[2,7]);
        Jy = filledgegapsel(Jy,25,3);
        Jxy = Jx & Jy;
        figure, imshowpair(imdilate(Jxy,strel('disk',5)),filtered,'falsecolor');
    end
end

%% Fissa il centro degli assi
% inv_BW_axis = uint8(obj_red.black_white|obj_yellow.black_white);
% col_BW_axis = uint8(obj_red.color_mask|obj_yellow.color_mask);
% BW_axis = 255*(inv_BW_axis-col_BW_axis);
% fuse_axis = obj_red.masked_rgb+obj_yellow.masked_rgb+repmat(BW_axis,[1 1 3]);
% for i=1:size(fuse_axis,1)
%     for k=2:size(fuse_axis,2)
%         prev = fuse_axis(i,k-1);
%         curr = fuse_axis(i,k);
%         %yellow 255 255 0
%         %white 255 255 255
%         %red 255 0 0
%         %black 0 0 0
%         if(prev(1) == 255 && prev(2) == 255 && prev(1) == 0)
%             %yellow
%         end
%     end
% end
% 
% filtered1 = Kuwahara(im2double(I(:,:,1)),4*4+1);
% filtered2 = Kuwahara(im2double(I(:,:,2)),4*4+1);
% filtered3 = Kuwahara(im2double(I(:,:,3)),4*4+1);
% filtered = cat(3,filtered1,filtered2,filtered3);
% filter = fspecial('prewit');%[-ones(5,1),zeros(5,1),ones(5,1)];
% [edge_magnitude, edge_orientation, Jx, Jy, Jxy] = coloredges(filtered);
% Jx = bwareaopen(Jx>0.4,100);
% Jx = medfilt2(Jx,[7,2]);
% Jx = filledgegapsel(Jx,1,25);
% Jx = bwareaopen(Jx,150);
% Jy = bwareaopen(Jy>0.4,100);
% Jy = medfilt2(Jy,[2,7]);
% Jy = filledgegapsel(Jy,25,3);
% Jxy = Jx & Jy;
% imshowpair(imdilate(Jxy,strel('disk',5)),filtered,'falsecolor');
% 
% 
% 
% Jy = filledgegaps(Jy,11);
% imshowpair(Jx,Jy,'falsecolor');
% 
% f1 = imfilter(im2double(filtered(:,:,1)),filter)+imfilter(im2double(filtered(:,:,2)),filter)+imfilter(im2double(filtered(:,:,3)),filter);
% f2 = imfilter(im2double(filtered(:,:,1)),filter')+imfilter(im2double(filtered(:,:,2)),filter')+imfilter(im2double(filtered(:,:,3)),filter');
% 
% f1 = bwareaopen(imfilter(im2double(filtered(:,:,2)),filter)>0.5,30);
% f1 = imdilate(f1,strel('line',20,0));
% f2 = bwareaopen(imfilter(im2double(filtered(:,:,2)),filter')>0.5,30);
% f2 = imdilate(f2,strel('line',20,90));
% f = cat(3,f2,f1,zeros(size(f1)));
% f1 = bwareaopen(imfilter(im2double(filtered(:,:,1)),filter)>0.5,30);
% f1 = imdilate(f1,strel('line',20,0));
% f2 = bwareaopen(imfilter(im2double(filtered(:,:,1)),filter')>0.5,30);
% f2 = imdilate(f2,strel('line',20,90));
% f = f+cat(3,f2,f1,zeros(size(f1)));
% f1 = bwareaopen(imfilter(im2double(filtered(:,:,3)),filter)>0.5,30);
% f1 = imdilate(f1,strel('line',20,0));
% f2 = bwareaopen(imfilter(im2double(filtered(:,:,3)),filter')>0.5,30);
% f2 = imdilate(f2,strel('line',20,90));
% f = f+cat(3,f2,f1,zeros(size(f1)));
% imshowpair(f,I,'falsecolor')
% 
% 
% 
% [imagePoints,boardSize] = detectCheckerboardPoints(filtered);
% scatter(imagePoints(:,1),imagePoints(:,2));
% L = watershed(filtered(:,:,1));
% imshow(label2rgb(L))
% 
% [featureVector,hogVisualization] = extractHOGFeatures(filtered,'CellSize',[5 5]);
% figure;
% imshow(filtered);
% hold on;
% plot(hogVisualization);
% imshow(filtered);
% 
% [cim, r, c] = noble(filtered(:,:,1),1);
%  hold on; plot(r,c);
% 
% 
% 
% B = cartoon(im2double(I));
% imshow(B);
% 
% [l, Am, C] = slic(cleanimage, 10000, 20, 1, 'median');
% show(drawregionboundaries(l, IM1, [255 255 255]));
% lc = spdbscan(l, C, Am, 4);
% show(drawregionboundaries(lc, cleanimage, [255 255 255]))



% BW = uint8(255*(inv_BW_red-BW_red));
% Red = RedImage+repmat(BW,[1 1 3]);
% BW = uint8(255*(inv_BW_yellow-BW_yellow));
% Yellow = YellowImage+repmat(BW,[1 1 3]);
% ry = Red+Yellow;
% 
% img = (circshift(Yellow,1)-circshift(Yellow,-1))+(circshift(Red,1)-circshift(Red,-1));
 
% imshowpair(Red,Yellow,'Montage');
% col = 15;
% for i=1:size(fuse,2)-col
%     imshow(fuse(:,i:col+i,:));
% end
 
% figure, imshow(fuse);
% [c r] = getpts(1);
% loc = int32([c r]);
% if size(loc)>1
%     loc = [loc(1,1) loc(1,2)];
% end

% I = im2double(I);
% r1 = (I(:,:,1)-I(:,:,2))./(I(:,:,1)+I(:,:,2)); %Red
% r5 = (I(:,:,3)-I(:,:,1))./(I(:,:,3)+I(:,:,1)); %Blue
% r2 = (I(:,:,1)-I(:,:,3))./(I(:,:,1)+I(:,:,3)); %-Red & Yellow
% r7 = r2-r1;                                    %Yellow
% r4 = (I(:,:,2)-I(:,:,1))./(I(:,:,2)+I(:,:,1)); %-Green & Blue
% r8 = r4-r5;                                    %Blue
% r3 = (I(:,:,2)-I(:,:,3))./(I(:,:,2)+I(:,:,3)); %Yellow
% r6 = (I(:,:,3)-I(:,:,2))./(I(:,:,3)+I(:,:,2));
 
% figure, imshow(r7)
 
% figure, imshow(I);
% Red = I(:,:,1);
% Green = I(:,:,2);
% Blue = I(:,:,3);
% corners = detectFASTFeatures(Red);
% imshow(I); hold on;
% plot(corners);
% corners = detectFASTFeatures(Green);
% plot(corners);
% corners = detectFASTFeatures(Blue);
% plot(corners);
 
% s = 30;
% Red_corners = corner(Red,'Harris');
% Green_corners = corner(Green,'Harris');
% Blue_corners = corner(Blue,'Harris');
% figure,imshow(I), hold on, scatter(Red_corners(:,1),Red_corners(:,2));
% scatter(Green_corners(:,1),Green_corners(:,2));
% scatter(Blue_corners(:,1),Blue_corners(:,2));
 
% Gray = rgb2gray(I);
% Gray_corners = corner(Gray,'Harris');
% scatter(Gray_corners(:,1),Gray_corners(:,2));
 
% rgbImage = imread('_DSC0072.JPG');
% % Extract individual color channels and cast to double to avoid clipping.
% R = double(rgbImage(:,:,1));
% G = double(rgbImage(:,:,2));
% B = double(rgbImage(:,:,3));
% blueRatio = uint8(((100 * B)./(1+R+G)) .* (256./(1+B+R+G)));
% % Display images
% fontSize = 28;
% subplot(2, 3, 1);
% imshow(rgbImage);
% title('Original RGB Image', 'FontSize', fontSize);
% subplot(2, 3, 2);
% imshow(uint8(R)); % Must be uint8 for display
% title('R Image', 'FontSize', fontSize);
% subplot(2, 3, 3);
% imshow(uint8(G)); % Must be uint8 for display
% title('G Image', 'FontSize', fontSize);
% subplot(2, 3, 4);
% imshow(uint8(B)); % Must be uint8 for display
% title('B Image', 'FontSize', fontSize);
% subplot(2, 3, 5);
% imshow(blueRatio);
% title('Blue Ratio Image', 'FontSize', fontSize);
% % Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % Give a name to the title bar.
% set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off') 


