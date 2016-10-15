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
figure, I = imcrop(orig);
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
cymk = rgb2cmyk(I);
R = im2double(I(:,:,1));
G = im2double(I(:,:,2));
B = im2double(I(:,:,3));
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
K = zeros(size(I));
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
order = [];
for l=1:size(obj_chess,1)
    for i = 1:size(obj_chess(l).chess,2)
        order = [order,[l;i;obj_chess(l).chess(i).centroid(1)]];
    end
end
min = order(3,1); % Inizializzazione
for l=1:size(order,2)-1
    for i=l:size(order,2)
        if(order(3,l)>order(3,i))
            tmp = order(:,l);
            order(:,l) = order(:,i);
            order(:,i) = tmp;
        end
    end
end
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
        flipped_cv = fliplr(checker_vector);
        for i = 1:size(obj_chess(l).chess,2)
            if(~strcmp(obj_chess(l).chess(i).position,positions(2)))
%                curr_color = obj_chess(l).name
%                curr_background = obj_chess(l).chess(i).background
                idx = find(order(3,:) == obj_chess(l).chess(i).centroid(1));
                if(idx<=2 || idx>=size(order,2)-1)
                    obj_chess(l).chess(i).type = types(2);
                else
                    obj_chess(l).chess(i).type = types(1);
                end
%                 if(idx>1 && idx<size(order,2))
%                     left_color = obj_chess(order(1,idx-1)).name;
%                     left_background = obj_chess(order(1,idx-1)).chess(order(2,idx-1)).background;
%                     right_color = obj_chess(order(1,idx+1)).name;
%                     right_background = obj_chess(order(1,idx+1)).chess(order(2,idx+1)).background;
%                 end
            end 
        end
    end
    
    
%     left_color = 
%     right_color
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


