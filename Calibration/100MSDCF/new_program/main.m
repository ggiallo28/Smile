close all; clear all; clc;
figure,imshow(imread('checkerboard.jpg'));
checker_vector = reshape([[0,0,0;255,0,255];[0,0,0;0,255,255];[0,0,0;255,255,0];[255,255,255;255,0,0];[255,255,255;0,255,0];[255,255,255;0,0,255]],[2,6,3]);
checker_center = [0.5*size(checker_vector,2),0.5*size(checker_vector,2)+1];
orig = imread(['foto/DSC0012',num2str(3),'.JPG']);
%% Normalizzazione
figure, O = imcrop(orig);
R = im2double(O(:,:,1));
G = im2double(O(:,:,2));
B = im2double(O(:,:,3));
L = (0.2126 * R) + (0.7152 * G) + (0.0722 * B); % 
K = zeros(size(O));
K(:,:,1) = R.*(2-L);
K(:,:,2) = G.*(2-L);
K(:,:,3) = B.*(2-L);
M = mean(max(max(K(:,:,1),max(K(:,:,2),K(:,:,3)))));
T = 0.3*mean([graythresh(K(:,:,1)),graythresh(K(:,:,2)),graythresh(K(:,:,3))]);
R1 = K(:,:,1); R1(R<T) = 0;
G1 = K(:,:,2); G1(G<T) = 0;
B1 = K(:,:,3); B1(B<T) = 0;
I = im2uint8(K);
I(repmat(R1==0 & G1 == 0 & B1 == 0,1,3)) = 0;
%% Segmentazione
obj_red = createMask(I,'red');
obj_green = createMask(I,'green');
obj_blue = createMask(I,'blue');
obj_yellow = createMask(I,'yellow');
obj_chess = [obj_red;obj_green;obj_blue;obj_yellow];
%% Plot Results
inv_BW = uint8(obj_red.black_white|obj_green.black_white|obj_blue.black_white|obj_yellow.black_white);
col_BW = uint8(obj_red.color_mask|obj_green.color_mask|obj_blue.color_mask|obj_yellow.color_mask);
dual_BW = uint8(obj_red.inv_color_mask|obj_green.inv_color_mask|obj_blue.inv_color_mask|obj_yellow.inv_color_mask);
BW = 255*(inv_BW-col_BW);
colors_fuse = obj_red.masked_rgb+obj_green.masked_rgb+obj_blue.masked_rgb+obj_yellow.masked_rgb;
fuse = colors_fuse+repmat(BW,[1 1 3]);
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
% TODO ANLCUNI PUNTI VENGONO MESSI MALE
%% Separa i riflessi
enlarged = imdilate(fuse,strel('disk',20));
CC = bwconncomp(rgb2gray(enlarged)~=0);
assert(CC.NumObjects == 3, 'Ci sono troppi riflessi');
maskCenter = false(size(rgb2gray(enlarged))); maskCenter(CC.PixelIdxList{2}) = 1;
maskLeft = false(size(rgb2gray(enlarged))); maskLeft(CC.PixelIdxList{1}) = 1;
maskRight = false(size(rgb2gray(enlarged))); maskRight(CC.PixelIdxList{3}) = 1;
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
    %% TODO Sistemare sta cosa non va bene fare assunzioni sul numero di riflessi che vedi
    if(left >= 2 && right == 1)
       idx = find(strcmp(arr,positions(1)));
       obj_chess(l).chess(idx(1)).type = types(2);
       obj_chess(l).chess(idx(2)).type = types(1);
       idx = find(strcmp(arr,positions(3)));
       obj_chess(l).chess(idx(1)).type = types(1);
    end
    if(right >= 2 && left == 1)
       idx = find(strcmp(arr,positions(3)));
       obj_chess(l).chess(idx(1)).type = types(1);
       obj_chess(l).chess(idx(2)).type = types(2);  
       idx = find(strcmp(arr,positions(1)));
       obj_chess(l).chess(idx(1)).type = types(1);
    end
    empty_type = false;
    for i = 1:size(obj_chess(l).chess,2)
        if(isempty(obj_chess(l).chess(i).type))
            empty_type = true;
        end
    end
    if(empty_type)
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
%% Seprazione Riflessi
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
    j = boundary(idx,idy,0.1); % Parametro
    maskI = poly2mask(idy(j),idx(j), size(mask,1), size(mask,2));
    DilatedMask = imdilate(maskI,strel('rectangle',[3,20]));
    idx = find(maskI == 1);
    [idx,idy]=ind2sub(size(mask),idx);
    maskI = false(size(maskI));
    maskI(min(idx):max(idx),min(idy):max(idy)) = DilatedMask(min(idx):max(idx),min(idy):max(idy));
    switch(k)
        case 1
            maskCI = maskI;
        case 2
            maskL1I = maskI;
        case 3
            maskL2I = maskI;
        case 4
            maskR1I = maskI;
        case 5
            maskR2I = maskI;
    end   
end
%% Calcolo Convexhull maschere
for i=1:size(obj_chess,1)
    for j=1:size(obj_chess(i).chess,2) %Inew = I.*repmat(M,[1,1,3]);
        cut_x = obj_chess(i).bbox_x(j,:);
        cut_y = obj_chess(i).bbox_y(j,:); 
        if(strcmp(obj_chess(i).chess(j).position,'Left') && strcmp(obj_chess(i).chess(j).type,'Secondary'))
            mask = maskL2I;
        end
        if(strcmp(obj_chess(i).chess(j).position,'Right') && strcmp(obj_chess(i).chess(j).type,'Primary'))
            mask = maskR1I;
        end
        if(strcmp(obj_chess(i).chess(j).position,'Left') && strcmp(obj_chess(i).chess(j).type,'Primary'))
            mask = maskL1I;
        end
        if(strcmp(obj_chess(i).chess(j).position,'Right') && strcmp(obj_chess(i).chess(j).type,'Secondary'))
            mask = maskR2I;
        end
        if(strcmp(obj_chess(i).chess(j).position,'Center'))
            mask = maskCI;
        end
        obj_chess(i).chess(j).ch_mask = false(size(fuse,1),size(fuse,2));      
        obj_chess(i).chess(j).ch_mask(cut_y(1):cut_y(2),cut_x(1):cut_x(2)) = mask(cut_y(1):cut_y(2),cut_x(1):cut_x(2));
    end
end

% for k=1:5
%     switch(k)
%         case 1
%             mask = maskC;
%         case 2
%             mask = maskL1;
%         case 3
%             mask = maskL2;
%         case 4
%             mask = maskR1;
%         case 5
%             mask = maskR2;
%     end
%     idx = find(mask == 1);
%     [idx,idy]=ind2sub(size(mask),idx);
%     j = boundary(idx,idy,0.2); % Parametro
%     mask = poly2mask(idy(j),idx(j), size(mask,1), size(mask,2));
%     white = imopen(logical(BW),strel('square',5));
%     black = imopen(dual_BW & ~col_BW & ~inv_BW,strel('square',5));
%     black2 = imdilate(black,strel('square',15));
%     black2 = bwareaopen(black2, 4000); % Parametro
%     black = black & black2;
%     masked = im2double(I);
%     masked(repmat(black,1,1,3)) = 0;
%     masked(repmat(logical(col_BW),1,1,3)) = 0;  
%     masked(repmat(white,1,1,3)) = 255;
%     masked = (masked + im2double(colors_fuse)).*repmat(mask,1,1,3);    
%     for i=1:size(masked,1)
%         for j=1:size(masked,2)
%             c = masked(i,j,:);        
%             if(sum(c) == 0 || sum(c) == 3)
%                 continue;
%             end
% %             isWhiteORblack = (abs(c(:,:,1)-c(:,:,2)) + abs(c(:,:,1)-c(:,:,3)) + abs(c(:,:,2)-c(:,:,3))) < 0.4;
% %             if(isWhiteORblack)
% %                 white = reshape([1 1 1],1,1,3); 
% %                 black = reshape([0 0 0],1,1,3); 
% %                 colors = [white, black];
% %             else
%                 colors = getColors(i,j,obj_chess,10,c);
% %             end           
%             v = zeros(1,size(colors,2));
%             for l=1:size(v,2)
%                 dist = colors(:,l,:)-c;
%                 dist = [dist(1,1,2) dist(1,1,1) dist(1,1,3)];
%                 v(l) = norm(dist);
%             end
%             idx = find(v==min(v));
%             c = colors(:,max(idx),:);
%             if(max(c)>0)
%                 c = c./max(c);
%             end
%             masked(i,j,:) = c;
%         end
%     end
%     
%     masked(:,:,1) = imclose(bwareaopen(masked(:,:,1),300),strel('square',10));
%     masked(:,:,2) = imclose(bwareaopen(masked(:,:,2),300),strel('square',10));
%     masked(:,:,3) = imclose(bwareaopen(masked(:,:,3),300),strel('square',10));
%     masked(:,:,1) = imopen(bwareaopen(masked(:,:,1),300),strel('square',10));
%     masked(:,:,2) = imopen(bwareaopen(masked(:,:,2),300),strel('square',10));
%     masked(:,:,3) = imopen(bwareaopen(masked(:,:,3),300),strel('square',10));
%     [edge_magnitude, edge_orientation, Jx, Jy, Jxy] = coloredges(masked);
%     edge_magnitude = bwareaopen(edge_magnitude>0.4,20);
%     edge_magnitude = filledgegaps(edge_magnitude,10)+edge(mask); % Capire se lasciare
%     [rj, cj, re, ce] = findendsjunctions(edge_magnitude, 1);
%     figure, imshow(masked);
%     pause
% end
%% Full Mask
cymk = rgb2cmyk(I);
Cyano = im2double(cymk(:,:,1));
Magenta = im2double(cymk(:,:,2));
Yellow = im2double(cymk(:,:,3));
TheColors = Cyano+Magenta+Yellow;
Key = im2double(cymk(:,:,4));
FullMask = maskCI+maskL1I+maskL2I+maskR1I+maskR2I;
BlurredMask = imgaussfilt(FullMask,10);
%% Fissa il centro degli assi
MINRGB = min(R,G); MINRGB = min(MINRGB,B);
MINRGB = imgaussfilt(MINRGB.*BlurredMask,2);
T = graythresh(MINRGB);
bw = im2bw(MINRGB,T);
RIGHTMASK = bwareaopen(bw, 300); % Parametro
for i=1:3
    switch(i)
        case 1
            mask = maskCI;
        case 2
            mask = maskL1I|maskL2I;
        case 3
            mask = maskR1I|maskR2I;
    end
    idx = find(RIGHTMASK.*mask == 1);
    [idx,idy]=ind2sub(size(RIGHTMASK),idx);
    j = boundary(idx,idy,0.1); % Parametro
    RIGHTMASK = RIGHTMASK | poly2mask(idy(j),idx(j), size(RIGHTMASK,1), size(RIGHTMASK,2));
end
LEFTMASK = FullMask-RIGHTMASK;
LEFTMASK = imclose(LEFTMASK,strel('square',5));
LEFTMASK = imopen(LEFTMASK,strel('square',5));
CC = bwconncomp(LEFTMASK,4); CC_dim = zeros(1,size(CC.PixelIdxList,2));
for i=1:size(CC.PixelIdxList,2)
    CC_dim(i) = size(CC.PixelIdxList{i},1);
end
LEFTMASK = false(size(LEFTMASK));
for i=1:3
    idx = find(CC_dim == max(CC_dim));
    LEFTMASK(CC.PixelIdxList{idx}) = 1;
    CC_dim(idx) = 0;
end
left = edge(LEFTMASK);
right = edge(RIGHTMASK);
center_axis = left & right;
imshow(center_axis);
%% Corner
cymk = rgb2cmyk(I);
Cyano = im2double(cymk(:,:,1));
Magenta = im2double(cymk(:,:,2));
Yellow = im2double(cymk(:,:,3));
TheColors = Cyano+Magenta+Yellow;
Key = im2double(cymk(:,:,4));

KEY = (1-Key).*FullMask;
ORIG = rgb2gray(im2double(O)).*FullMask;

x_filter = fspecial('prewit');
y_filter = x_filter';
filtered_x = 0.5*(absimfilter(KEY,x_filter) + absimfilter(ORIG,x_filter));
filtered_y = 0.5*(absimfilter(KEY,y_filter) + absimfilter(ORIG,y_filter));
img_filtered = sqrt(filtered_y.^2 + filtered_x.^2);
bw = im2bw(img_filtered,graythresh(img_filtered));
bw = bwareaopen(bw, 100); % Parametro
for k=1:5
    switch(k)
        case 1
            mask = maskCI;
        case 2
            mask = maskL1I;
        case 3
            mask = maskL2I;
        case 4
            mask = maskR1I;
        case 5
            mask = maskR2I;
    end
    bw_edge = filledgegaps(bw.*mask+edge(mask),1);
    CC = bwconncomp(~bw_edge,4);
    for i=1:size(CC.PixelIdxList,2)
        size(CC.PixelIdxList{i},1)
        if(size(CC.PixelIdxList{i},1)<40)
            bw_edge(CC.PixelIdxList{i}) = 1;
        end
    end
    bw_edge = filledgegaps(bw_edge,1);    
    [rj, cj, re, ce] = findendsjunctions(bw_edge, 1,k);   
    figure, imshow([bw_edge;ORIG.*mask]);
end



imshow(img_filtered);


[rj, cj, re, ce] = findendsjunctions(img_filtered, 1);
img_filtered(re,ce) = 0;





imshow([filtered_x>0.05;MINRGB]);
img_filtered = medfilt2(img_filtered,[10,10]);



T = adaptthresh(TheColors,1,'ForegroundPolarity','bright');
BW = imbinarize(TheColors,T);
BW = imopen(BW & FullMask,strel('square',15));
imshow(BW)





figure
imshow(BW)
[Gx, Gy] = imgradientxy(R,'central');
figure,
imshowpair(Gx,Gy,'falsecolor');





