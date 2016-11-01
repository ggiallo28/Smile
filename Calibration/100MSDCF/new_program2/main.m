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
K(:,:,1) = imadjust(R.*(2-L));
K(:,:,2) = imadjust(G.*(2-L));
K(:,:,3) = imadjust(B.*(2-L));
I = im2uint8(K);
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
label = cell(4,k-1); k =1;
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
for j=1:size(order,2)
    idx_chess_vector = order(1,j);
    idx_color_chess = order(2,j);
    label(4,j) = obj_chess(idx_chess_vector).chess(idx_color_chess).type;
end
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

%% Full Mask
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
left_center_axis = center_axis.*(maskL1I|maskL2I);
idx = find(left_center_axis == 1);
[idy,idx] = ind2sub(size(RIGHTMASK),idx);
[left_fitresult, left_gof] = createLine(idx, idy);

right_center_axis = center_axis.*(maskR1I|maskR2I);
idx = find(right_center_axis == 1);
[idy,idx] = ind2sub(size(RIGHTMASK),idx);
[right_fitresult, right_gof] = createLine(idx, idy);

mid_center_axis = center_axis.*(maskCI);
idx = find(mid_center_axis == 1);
[idy,idx] = ind2sub(size(RIGHTMASK),idx);
[mid_fitresult, mid_gof] = createLine(idx, idy);

imshow(I); hold on; 
plot(left_fitresult,'b'); 
plot(right_fitresult,'r'); 
plot(mid_fitresult, 'y');
legend('left axis', 'right axis', 'center axis');
%% Corner
cymk = rgb2cmyk(I);
Cyano = im2double(cymk(:,:,1));
Magenta = im2double(cymk(:,:,2));
Yellow = im2double(cymk(:,:,3));
TheColors = Cyano+Magenta+Yellow;
Key = im2double(cymk(:,:,4));

KEY = (1-Key);
ORIG = rgb2gray(im2double(O));
[key_filtered, key_filtered_x, key_filtered_y] = absimfilter(KEY,fspecial('prewit'));
[orig_filtered, orig_filtered_x, orig_filtered_y] = absimfilter(ORIG,fspecial('prewit'));
% 
% [eeee, eeeee, eeeeee] = absimfilter(L,fspecial('sobel'));
% bw = im2bw(eeee,0.3*graythresh(eeee.*FullMask));
% bw = imdilate(bw,strel('square',2));
% bw = imerode(bw,strel('square',5));

img_filtered = 0.5*(key_filtered+orig_filtered);

% ORIG_Kywahara = Kuwahara(ORIG,(4*1+1));
% ORIG_edge = edge(ORIG_Kywahara,'canny').*FullMask;
% ORIG_edge = imdilate(ORIG_edge,strel('square',4));
% ORIG_edge = imopen(ORIG_edge,strel('square',4));
% ORIG_edge = bwareaopen(ORIG_edge, 4000); % Parametro
% ORIG_edge = imerode(ORIG_edge,strel('square',4));
% ORIG_edge = bwareaopen(ORIG_edge, 4000); % Parametro
% ORIG_edge = imclose(ORIG_edge,strel('square',4));
% img_filtered = ORIG_edge.*FullMask;
% 
% imshow(imfilter(1-Key,fspecial('unsharp')))

bw = im2bw(img_filtered,graythresh(img_filtered.*FullMask)); % Se si cambia la soglia qua non funge più un cazzo
bw = bwareaopen(bw, 100); % Parametro
grid_image = false(size(bw));
horizontalCell = cell(5,20);
verticalCell = cell(5,20);
pointsArray = struct();
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
    mask = imdilate(mask,strel('rectangle',[2 12]));
    mask = imerode(mask,strel('rectangle',[14 2]));
    bw_edge = filledgegaps(bw.*mask,1);
    CC = bwconncomp(~bw_edge,4);
    for i=1:size(CC.PixelIdxList,2)
        if(size(CC.PixelIdxList{i},1)<100)
            bw_edge(CC.PixelIdxList{i}) = 1;
        end
    end
    bw_edge = filledgegaps(bw_edge,1);     
    horizontalEdgeImage = imfilter(bw_edge, [-1 0 1]');
    horizontalEdgeImage = bwareaopen(horizontalEdgeImage, 10); % Parametro
    horizontalEdgeImage = filledgegaps(horizontalEdgeImage,9);
    horizontalEdgeImage = bwareaopen(horizontalEdgeImage, 100); % Parametro
    CCh = bwconncomp(horizontalEdgeImage,8);
    horizontal_image = false(size(horizontalEdgeImage));
    xxxd = [];
    yyyd = [];
    for i = 1:size(CCh.PixelIdxList,2)
        [idy,idx] = ind2sub(size(I),CCh.PixelIdxList{i});
        xxxd = [xxxd;-1;idx];
        yyyd = [yyyd;-1;idy];
        [horizonta_fitresult, horiz_gof] = createCircle(idx, idy);
        coeffs = coeffvalues(horizonta_fitresult);
        x = 1:size(horizontalEdgeImage,2);
        y = round(polyval(coeffs,x));
        x(y<1 | y>size(horizontalEdgeImage,1)) = [];
        y(y<1 | y>size(horizontalEdgeImage,1)) = [];
        for j = 1:size(x,2)
            horizontal_image(y(j),x(j)) = 1;
        end
        plot(x,y,'k'); 
        horizontalCell{k,i} = horizonta_fitresult;
    end
    grid_image = grid_image | horizontal_image;
    tmp_grid_image = horizontal_image;
    horizontal_image = imdilate(horizontal_image,strel('disk',2));
    verticalEdgeImage = imfilter(bw_edge, [-1 0 1]);
    verticalEdgeImage = bwareaopen(verticalEdgeImage, 20); % Parametro
    verticalEdgeImage = (verticalEdgeImage - horizontal_image)>0;
    verticalEdgeImage = filledgegapsel(verticalEdgeImage,3,8);
    verticalEdgeImage = bwareaopen(verticalEdgeImage, 100); % Parametro
    CCv = bwconncomp(verticalEdgeImage,8);
    vertical_image = false(size(verticalEdgeImage));
    for i = 1:size(CCv.PixelIdxList,2)
        [idy,idx] = ind2sub(size(I),CCv.PixelIdxList{i});
        [vertical_fitresult, vert_gof] = createLine(idy, idx);
        coeffs = coeffvalues(vertical_fitresult);
        y = 1:0.01:size(verticalEdgeImage,1);
        x = polyval(coeffs,y);
        for j = 1:size(x,2)
            vertical_image(round(y(j)),round(x(j))) = 1;
        end
        plot(x,y,'k');
        verticalCell{k,i} = vertical_fitresult;
    end
    grid_image = grid_image | vertical_image;   
%     Se la distanza tra due rette e minore di una soglia allora sono la
%     stessa retta, riusa i punti per farne una nuova, oppure se si
%     intersecano sull'immagine e sicuro un errore
%     [idy,idx] = ind2sub(size(I),CCo.PixelIdxList{1}); [idy1,idx1] = ind2sub(size(I),CCo.PixelIdxList{2});
%     [oriz_fitresult, oriz_gof] = createFit([idx;idx1], [idy;idy1]);
%     imshow(I); hold on; plot(oriz_fitresult,'b'); 
%     pause 
% usare prima le rette per unire pezzi spuri di contorni e poi interpolare
     mask = imdilate(mask,strel('rectangle',[12 12]));
     tmp_grid_image = (vertical_image | tmp_grid_image).*mask;
     [y_rj_matrix, x_cj_matrix] = findIntersections(horizontalCell(k,:), verticalCell(k,:), mask);
%     [rj, cj, ~, ~] = findendsjunctions(bw_edge.*mask, 0,k); % Eventualmente fare la media con questi per rendere il tutto più "preciso"
     scatter(x_cj_matrix(:),y_rj_matrix(:));
     pointsArray.x_points{k} = x_cj_matrix;
     pointsArray.y_points{k} = y_rj_matrix;
end
% positions = [{'Left'} {'Center'} {'Right'}];
% types = [{'Primary'} {'Secondary'} {'Real'}];
%% Divido i punti per ogni componente
% Scansioniamo le tessere, facciamo intersezione con le maschere di ogni
% singola tessera, i punti che ricadono dell'intersezione vengono assegnati
% alla tessera associata alla specifica maschera.
% Per fare questo usiamo il vettore delle label e order creato
% precedentemente per poter facilemente recuperare gli indici di riga e
% colonna associati ad ogni singolo elemento.
for k=1:5
    switch(k)
        case 1
            mask = maskCI;
            idx = find(strcmp(label(3,:),positions(2)));
        case 2
            mask = maskL1I;
            idx = find(strcmp(label(3,:),positions(1)) & strcmp(label(4,:),types(1)));
        case 3
            mask = maskL2I;
            idx = find(strcmp(label(3,:),positions(1)) & strcmp(label(4,:),types(2)));
        case 4
            mask = maskR1I;
            idx = find(strcmp(label(3,:),positions(3)) & strcmp(label(4,:),types(1)));
        case 5
            mask = maskR2I;
            idx = find(strcmp(label(3,:),positions(3)) & strcmp(label(4,:),types(2)));
    end
    for i=1:size(idx,2)
        idx_chess_vector = order(1,idx(i));
        idx_color_chess = order(2,idx(i));
        x_cj_matrix = pointsArray.x_points{k};
        y_rj_matrix = pointsArray.y_points{k};
        ch_mask = obj_chess(idx_chess_vector).chess(idx_color_chess).ch_mask;
        w = find(sum(ch_mask) ~=0);
        w_min_x = min(w); w_max_y = max(w);
        check_matrix = false(size(x_cj_matrix));
        for j=1:size(x_cj_matrix,1)
            for q=1:size(x_cj_matrix,2)
                if(x_cj_matrix(j,q)>=(w_min_x-10) && x_cj_matrix(j,q)<=(w_max_y+10)) % 10 pixel di tolleranza: Parametro, speriamo vada bene
                    check_matrix(j,q) = true;
                end
            end
        end
        num_cols = mode(sum(check_matrix,2));
        sum_cols = sum(check_matrix);
        check_matrix = false(size(check_matrix)); % Fix, dobbiamo sempre prendere tutta una colonna, scegliamo il valore più frequente e completiamo
        for j=1:num_cols
            id = find(sum_cols ==max(sum_cols), 1 ); % prendi solo il primo
            check_matrix(:,id) = 1;
            sum_cols(id) = 0;
        end    
        y_rj_matrix_copy = y_rj_matrix.*check_matrix;
        x_cj_matrix_copy = x_cj_matrix.*check_matrix;
        id = find(y_rj_matrix_copy(1,:) == 0);
        y_rj_matrix_copy(:,id) = [];
        x_cj_matrix_copy(:,id) = [];
        
        obj_chess(idx_chess_vector).chess(idx_color_chess).intersections_x = x_cj_matrix_copy;
        obj_chess(idx_chess_vector).chess(idx_color_chess).intersections_y = y_rj_matrix_copy;   
        imshow([I;repmat(255.*ch_mask,1,1,3)]); hold on;
        scatter(x_cj_matrix_copy(:),y_rj_matrix_copy(:));
    end    
end
%% Padding delle matrici e allineamento
% Per ogni combinazione colore/sfondo controlliamo quant'e il massimo
% numero di colorre, aggiungiamo zeri alle matrici che hanno un numero
% inferiore di colonne
for l=1:size(obj_chess,1)
    v = zeros(1,size(obj_chess(l).chess,2));
    for q=1:size(obj_chess(l).chess,2)
       v(q) = size(obj_chess(l).chess(q).intersections_x,2); 
    end
    max_v = max(v);
    for q=1:size(obj_chess(l).chess,2)
        num_col = size(obj_chess(l).chess(q).intersections_x,2);
        num_row = size(obj_chess(l).chess(q).intersections_x,1);
        if(num_col<max_v)
            obj_chess(l).chess(q).intersections_x = ...
                [obj_chess(l).chess(q).intersections_x, zeros(num_row,max_v-num_col)];
            obj_chess(l).chess(q).intersections_y = ...
                [obj_chess(l).chess(q).intersections_y, zeros(num_row,max_v-num_col)];
        end
        % Allineo a sinistra o a destra in funzione del fatto che ho un elemento a sinistra o a destra.
        % Ciò mi è utile per tenere tutti i punti compatti verso il centro
        % di ciò che è visibile sulla checherboard.
        % Per capire se ho un elemento a sinistra controllo che a sinistra
        % ci sia una griglia dello stesso tipo e nella stessa posizione
        % della corrente.
        % Per capire se ho un elemento a destra faccio la medesima cosa
        % guardando a destra.
        % Se abbiamo una griglia sinistra e una griglia destra sto vedendo
        % tutta la griglia corrente, quindi non ha importanza in che
        % direzione shifto.
        % Una volta che gli elementi sono stati allineati in questo modo è
        % sufficiente flippare le matrici associate ai riflessi primari,
        % unico caso in cui abbiamo un'inversione.
        % Gli indici di riga/colonna rappresentano le associazioni, abbiamo
        % degli zeri dove l'associazione non esiste poichè il punto non
        % risulta visibile.
        % LABEL: Colore Foreground, Colore Background, Posizione, Tipo
        % ORDER: Indice nel vettore obj_chess, Indice Chess, Centroide
        idl = find(order(1,:)==l);
        idq = find(order(2,:)==q);
        idx_order = intersect(idl,idq);
        if(hasLeft(label,idx_order))
          obj_chess(l).chess(q).intersections_x = ...
              allignleftdouble(obj_chess(l).chess(q).intersections_x);
          obj_chess(l).chess(q).intersections_y = ...
              allignleftdouble(obj_chess(l).chess(q).intersections_y);
        end
        if(hasRight(label,idx_order))
           obj_chess(l).chess(q).intersections_x = ...
               allignrightdouble(obj_chess(l).chess(q).intersections_x);
           obj_chess(l).chess(q).intersections_y = ...
               allignrightdouble(obj_chess(l).chess(q).intersections_y);
        end
        % positions = [{'Left'} {'Center'} {'Right'}];
        % types = [{'Primary'} {'Secondary'} {'Real'}];
        if(strcmp(label(4,idx_order),types(1)))
           obj_chess(l).chess(q).intersections_x = ...
               fliplr(obj_chess(l).chess(q).intersections_x);
           obj_chess(l).chess(q).intersections_y = ...
               fliplr(obj_chess(l).chess(q).intersections_y);
        end
    end
end

%% Mapping sulla superficie del Manifold
axis_line = reshape([192,192,192;192,192,192],2,1,3);
checher_vector_with_axis = uint8([checker_vector(:,1:3,:),axis_line,checker_vector(:,4:6,:)]); 
figure, imshow(checher_vector_with_axis);

% La checker è suddivisa verticalmente in 4 quadrati, ciò significa 5 punti compresi quelli di intersezione
% Ogni quadrante rappresenta 45°, quindi ad ogni punto corrisponde un angolo di 9°
for l=1:size(obj_chess,1)
    for i=1:size(obj_chess(l).chess(1).intersections_x,1)
        for j=1:size(obj_chess(l).chess(1).intersections_x,2)
            imshow(I); hold on;
            vect_x = [];
            vect_y = [];
            for k =1:size(obj_chess(l).chess,2)
                axis_distance = axisdistance(obj_chess(l).name,obj_chess(l).chess(k).background,checher_vector_with_axis);
                adjusted_distance = axis_distance-1*sign(axis_distance); % Se positivo togliamo 1, se negativo aggiungiamo 1, le tessere gialle e rosse stanno a distanza zero.
                if(sign(axis_distance) >0)
                    adjusted_offset = (j-1)*11.25;
                else
                    adjusted_offset = (j-size(obj_chess(l).chess(k).intersections_x,2))*11.25;
                end
                angle = adjusted_distance*45 + adjusted_offset
                h = i*2.6
                vect_x = [vect_x, obj_chess(l).chess(k).intersections_x(i,j)];
                vect_y = [vect_y, obj_chess(l).chess(k).intersections_y(i,j)];
            end
            scatter(vect_x(:),vect_y(:)); hold off
            pause
        end
    end
end

%% Ora come associo i punti se levo il cilindro? Con il codice che già hanno ma usando la griglia.
%% TODO: Il prof ha detto di stimare la differenza di colore negli histogrammi del rosso al centro e del rosso a destra

%% Come stimo l'angolo degli specchi?