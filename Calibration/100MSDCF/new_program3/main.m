close all; clear all; clc;
figure,imshow(imread('checkerboard.jpg'));
checker_vector = reshape([[0,0,0;255,0,255];[0,0,0;0,255,255];[0,0,0;255,255,0];[255,255,255;255,0,0];[255,255,255;0,255,0];[255,255,255;0,0,255]],[2,6,3]);
checker_center = [0.5*size(checker_vector,2),0.5*size(checker_vector,2)+1];
orig = imread(['foto/DSC0012',num2str(7),'.JPG']);
orig_bg = imread(['foto/DSC0012',num2str(8),'.JPG']);
%% Normalizzazione
figure, [O, BB] = imcrop(orig);
R = im2double(O(:,:,1));
G = im2double(O(:,:,2));
B = im2double(O(:,:,3));
L = (0.2126 * R) + (0.7152 * G) + (0.0722 * B); 
K = zeros(size(O));
K(:,:,1) = R.*(2-L);
K(:,:,2) = G.*(2-L);
K(:,:,3) = B.*(2-L);
I = im2uint8(K);
O_BG = imcrop(orig_bg,BB);
Rt = im2double(O_BG(:,:,1));
Gt = im2double(O_BG(:,:,2));
Bt = im2double(O_BG(:,:,3));
L = (0.2126 * Rt) + (0.7152 * Gt) + (0.0722 * Bt); 
K = zeros(size(O));
K(:,:,1) = Rt.*(2-L);
K(:,:,2) = Gt.*(2-L);
K(:,:,3) = Bt.*(2-L);
I_BG = im2uint8(K);
%% Parametri
confidence = 2.8;
[r,c] = size(I); Threshold = round(r*c/7000); % Soglia dimensione blob normalizzata alla dimensione dell'immagine
mpd = 15;
windowSize = 5;
%% Background Subtraction
lab_image_bg = rgb2lab(I_BG);
lab_imgae_fg = rgb2lab(I);
a_band = lab_image_bg(:,:,2);
mu_a = mean2(a_band);
si_a = std2(a_band); 
ll1 = a_band(:); ll1(ll1<=mu_a) = [];
ll2 = a_band(:); ll2(ll2>=mu_a) = [];
% Controlliamo quanti pixel ci sono a sinistra e a destra della media per aggiustare la soglia
p = size(ll1,1)/size(ll2,1);

channel2Min = mu_a-confidence*p*si_a;
channel2Max = mu_a+confidence*(2-p)*si_a;
BW1 = (lab_imgae_fg(:,:,2) >= channel2Min ) & (lab_imgae_fg(:,:,2) <= channel2Max);
maskedRGBImage = I; maskedRGBImage(repmat(bwareaopen(imopen(bwareaopen(~BW1,Threshold),strel('square',10)),round(Threshold/2)),[1 1 3])) = 0;
% Ripetiamo quello che abbiamo fatto prima sulla banda b
b_band = lab_image_bg(:,:,3); 
mu_b = mean2(b_band);
si_b = std2(b_band);
ll1 = b_band(:); ll1(ll1<=mu_b) = [];
ll2 = b_band(:); ll2(ll2>=mu_b) = [];
% Controlliamo quanti pixel ci sono a sinistra e a destra della media per aggiustare la soglia
p = size(ll1,1)/size(ll2,1);

channel2Min = mu_b-confidence*p*si_b;
channel2Max = mu_b+confidence*(2-p)*si_b;
BW2 = (lab_imgae_fg(:,:,3) >= channel2Min ) & (lab_imgae_fg(:,:,3) <= channel2Max);
maskedRGBImage(repmat(bwareaopen(imopen(bwareaopen(~BW2,Threshold),strel('square',10)),round(Threshold/2)),[1 1 3])) = 0;

result = I-maskedRGBImage;
bw = im2bw(result,0);
stats = regionprops(bw,'BoundingBox','Area','Image');
CC = bwconncomp(bw);
for i=1:size(stats)
    if(stats(i).BoundingBox(3) > 2*stats(i).BoundingBox(4))
        if(stats(i).Area < 10000)
            bw(CC.PixelIdxList{i}) = 0;
        end
    end
end
figure, imshow(maskedRGBImage);
figure, imshow((I-maskedRGBImage).*repmat(uint8(bw),[1 1 3]));
inImg = (I-maskedRGBImage).*repmat(uint8(bw),[1 1 3]);
%% Histogram Peak Finding
Rp = inImg(:,:,1); Rp = Rp(bw);
Gp = inImg(:,:,2); Gp = Gp(bw);
Bp = inImg(:,:,3); Bp = Bp(bw);
RGB = cat(3,Rp,Gp,Bp);
hsv = rgb2hsv(RGB);
[output_peak, output_minima_low, output_minima_high, output_minima_mid, hist_size] = findlocalminima(hsv(:,:,1),mpd,windowSize,0,1);
% Non sono in grado d distinguere tra il viola e il rosso, tra l'azzurro e il blu, quindi 4 cluster invece che 6
figure, imshow(imread('hsv.jpg'));
%% Segmentazione RGB, hist_size, output_minima_mid, Threshold, band
obj_red = createMask(I, bw, hist_size, output_minima_mid, Threshold, 'red');
obj_green = createMask(I, bw, hist_size, output_minima_mid, Threshold, 'green');
obj_blue = createMask(I, bw, hist_size, output_minima_mid, Threshold, 'blue');
obj_yellow = createMask(I, bw, hist_size, output_minima_mid, Threshold, 'yellow');
obj_chess = [obj_red;obj_green;obj_blue;obj_yellow];
% Se il contorno di sopra o sotto è troppo piccolo usare bbox
% Per eliminare merde laterali puoi calcolare la maskera con la fit square,
% erodere e risegmentare
% dilatare un pochino il bbox perchè nella fit square può succedere che si
% vada a finire fuori
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
        for j=1:size(obj.chess(i).center_x,2)
            for k=1:size(obj.chess(i).center_x,1)
                scatter(obj.chess(i).center_x(k,j),obj.chess(i).center_y(k,j)) % Riferimento Y verso il basso
            end
            plot(obj.chess(i).v_lines_centroid{j});
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
    % Se sono rimaste tessere senza nessuna labels associata
    if(empty_type)
        for i = 1:size(obj_chess(l).chess,2)
            % Se non si tratta di una tessera centrale, non ha senso analzizarla di sicuro non è ne primaria ne secondaria
            if(~strcmp(obj_chess(l).chess(i).position,positions(2)))
                % Recovery del colore, non siamo in grado di distinguere tra il viola/rosso blu/azzurro
                curr_color = name2code(obj_chess(l).name);
                curr_background = name2code(obj_chess(l).chess(i).background);
                idc = find(sum(checker_vector(1,:,:),3)==sum(curr_background));
                if(find(sum(checker_vector(2,idc(1):idc(end),:),3)==sum(curr_color(1,:))))
                    curr_color = curr_color(1,:);
                else
                	curr_color = curr_color(2,:);
                end
                % Cerchiamo la tessera corrente nel vettore di ordine
                idx = find(order(3,:) == obj_chess(l).chess(i).centroid(1));
                % Cerchiamo la posizione della tessera centrale
                labs = find(strcmp(label(3,:),positions(2)));
                % Se idx è precedente della prima tessera centrale e non è la prima tessera
                % Allora guardiamo a sinistra e preleviamo il colore della tessera che sta a sinsitra della corrente 
                
                % DOMANDA: Se è la prima tessera idx = 1 a non avere associazione? Questo non può capitare
                % devo vedere almeno due colori per riflesso, quindi della tessera con idx=1 devo necessariamnte vedere due rifelssi
                % a sinistra o a destra, idem per l'ultima. Il caso una a destra e uno a sinistra è gestito precedentemente.
                if(idx < labs(1) && idx~=1)
                    % Recovery del colore, non siamo in grado di distinguere tra il viola/rosso blu/azzurro
                    color = name2code(obj_chess(order(1,idx-1)).name);
                    background = name2code(obj_chess(order(1,idx-1)).chess(order(2,idx-1)).background);
                    idc = find(sum(checker_vector(1,:,:),3)==sum(background));
                    if(find(sum(checker_vector(2,idc(1):idc(end),:),3)==sum(color(1,:))))
                        color = color(1,:);
                    else
                        color = color(2,:);
                    end
                    % Sfruttiamo il checker vector, se i colori sono invertiti rispetto all'ordine originale allora
                    % abbiamo un riflesso primario, altrimenti abbiamo un riflesso secondario. 
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
%% Creazione Maschere per separazione riflessi
maskC  = false(size(fuse,1),size(fuse,2));
maskL1 = false(size(fuse,1),size(fuse,2));
maskL2 = false(size(fuse,1),size(fuse,2));
maskR1 = false(size(fuse,1),size(fuse,2));
maskR2 = false(size(fuse,1),size(fuse,2));

for i=1:size(obj_chess,1)
    for j=1:size(obj_chess(i).chess,2) %Inew = I.*repmat(M,[1,1,3]);
        if(strcmp(obj_chess(i).chess(j).position,'Left') && strcmp(obj_chess(i).chess(j).type,'Secondary'))
            maskL2 = (maskL2 | obj_chess(i).chess(j).mask);
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
%% Rifinitura Maschere
% Conto quante linee ci devono essere in funzione dei colori
idxLinesCenter = find(strcmp(label(3,:),positions(2)));
numLinesCenter = size(idxLinesCenter,2)+1;
idxLinesLeftSec = find(strcmp(label(3,:),positions(1)) & strcmp(label(4,:),types(2)));
numLinesLeftSec = size(idxLinesLeftSec,2)+1;
idxLinesLeftPri = find(strcmp(label(3,:),positions(1)) & strcmp(label(4,:),types(1)));
numLinesLeftPri = size(idxLinesLeftPri,2)+1;
idxLinesRightSec = find(strcmp(label(3,:),positions(3)) & strcmp(label(4,:),types(2)));
numLinesRightSec = size(idxLinesRightSec,2)+1;
idxLinesRightPri = find(strcmp(label(3,:),positions(3)) & strcmp(label(4,:),types(1)));
numLinesRightPri = size(idxLinesRightPri,2)+1;

LinesCenter = findLines(idxLinesCenter,  obj_chess, order);
LinesCenter = mergeNearestLines(LinesCenter, numLinesCenter, size(maskC));
gap  = false(size(fuse,1),size(fuse,2));
for i=1:size(LinesCenter,2)
    gap = gap | line2image(LinesCenter{i},size(maskC));
end
gap = imdilate(gap,strel('disk',3));
maskC = maskC & ~gap;
maskC = imopen(maskC,strel('square',3));

LinesLeftSec = findLines(idxLinesLeftSec,  obj_chess, order);
LinesLeftSec = mergeNearestLines(LinesLeftSec, numLinesLeftSec, size(maskC));
gap  = false(size(fuse,1),size(fuse,2));
for i=1:size(LinesLeftSec,2)
    gap = gap | line2image(LinesLeftSec{i},size(maskC));
end
gap = imdilate(gap,strel('disk',3));
maskL2 = maskL2 & ~gap;
maskL2 = imopen(maskL2,strel('square',3));

LinesLeftPri = findLines(idxLinesLeftPri,  obj_chess, order);
LinesLeftPri = mergeNearestLines(LinesLeftPri, numLinesLeftPri, size(maskC));
gap  = false(size(fuse,1),size(fuse,2));
for i=1:size(LinesLeftPri,2)
    gap = gap | line2image(LinesLeftPri{i},size(maskC));
end
gap = imdilate(gap,strel('disk',3));
maskL1 = maskL1 & ~gap;
maskL1 = imopen(maskL1,strel('square',3));

LinesRightSec = findLines(idxLinesRightSec,  obj_chess, order);
LinesRightSec = mergeNearestLines(LinesRightSec, numLinesRightSec, size(maskC));
gap  = false(size(fuse,1),size(fuse,2));
for i=1:size(LinesRightSec,2)
    gap = gap | line2image(LinesRightSec{i},size(maskC));
end
gap = imdilate(gap,strel('disk',3));
maskR2 = maskR2 & ~gap;
maskR2 = imopen(maskR2,strel('square',3));

LinesRightPri = findLines(idxLinesRightPri,  obj_chess, order);
LinesRightPri = mergeNearestLines(LinesRightPri, numLinesRightPri, size(maskC));
gap  = false(size(fuse,1),size(fuse,2));
for i=1:size(LinesRightPri,2)
    gap = gap | line2image(LinesRightPri{i},size(maskC));
end
gap = imdilate(gap,strel('disk',3));
maskR1 = maskR1 & ~gap;
maskR1 = imopen(maskR1,strel('square',3));
%% Ripulisco maschere da oggetti piccoli
statsL2 = cell2mat(struct2cell(regionprops(maskC,'Area')));
maskL2 = bwareaopen(maskL2,round(0.1*(mean(statsL2)-std(statsL2))));
figure, imshow(maskL2);
statsL1 = cell2mat(struct2cell(regionprops(maskL1,'Area')));
maskL1 = bwareaopen(maskL1,round(0.1*(mean(statsL1)-std(statsL1))));
figure, imshow(maskL1);
statsC = cell2mat(struct2cell(regionprops(maskC,'Area')));
maskC = bwareaopen(maskC,round(0.1*(mean(statsC)-std(statsC))));
figure, imshow(maskC);
statsR1 = cell2mat(struct2cell(regionprops(maskC,'Area')));
maskR1 = bwareaopen(maskR1,round(0.1*(mean(statsR1)-std(statsR1))));
figure, imshow(maskR1);
statsR2 = cell2mat(struct2cell(regionprops(maskR2,'Area')));
maskR2 = bwareaopen(maskR2,round(0.1*(mean(statsR2)-std(statsR2))));
figure, imshow(maskR2);
%% Calcolo convexhull dei riflessi: controllare se è necessario fare sta cosa
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
%     idx = find(mask == 1);
%     [idx,idy]=ind2sub(size(mask),idx);
%     j = boundary(idx,idy,0.1); % Parametro
%     maskI = poly2mask(idy(j),idx(j), size(mask,1), size(mask,2));
%     DilatedMask = imdilate(maskI,strel('rectangle',[3,20]));
%     idx = find(maskI == 1);
%     [idx,idy]=ind2sub(size(mask),idx);
%     maskI = false(size(maskI));
%     maskI(min(idx):max(idx),min(idy):max(idy)) = DilatedMask(min(idx):max(idx),min(idy):max(idy));
    maskI = bwconvhull(mask);
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
%% Calcolo Convexhull tessere singole
for i=1:size(obj_chess,1)
    for j=1:size(obj_chess(i).chess,2)
        cut_x = obj_chess(i).bbox_x(j,:);
        cut_y = obj_chess(i).bbox_y(j,:); 
        obj_chess(i).chess(j).ch_mask = false(size(fuse,1),size(fuse,2));
        mask = imdilate(bwconvhull(obj_chess(i).chess(j).mask),strel('square',3));
        obj_chess(i).chess(j).ch_mask(cut_y(1):cut_y(2),cut_x(1):cut_x(2)) = mask(cut_y(1):cut_y(2),cut_x(1):cut_x(2));     
        figure, imshowpair(obj_chess(i).chess(j).ch_mask,rgb2gray(I),'falsecolor');
    end
end
%% Full Mask
FullMask = maskCI+maskL1I+maskL2I+maskR1I+maskR2I;
BlurredMask = imgaussfilt(FullMask,10);
%% Fissa il centro degli assi
R = im2double(I(:,:,1)); G = im2double(I(:,:,2)); B = im2double(I(:,:,3));
MINRGB = min(R,G); MINRGB = min(MINRGB,B);
MINMASK = false(size(MINRGB));
for i=1:size(obj_chess,1)
    for j=1:size(obj_chess(i).chess,2) 
        if(~strcmp(obj_chess(i).name,'yellow'))
            if(~strcmp(obj_chess(i).chess(j).background,'Black'))
                tmp = obj_chess(i).chess(j).mask & obj_chess(i).color_mask;
                MINMASK = MINMASK | tmp;
            end
        end
    end
end
MINRGB = imgaussfilt((MINRGB+MINMASK).*BlurredMask,2);
T = graythresh(MINRGB);
bw = im2bw(MINRGB,T);
RIGHTMASK = bwareaopen(bw, 2*Threshold); % Parametro
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
LEFTMASK(LEFTMASK<0) = 0;
LEFTMASK = imopen(LEFTMASK,strel('square',10));
LEFTMASK = imclose(LEFTMASK,strel('square',5));
CC = bwconncomp(LEFTMASK,4); CC_dim = zeros(1,size(CC.PixelIdxList,2));
for i=1:size(CC.PixelIdxList,2)
    CC_dim(i) = size(CC.PixelIdxList{i},1);
end
LEFTMASK = false(size(LEFTMASK));
for i=1:CC.NumObjects
    idx = find(CC_dim == max(CC_dim));
    LEFTMASK(CC.PixelIdxList{idx}) = 1;
    CC_dim(idx) = 0;
end
left = imdilate(edge(LEFTMASK),strel('disk',3));
right = imdilate(edge(RIGHTMASK),strel('disk',3));
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

color_filtered = abs(imfilter(im2double(TheColors+left_center_axis+right_center_axis+mid_center_axis),fspecial('prewit')));
key_filtered = absimfilter(KEY,fspecial('prewit'));
orig_filtered = absimfilter(ORIG,fspecial('prewit'));
img_filtered = (key_filtered+orig_filtered+color_filtered)./3;
bw = im2bw(img_filtered,graythresh(img_filtered.*FullMask)); % Se si cambia la soglia qua non funge più un cazzo
bw = imclose(bw,strel('square',5));
bwV = bwareaopen(bw, 100); % Parametro
for i=1:5
    switch (i)
        case 1
            % Left Secondary Reflection
            idl = find(strcmp(label(3,:),'Left') & strcmp(label(4,:),'Secondary'));
            grid_pos = 'Left';
            grid_typ = 'Secondary';
            maskI = maskL2I;
            mask = maskL2;
        case 2
            % Right Secondary Reflection
            idl = find(strcmp(label(3,:),'Right') & strcmp(label(4,:),'Secondary'));
            grid_pos = 'Right';
            grid_typ = 'Secondary';
            maskI = maskR2I;
            mask = maskR2;
        case 3
            % Left Primary Reflection
            idl = find(strcmp(label(3,:),'Left') & strcmp(label(4,:),'Primary'));
            grid_pos = 'Left';
            grid_typ = 'Primary';
            maskI = maskL1I;
            mask = maskL1;
        case 4
            % Right Primary Reflection
            idl = find(strcmp(label(3,:),'Right') & strcmp(label(4,:),'Primary'));
            grid_pos = 'Right';
            grid_typ = 'Primary';
            maskI = maskR1I;
            mask = maskR1;
        case 5
            % Center Real Reflection
            idl = find(strcmp(label(3,:),'Center') & strcmp(label(4,:),'Real'));
            grid_pos = 'Center';
            grid_typ = 'Real';
            maskI = maskCI;
            mask = maskC;
    end   
    id_obj_chess = order(1,idl);
    id_chess = order(2,idl);
    gapGRID = line2image(obj_chess(id_obj_chess(1)).chess(id_chess(1)).v_lines{1},size(mask));
    gapGRID = gapGRID | line2image(obj_chess(id_obj_chess(end)).chess(id_chess(end)).v_lines{end},size(mask));
    [GRIDhh, GRIDh] = getHImage(grid_pos,grid_typ,obj_chess, order, label, size(maskI), maskI);
    GRIDvv = getVImage(bwV, maskI, GRIDh);
    GRID = imerode(maskI,strel('rectangle',[30,10])) & ~mask & ~imdilate(gapGRID,strel('disk',12)) | GRIDvv;
    % ordina i blob in base alla distanza dalla retta destra e sinistra,
    % misura la distanza media, usa questa per definire di quanto allargare
    % le rette che tagliano il bordo più esterno
    GRIDv = abs(imfilter(im2double(GRID),[-1 0 1]));
        
    GRIDh = imdilate(GRIDh,strel('disk',10)).*maskI;
    stats = cell2mat(struct2cell(regionprops(GRIDh,'Area')));
    GRIDh = imopen(bwareaopen(GRIDh,round(0.1*(mean(stats)-std(stats)))),strel('disk',10));

    GRIDv = GRIDv - GRIDh;
    GRIDv(GRIDv<0) = 0;
    stats = cell2mat(struct2cell(regionprops(bwconncomp(GRIDv,8),'Area')));
    GRIDv = bwareaopen(GRIDv,50); % Aggiustare questa soglia
    
    figure, imshow(GRIDv);
    
    props_center = regionprops(GRIDv,'Centroid'); % [X; Y]
    props_center = reshape(cell2mat(struct2cell(props_center)),2,size(props_center,1));
    props = [1:size(props_center,2);props_center];
    [props(2,:), props(1,:)] = sort(props(2,:));
    props(3,:) = props(3, props(1,:));
    CC = bwconncomp(GRIDv,8); usu_center = [zeros(1,size(props_center,2));props_center]; usu = CC.NumObjects;
    figure, imshow(I); hold on; image_line = false(size(GRIDv));
    for j=1:size(idl,2)
        for l=1:size(obj_chess(id_obj_chess(j)).chess(id_chess(j)).v_lines_centroid,2)
            line = obj_chess(id_obj_chess(j)).chess(id_chess(j)).v_lines_centroid{l};
            lines = [];
            if ~isempty(line)
                for k=1:CC.NumObjects
                    if(usu_center(1,k)==0)
                        if(isLeft(line,usu_center(2,k),usu_center(3,k)))
                            usu_center(1,k) = 1;
                            usu = usu - 1;
                            lines = [lines; CC.PixelIdxList{k}];
                        end
                    end
                end
                if(~isempty(lines))
                    [liney,linex] = ind2sub(size(GRIDv),lines); 
                    [linefitresult, ~] = createLineInv(liney, linex, size(GRIDv));
                    coeffs = coeffvalues(linefitresult);
                    x = 1:0.001:size(GRIDv,2);
                    y = floor(polyval(coeffs,x));
                    x = floor(x);
                    x(y<1 | y>size(GRIDv,1)) = [];
                    y(y<1 | y>size(GRIDv,1)) = [];
                    for count = 1:size(x,2)
                        image_line(y(count),x(count)) = 1;
                    end
                    plot(linefitresult);
                end
            end
        end
    end
    lines = [];
    if usu > 0
        for k=1:CC.NumObjects
            if(usu_center(1,k)==0)
                if(~isLeft(line,usu_center(2,k),usu_center(3,k)))
                    usu_center(1,k) = 1;
                    usu = usu - 1;
                    lines = [lines; CC.PixelIdxList{k}];
                end
            end
        end
        if(~isempty(lines))
            [liney,linex] = ind2sub(size(GRIDv),lines); 
            [linefitresult, ~] = createLineInv(liney, linex, size(GRIDv));
            coeffs = coeffvalues(linefitresult);
            x = 1:0.001:size(GRIDv,2);
            y = floor(polyval(coeffs,x));
            x = floor(x);
            x(y<1 | y>size(GRIDv,1)) = [];
            y(y<1 | y>size(GRIDv,1)) = [];
            for count = 1:size(x,2)
                image_line(y(count),x(count)) = 1;
            end
            plot(linefitresult);
        end
    end
% Orizzontale    
    GRIDv2 = imdilate(image_line,strel('disk',10));
    gray = im2double(rgb2gray(I));
    sss = gray.*GRIDv2.*imdilate(GRIDh,strel('disk',10));
    sss = sss.*bwareaopen(sss>0,10);
    CC = bwconncomp(sss>0,8);
    sss_image = zeros(size(sss,1),size(sss,2),3); P = [];
    for l=1:CC.NumObjects
        [cutR,cutC] = ind2sub(size(GRIDv),CC.PixelIdxList{l}); 
        tmp = sss(min(cutR):max(cutR),min(cutC):max(cutC));
        tt = 0.3*size(tmp,1)*size(tmp,2);
        condition = true;
        th = 0.9; store = cell(2,1);
        dot = findDots(tmp);
        %figure, imshow(tmp), hold on, scatter(dot(:,1),dot(:,2))
        P = [P;min(cutR)+dot(:,2),min(cutC)+dot(:,1)];
        if(isempty(dot))
            while(condition && th>0)
                bw_tmp = im2bw(tmp,th);
                th = th-0.01;
                CC_tmp = bwconncomp(bw_tmp,4);        
                if(CC_tmp.NumObjects ==1 && size(CC_tmp.PixelIdxList{1},1)>tt)
                    condition = false;
                    [blob1y,blob1x] = ind2sub(size(tmp),store{1});
                    [blob2y,blob2x] = ind2sub(size(tmp),store{2});
                    D = zeros(5,size(blob1y,1)*(size(blob2y,1)-1)); idist = 1;
                    for k=1:size(blob1y)
                        for j=1:size(blob2y)
                            D(1,idist) = blob1x(k);
                            D(2,idist) = blob1y(k);
                            D(3,idist) = blob2x(j);
                            D(4,idist) = blob2y(j);
                            D(5,idist) = pdist([blob1x(k) blob1y(k); blob2x(j) blob2y(j)],'euclidean');
                            idist = idist+1;
                        end
                    end
                    if(isempty(D))
                         condition = true;
                    end
                end
                if(CC_tmp.NumObjects == 2)
                   store(1) = CC_tmp.PixelIdxList(1);
                   store(2) = CC_tmp.PixelIdxList(2);          
                end 
            end
            dot = findDots(bw_tmp);
            if(isempty(dot))
                [D(5,:), idd] = sort(D(5,:));
                D(1,:) = D(1,idd);
                D(2,:) = D(2,idd);
                D(3,:) = D(3,idd);
                D(4,:) = D(4,idd);
                XX = round(mean([D(1,1:5),D(3,1:5)]));
                YY = round(mean([D(2,1:5),D(4,1:5)]));
                %figure, imshow(bw_tmp), hold on, scatter(XX,YY);
                sss_image(min(cutR):max(cutR),min(cutC):max(cutC),1) = bw_tmp;
                sss_image(min(cutR)+YY,min(cutC)+XX,2) = 1;
            else
                YY = dot(2);
                XX = dot(1);
            end
            P = [P;min(cutR)+YY,min(cutC)+XX];
        end
    end
    figure, imshow(I), hold on, scatter(P(:,2),P(:,1));
end

% CL2v = GRIDv;
% 
% line1 = [CC.PixelIdxList{props(1,1)};CC.PixelIdxList{props(1,2)};CC.PixelIdxList{props(1,3)};CC.PixelIdxList{props(1,4)};CC.PixelIdxList{props(1,5)};CC.PixelIdxList{props(1,6)};CC.PixelIdxList{props(1,7)}];
% [line1y,line1x] = ind2sub(size(CL2v),line1); 
% [line1_fitresult, ~] = createLineInv(line1y, line1x, size(CL2v));
% coeffs = coeffvalues(line1_fitresult);
% dim = size(CL2v);
% x = 1:0.001:size(CL2v,1);
% y = floor(polyval(coeffs,x));
% x = floor(x);
% x(y<1 | y>dim(2)) = [];
% y(y<1 | y>dim(2)) = [];
% image_line = false(dim);
% for j = 1:size(x,2)
%     image_line(x(j),y(j)) = 1;
% end
% imshow(I); hold on;
% plot(y,x,'r');
% line2 = [CC.PixelIdxList{props(1,11)};CC.PixelIdxList{props(1,12)};CC.PixelIdxList{props(1,13)};CC.PixelIdxList{props(1,14)};CC.PixelIdxList{props(1,15)};CC.PixelIdxList{props(1,16)};CC.PixelIdxList{props(1,17)};CC.PixelIdxList{props(1,18)};CC.PixelIdxList{props(1,19)};CC.PixelIdxList{props(1,20)}];
% [line1y,line1x] = ind2sub(size(CL2v),line2); 
% [line1_fitresult, ~] = createLineInv(line1y, line1x);
% coeffs = coeffvalues(line1_fitresult);
% x = 1:0.001:size(CL2v,1);
% y = floor(polyval(coeffs,x));
% x = floor(x);
% x(y<1 | y>dim(2)) = [];
% y(y<1 | y>dim(2)) = [];
% for j = 1:size(x,2)
%     image_line(x(j),y(j)) = 1;
% end
% plot(y,x,'r');
% line3 = [CC.PixelIdxList{props(1,21)};CC.PixelIdxList{props(1,22)};CC.PixelIdxList{props(1,23)};CC.PixelIdxList{props(1,24)};CC.PixelIdxList{props(1,25)};CC.PixelIdxList{props(1,26)};CC.PixelIdxList{props(1,27)};CC.PixelIdxList{props(1,28)};CC.PixelIdxList{props(1,29)};CC.PixelIdxList{props(1,30)}];
% [line1y,line1x] = ind2sub(size(CL2v),line3); 
% [line1_fitresult, ~] = createLineInv(line1y, line1x);
% coeffs = coeffvalues(line1_fitresult);
% x = 1:0.001:size(CL2v,1);
% y = floor(polyval(coeffs,x));
% x = floor(x);
% x(y<1 | y>dim(2)) = [];
% y(y<1 | y>dim(2)) = [];
% for j = 1:size(x,2)
%     image_line(x(j),y(j)) = 1;
% end
% hold on;
% plot(y,x,'r');
% line4 = [CC.PixelIdxList{props(1,31)};CC.PixelIdxList{props(1,32)};CC.PixelIdxList{props(1,33)};CC.PixelIdxList{props(1,34)};CC.PixelIdxList{props(1,35)};CC.PixelIdxList{props(1,36)};CC.PixelIdxList{props(1,37)};CC.PixelIdxList{props(1,38)};CC.PixelIdxList{props(1,39)}; CC.PixelIdxList{props(1,40)}];
% [line1y,line1x] = ind2sub(size(CL2v),line4); 
% [line1_fitresult, ~] = createLineInv(line1y, line1x);
% coeffs = coeffvalues(line1_fitresult);
% x = 1:0.001:size(CL2v,1);
% y = floor(polyval(coeffs,x));
% x = floor(x);
% x(y<1 | y>dim(2)) = [];
% y(y<1 | y>dim(2)) = [];
% for j = 1:size(x,2)
%     image_line(x(j),y(j)) = 1;
% end
% hold on;
% plot(y,x,'r');

CL2v2 = imdilate(image_line,strel('disk',15));
gray = im2double(rgb2gray(I));
sss = gray.*CL2v2.*imerode(maskL2I,strel('rectangle',[50,2])).*imdilate(CL2h,strel('disk',10));
sss = sss.*bwareaopen(sss>0,10);
CC = bwconncomp(sss>0,8);
sss_image = zeros(size(sss,1),size(sss,2),3); P = [];
for i=1:CC.NumObjects
    [cutR,cutC] = ind2sub(size(CL2v),CC.PixelIdxList{i}); 
    tmp = sss(min(cutR):max(cutR),min(cutC):max(cutC));
    tt = 0.3*size(tmp,1)*size(tmp,2);
    condition = true;
    th = 0.9; store = cell(2,1);
    while(condition)
        bw_tmp = im2bw(tmp,th);
        th = th-0.01;
        CC_tmp = bwconncomp(bw_tmp,8);        
        if(CC_tmp.NumObjects ==1 && size(CC_tmp.PixelIdxList{1},1)>tt)
            condition = false;
            [blob1y,blob1x] = ind2sub(size(tmp),store{1});
            [blob2y,blob2x] = ind2sub(size(tmp),store{2});
            D = zeros(5,size(blob1y,1)*(size(blob2y,1)-1)); idist = 1;
            for k=1:size(blob1y)
                for j=1:size(blob2y)
                    D(1,idist) = blob1x(k);
                    D(2,idist) = blob1y(k);
                    D(3,idist) = blob2x(j);
                    D(4,idist) = blob2y(j);
                    D(5,idist) = pdist([blob1x(k) blob1y(k); blob2x(j) blob2y(j)],'euclidean');
                    idist = idist+1;
                end
            end
        end
        if(CC_tmp.NumObjects == 2)
           store(1) = CC_tmp.PixelIdxList(1);
           store(2) = CC_tmp.PixelIdxList(2);          
        end 
    end
    [D(5,:), idd] = sort(D(5,:));
    D(1,:) = D(1,idd);
    D(2,:) = D(2,idd);
    D(3,:) = D(3,idd);
    D(4,:) = D(4,idd);
    XX = round(mean([D(1,1:5),D(3,1:5)]));
    YY = round(mean([D(2,1:5),D(4,1:5)]));
    figure, imshow(bw_tmp), hold on, scatter(XX,YY);
    sss_image(min(cutR):max(cutR),min(cutC):max(cutC),1) = bw_tmp;
    sss_image(min(cutR)+YY,min(cutC)+XX,2) = 1;
    P = [P;min(cutR)+YY,min(cutC)+XX];
end

figure, imshow(I); hold on; scatter(P(:,2),P(:,1))
C = regionprops(sss_image,'Extrema');
imshow(sss_image);
hold on
CC_e = 1.0e+03 *[
   
    4.3775    0.1245
    4.3805    0.1245
    4.3955    0.1465
    4.3955    0.1625
    4.3865    0.1635
    4.3825    0.1635
    4.3645    0.1415
    4.3645    0.1265];
scatter(CC_e(:,1),CC_e(:,2));
hist = imhist(sss(CC.PixelIdxList{1}));
hist = [1:256;hist'];


T = adaptthresh(sss,0,'ForegroundPolarity','bright');
BW = imbinarize(sss,T);


CL2h2 = imfilter(CL2,[-1 0 1]') | imfilter(CL2,[1 0 -1]');
CL2v2 = imdilate(image_line,strel('disk',15));
CL2h2 = CL2 - CL2v2 - CL2v;
CL2h2(CL2h2<0) = 0;
CL2h2 = bwareaopen(CL2h2,20);
CL2h2 = imopen(CL2h2,strel('disk',5));

props = regionprops(CL2h2,'Centroid');
props = reshape(cell2mat(struct2cell(props)),2,size(props,1));
props = [1:size(props,2);props];
[props(3,:), props(1,:)] = sort(props(3,:));
props(2,:) = props(2, props(1,:));
CC = bwconncomp(CL2h2,8);

line1 = [CC.PixelIdxList{props(1,1)};CC.PixelIdxList{props(1,2)};CC.PixelIdxList{props(1,3)};CC.PixelIdxList{props(1,4)}];
[line1y,line1x] = ind2sub(size(CL2v),line1); 
[line1_fitresult, ~] = create2poly(line1x, line1y);
coeffs = coeffvalues(line1_fitresult);
dim = size(CL2v);
x = 1:0.001:size(CL2v,1);
y = polyval(coeffs,x);
plot(x,y,'r');
line2 = [CC.PixelIdxList{props(1,9)};CC.PixelIdxList{props(1,10)};CC.PixelIdxList{props(1,11)};CC.PixelIdxList{props(1,12)};CC.PixelIdxList{props(1,13)};CC.PixelIdxList{props(1,14)};CC.PixelIdxList{props(1,15)};CC.PixelIdxList{props(1,16)}];
[line1y,line1x] = ind2sub(size(CL2v),line2); 
[line1_fitresult, ~] = createLineInv(line1y, line1x);
coeffs = coeffvalues(line1_fitresult);
x = 1:0.001:size(CL2v,1);
y = polyval(coeffs,x);
plot(y,x,'r');
line3 = [CC.PixelIdxList{props(1,17)};CC.PixelIdxList{props(1,18)};CC.PixelIdxList{props(1,19)};CC.PixelIdxList{props(1,20)};CC.PixelIdxList{props(1,21)};CC.PixelIdxList{props(1,22)};CC.PixelIdxList{props(1,23)};CC.PixelIdxList{props(1,24)}];
[line1y,line1x] = ind2sub(size(CL2v),line3); 
[line1_fitresult, ~] = createLineInv(line1y, line1x);
coeffs = coeffvalues(line1_fitresult);
x = 1:0.001:size(CL2v,1);
y = polyval(coeffs,x);
hold on;
plot(y,x,'r');
line4 = [CC.PixelIdxList{props(1,25)};CC.PixelIdxList{props(1,26)};CC.PixelIdxList{props(1,27)};CC.PixelIdxList{props(1,28)};CC.PixelIdxList{props(1,29)};CC.PixelIdxList{props(1,30)};CC.PixelIdxList{props(1,31)};CC.PixelIdxList{props(1,32)}];
[line1y,line1x] = ind2sub(size(CL2v),line4); 
[line1_fitresult, ~] = createLineInv(line1y, line1x);
coeffs = coeffvalues(line1_fitresult);
x = 1:0.001:size(CL2v,1);
y = polyval(coeffs,x);
hold on;
plot(y,x,'r');


tmp = CL2h2;
for i=1:10
    tmp(CC.PixelIdxList{props(1,i)}) = 0;
    imshow(tmp)
    i
    pause
end


%% Corner
cymk = rgb2cmyk(I);
Cyano = im2double(cymk(:,:,1));
Magenta = im2double(cymk(:,:,2));
Yellow = im2double(cymk(:,:,3));
TheColors = Cyano+Magenta+Yellow;
Key = im2double(cymk(:,:,4));

KEY = (1-Key);
ORIG = rgb2gray(im2double(O));
[color_filtered, color_filtered_x, color_filtered_y] = absimfilter(TheColors+left_center_axis+right_center_axis+mid_center_axis,fspecial('prewit'));
[key_filtered, key_filtered_x, key_filtered_y] = absimfilter(KEY,fspecial('prewit'));
[orig_filtered, orig_filtered_x, orig_filtered_y] = absimfilter(ORIG,fspecial('prewit'));
% 
% [eeee, eeeee, eeeeee] = absimfilter(L,fspecial('sobel'));
% bw = im2bw(eeee,0.3*graythresh(eeee.*FullMask));
% bw = imdilate(bw,strel('square',2));
% bw = imerode(bw,strel('square',5));

img_filtered = (key_filtered+orig_filtered+color_filtered)./3;

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
bw = imclose(bw,strel('square',5));
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
    for i = 1:size(CCh.PixelIdxList,2)
        [idy,idx] = ind2sub(size(I),CCh.PixelIdxList{i});  % Ma è giusta sta size? I è a 3 dimensioni    
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





    % Allineamento a sinistra/destra delle matrici: operazione necessaria
    % prima del flip per generare le corrispondenze.
    % Flip a Destra: I valori nella matrice vanno allineati a sinistra.
    % Flip a Sinistra: I valori nella matrice vanno allineati a destra.
    % L'allineamento non può essere permanente, ma deve essere temporanemo
    % dipende dalla direzione in cui si flippa 
    % si portrebbe pensare di avere due matrici una allieata a sininstra
    % per le rotazioni a destra, una allineata a destra per le rotazioni a
    % sinistra
    % inutile memorizzarle tutte e due, ne memorizziamo una l'altra si
    % ricava e si usa temporaneamente per creare le associazioni
    % per le associazioni utilizziamo un grafo, uno per ogni coppia
    % colore/sfondo

% for l=1:size(obj_chess,1)
%     idx_color = find(strcmp(label(1,:),obj_chess(l).name));
%     idx_bgwhite = find(strcmp(label(2,:),'White'));
%     idx_bgblack = find(strcmp(label(2,:),'Black'));
%     idx_center = find(strcmp(label(3,:),positions(2)));
%     idx_bgwhite = intersect(idx_bgwhite,idx_color);
%     idx_bgblack = intersect(idx_bgblack,idx_color);
%     idx_center_bgwhite = intersect(idx_center,idx_bgwhite);
%     idx_center_bgblack = intersect(idx_center,idx_bgblack);
% 
%     % Associa la centrale con le laterali
%     for i=1:size(idx_bgwhite,2)
%         if(idx_bgwhite(i) ~= idx_center_bgwhite)
%             label_copy = label(:,idx_bgwhite(i)); % Colore Foreground, Colore Background, Posizione, Tipo
%             order_copy = order(:,idx_bgwhite(i)); % Indice nel vettore obj_chess, Indice Chess, Centroide
%             logical_matrix = logical(obj_chess(l).chess(order_copy(2)).intersections_x);
%             logical_center = logical(obj_chess(l).chess(order(2,idx_center_bgwhite)).intersections_x);
%                         hasLeft(label,idx_bgwhite(i));
%             hasRight(label,idx_bgwhite(i));
%             
%             
%             string = strcat(label_copy(3),label_copy(4));
%             if iscell(string)
%                 string = string{1};
%             end
%             string = lower(string) 
%             obj_chess(l).name
%             
%             
%             switch string
%                 case 'leftsecondary'
%                     [logical_matrix,shift_matrix] = allignright(logical_matrix);
%                     index_logical_matrix = 1:size(logical_matrix,2);
%                     [logical_center,shift_center] = allignleft(logical_center);
%                     index_logical_center = index_logical_matrix;
%                     overlap_matrix = logical_matrix & logical_center;
%                     
%                     index_logical_center(overlap_matrix(1,:) == 0) = [];
%                     index_logical_matrix(overlap_matrix(1,:) == 0) = [];
%                     
%                     index_logical_center = index_logical_center - shift_center;
%                     index_logical_matrix = index_logical_matrix - shift_matrix;
%                     corrispondency = [index_logical_center;index_logical_matrix]
%                 case 'leftprimary'
%                     [logical_matrix,shift_matrix] = allignright(logical_matrix);
%                     [logical_center,shift_center] = allignleft(logical_center);
%                     logical_center = fliplr(logical_center);                    
%                     index_logical_matrix = 1:size(logical_matrix,2);
%                     index_logical_center = fliplr(1:size(logical_center,2));                
%                     overlap_matrix = logical_matrix & logical_center;
%                     
%                     index_logical_center(overlap_matrix(1,:) == 0) = [];
%                     index_logical_matrix(overlap_matrix(1,:) == 0) = [];
%                     
%                     index_logical_center = index_logical_center + shift_center;
%                     index_logical_matrix = index_logical_matrix - shift_matrix;
%                     corrispondency = [index_logical_center;index_logical_matrix]
%                 case 'rightsecondary'
%                     [logical_matrix,shift_matrix] = allignleft(logical_matrix);
%                     index_logical_matrix = 1:size(logical_matrix,2);
%                     [logical_center,shift_center] = allignleft(logical_center);
%                     index_logical_center = index_logical_matrix;
%                     overlap_matrix = logical_matrix & logical_center;
%                     
%                     index_logical_center(overlap_matrix(1,:) == 0) = [];
%                     index_logical_matrix(overlap_matrix(1,:) == 0) = [];
%                     
%                     index_logical_center = index_logical_center - shift_center;
%                     index_logical_matrix = index_logical_matrix - shift_matrix; % non dovrebbe essere con il +?
%                     corrispondency = [index_logical_center;index_logical_matrix]    
%                 case 'rightprimary'
%                     [logical_matrix,shift_matrix] = allignleft(logical_matrix);
%                     [logical_center,shift_center] = allignright(logical_center);
%                     logical_center = fliplr(logical_center);                    
%                     index_logical_matrix = 1:size(logical_matrix,2);
%                     index_logical_center = fliplr(1:size(logical_center,2));                  
%                     overlap_matrix = logical_matrix & logical_center;
%                     
%                     index_logical_center(overlap_matrix(1,:) == 0) = [];
%                     index_logical_matrix(overlap_matrix(1,:) == 0) = [];
%                     
%                     index_logical_center = index_logical_center + shift_center;
%                     index_logical_matrix = index_logical_matrix - shift_matrix;
%                     corrispondency = [index_logical_center;index_logical_matrix]
%             end
%         end 
%         pause
%     end
%     for i=1:size(idx_bgblack,2)
%         if(idx_bgblack(i) ~= idx_center_bgblack)
%             label_copy = label(:,idx_bgblack(i)); % Colore Foreground, Colore Background, Posizione, Tipo
%             order_copy = order(:,idx_bgblack(i)); % Indice nel vettore obj_chess, Indice Chess, Centroide
%             logical_matrix = logical(obj_chess(l).chess(order_copy(2)).intersections_x);
%             logical_center = logical(obj_chess(l).chess(order(2,idx_center_bgblack)).intersections_x);
%             string = strcat(label_copy(3),label_copy(4));
%             if iscell(string)
%                 string = string{1};
%             end
%             string = lower(string) 
%             obj_chess(l).name 
%             switch string
%                 case 'leftsecondary'
%                     [logical_matrix,shift_matrix] = allignright(logical_matrix);
%                     index_logical_matrix = 1:size(logical_matrix,2);
%                     [logical_center,shift_center] = allignright(logical_center);
%                     index_logical_center = index_logical_matrix;
%                     overlap_matrix = logical_matrix & logical_center;
%                     
%                     index_logical_center(overlap_matrix(1,:) == 0) = [];
%                     index_logical_matrix(overlap_matrix(1,:) == 0) = [];
%                     
%                     index_logical_center = index_logical_center - shift_center;
%                     index_logical_matrix = index_logical_matrix - shift_matrix; % non dovrebbe essere con il +?
%                     corrispondency = [index_logical_center;index_logical_matrix]
%                 case 'leftprimary'
%                     [logical_matrix,shift_matrix] = allignleft(logical_matrix);
%                     [logical_center,shift_center] = allignright(logical_center);
%                     logical_center = fliplr(logical_center);                    
%                     index_logical_matrix = 1:size(logical_matrix,2);
%                     index_logical_center = fliplr(1:size(logical_center,2));                  
%                     overlap_matrix = logical_matrix & logical_center;
%                     
%                     index_logical_center(overlap_matrix(1,:) == 0) = [];
%                     index_logical_matrix(overlap_matrix(1,:) == 0) = [];
%                     
%                     index_logical_center = index_logical_center + shift_center;
%                     index_logical_matrix = index_logical_matrix - shift_matrix;
%                     corrispondency = [index_logical_center;index_logical_matrix]
%                 case 'rightsecondary'
%                     [logical_matrix,shift_matrix] = allignleft(logical_matrix);
%                     index_logical_matrix = 1:size(logical_matrix,2);
%                     [logical_center,shift_center] = allignleft(logical_center);
%                     index_logical_center = index_logical_matrix;
%                     overlap_matrix = logical_matrix & logical_center;
%                     
%                     index_logical_center(overlap_matrix(1,:) == 0) = [];
%                     index_logical_matrix(overlap_matrix(1,:) == 0) = [];
%                     
%                     index_logical_center = index_logical_center - shift_center;
%                     index_logical_matrix = index_logical_matrix - shift_matrix; % non dovrebbe essere con il +?
%                     corrispondency = [index_logical_center;index_logical_matrix]                    
%                 case 'rightprimary'
%                     [logical_matrix,shift_matrix] = allignleft(logical_matrix);
%                     [logical_center,shift_center] = allignright(logical_center);
%                     logical_center = fliplr(logical_center);                    
%                     index_logical_matrix = 1:size(logical_matrix,2);
%                     index_logical_center = fliplr(1:size(logical_center,2));                  
%                     overlap_matrix = logical_matrix & logical_center;
%                     
%                     index_logical_center(overlap_matrix(1,:) == 0) = [];
%                     index_logical_matrix(overlap_matrix(1,:) == 0) = [];
%                     
%                     index_logical_center = index_logical_center + shift_center;
%                     index_logical_matrix = index_logical_matrix - shift_matrix;
%                     corrispondency = [index_logical_center;index_logical_matrix]
%             end
%         end
%         pause
%     end
% end
% 
% 
% 
% 

%--------------------------------------------------------------------------

