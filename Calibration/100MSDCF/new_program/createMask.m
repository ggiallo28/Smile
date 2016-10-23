function [obj] = createMask(RGB, band)
%createMask  Threshold RGB image using auto-generated code from colorThresholder app.
%  [BW,MASKEDRGBIMAGE] = createMask(RGB) thresholds image RGB using
%  auto-generated code from the colorThresholder App. The colorspace and
%  minimum/maximum values for each channel of the colorspace were set in the
%  App and result in a binary mask BW and a composite image maskedRGBImage,
%  which shows the original RGB image values under the mask BW.
warning off
% Auto-generated by colorThresholder app on 08-Oct-2016
%------------------------------------------------------
% Convert RGB image to chosen color space
I = rgb2hsv(RGB);

switch band
    case 'red' 
        RGB = im2double(RGB);
        cform = makecform('srgb2lab', 'AdaptedWhitePoint', whitepoint('D65'));
        I = applycform(RGB,cform);
%         % Define thresholds for channel 1 based on histogram settings
%         channel1Min = 10.876;
%         channel1Max = 75.798;
%         % Define thresholds for channel 2 based on histogram settings
%         channel2Min = 6.705;
%         channel2Max = 76.596;
%         % Define thresholds for channel 3 based on histogram settings
%         channel3Min = 3.081;
%         channel3Max = 70.703;
%         BW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
%              (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
%              (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);
        % Define thresholds for channel 1 based on histogram settings
        channel1Min = 10.171;
        channel1Max = 75.909;

        % Define thresholds for channel 2 based on histogram settings
        channel2Min = 6.703;
        channel2Max = 77.961;

        % Define thresholds for channel 3 based on histogram settings
        channel3Min = -20.494;
        channel3Max = 70.297;

        % Create mask based on chosen histogram thresholds
        BW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
            (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
            (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);
        % Initialize output masked image based on input image.
        BW = bwareaopen(BW, 300); % Parametro
        BW2 =  imdilate(BW,strel('square',10));
        BW2 = bwareaopen(BW2, 2000); % Parametro
        BW = BW & BW2;
        BW = imerode(BW,strel('square',5));
        BW = imdilate(BW,strel('square',5));
        BW = bwareaopen(BW, 300); % Parametro
        CC = bwconncomp(BW);
        for i=1:size(CC.PixelIdxList,2)
            BW_TMP = false(size(BW));
            BW_TMP(CC.PixelIdxList{i}) = 1;
            s = regionprops(BW_TMP,'BoundingBox');
            img = imcrop(BW_TMP,s.BoundingBox);
            white = sum(sum(img));
            black = size(img,1)*size(img,2);
            %figure, imshow(img);
            if(white<0.3*black) % Portato a 0.3 prima era 0.2
               BW(CC.PixelIdxList{i}) = 0; 
            end
        end
        maskedRGBImage = RGB;
        % Set background pixels where BW is false to zero.
        maskedRGBImage(repmat(~BW,[1 1 3])) = 0;
        R = maskedRGBImage(:,:,1);
        G = maskedRGBImage(:,:,2);
        B = maskedRGBImage(:,:,3);
        hold off
        R(BW) = 255;
        G(BW) = 0;
        B(BW) = 0;
        maskedRGBImage(:,:,1)=R;
        maskedRGBImage(:,:,2)=G;
        maskedRGBImage(:,:,3)=B;
        RGB = im2uint8(RGB);
    case 'green'
        I = rgb2hsv(RGB);
        % Define thresholds for channel 1 based on histogram settings
        channel1Min = 0.220;
        channel1Max = 0.463;
        % Define thresholds for channel 2 based on histogram settings
        channel2Min = 0.261;
        channel2Max = 1.000;
        % Define thresholds for channel 3 based on histogram settings
        channel3Min = 0.1;
        channel3Max = 1.000;
        % Create mask based on chosen histogram thresholds
        BW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
             (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
             (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);
        % Initialize output masked image based on input image.
        BW = bwareaopen(BW, 300); % Parametro
        BW2 =  imdilate(BW,strel('square',10));
        BW2 = bwareaopen(BW2, 2000); % Parametro
        BW = BW & BW2;
        BW = imerode(BW,strel('square',5));
        BW = imdilate(BW,strel('square',5));
        BW = bwareaopen(BW, 300); % Parametro
        CC = bwconncomp(BW);
        for i=1:size(CC.PixelIdxList,2)
            BW_TMP = false(size(BW));
            BW_TMP(CC.PixelIdxList{i}) = 1;
            s = regionprops(BW_TMP,'BoundingBox');
            img = imcrop(BW_TMP,s.BoundingBox);
            white = sum(sum(img));
            black = size(img,1)*size(img,2);
            %figure, imshow(img);
            if(white<0.2*black)
               BW(CC.PixelIdxList{i}) = 0; 
            end
        end
        maskedRGBImage = RGB;
        % Set background pixels where BW is false to zero.
        maskedRGBImage(repmat(~BW,[1 1 3])) = 0;
        R = maskedRGBImage(:,:,1);
        G = maskedRGBImage(:,:,2);
        B = maskedRGBImage(:,:,3);
        R(BW) = 0;
        G(BW) = 255;
        B(BW) = 0;
        maskedRGBImage(:,:,1)=R;
        maskedRGBImage(:,:,2)=G;
        maskedRGBImage(:,:,3)=B;
    case 'blue'
        RGB = im2double(RGB);
        cform = makecform('srgb2lab', 'AdaptedWhitePoint', whitepoint('D65'));
        I = applycform(RGB,cform);
        % Define thresholds for channel 1 based on histogram settings
        channel1Min = 0;
        channel1Max = 83.040;

        % Define thresholds for channel 2 based on histogram settings
        channel2Min = -40.254;
        channel2Max = 23.068;

        % Define thresholds for channel 3 based on histogram settings
        channel3Min = -54.402;
        channel3Max = -7.623;
        % Create mask based on chosen histogram thresholds
        BW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
            (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
            (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);
        % Initialize output masked image based on input image.
        BW = bwareaopen(BW, 300); % Parametro
        BW2 =  imdilate(BW,strel('square',10));
        BW2 = bwareaopen(BW2, 2000); % Parametro
        BW = BW & BW2;
        BW = imerode(BW,strel('square',5));
        BW = imdilate(BW,strel('square',5));
        BW = bwareaopen(BW, 300); % Parametro
        CC = bwconncomp(BW);
        for i=1:size(CC.PixelIdxList,2)
            BW_TMP = false(size(BW));
            BW_TMP(CC.PixelIdxList{i}) = 1;
            s = regionprops(BW_TMP,'BoundingBox');
            img = imcrop(BW_TMP,s.BoundingBox);
            white = sum(sum(img));
            black = size(img,1)*size(img,2);
            %figure, imshow(img);
            if(white<0.2*black)
               BW(CC.PixelIdxList{i}) = 0; 
            end
        end
        maskedRGBImage = RGB;
        % Set background pixels where BW is false to zero.
        maskedRGBImage(repmat(~BW,[1 1 3])) = 0;
        R = maskedRGBImage(:,:,1);
        G = maskedRGBImage(:,:,2);
        B = maskedRGBImage(:,:,3);
        R(BW) = 0;
        G(BW) = 0;
        B(BW) = 255;
        maskedRGBImage(:,:,1)=R;
        maskedRGBImage(:,:,2)=G;
        maskedRGBImage(:,:,3)=B;
        RGB = im2uint8(RGB);
    case 'yellow'        
        I = rgb2hsv(RGB);

        % Define thresholds for channel 1 based on histogram settings
        channel1Min = 0.120;
        channel1Max = 0.216;

        % Define thresholds for channel 2 based on histogram settings
        channel2Min = 0.405;
        channel2Max = 0.962;

        % Define thresholds for channel 3 based on histogram settings
        channel3Min = 0.000;
        channel3Max = 1.000;

        % Create mask based on chosen histogram thresholds
        BW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
            (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
            (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);
        BW = bwareaopen(BW, 300); % Parametro
        BW2 =  imdilate(BW,strel('square',10));
        BW2 = bwareaopen(BW2, 2000); % Parametro
        BW = BW & BW2;
        BW = imerode(BW,strel('square',5));
        BW = imdilate(BW,strel('square',5));
        BW = bwareaopen(BW, 300); % Parametro
        CC = bwconncomp(BW);
        for i=1:size(CC.PixelIdxList,2)
            BW_TMP = false(size(BW));
            BW_TMP(CC.PixelIdxList{i}) = 1;
            s = regionprops(BW_TMP,'BoundingBox');
            img = imcrop(BW_TMP,s.BoundingBox);
            white = sum(sum(img));
            black = size(img,1)*size(img,2);
            %figure, imshow(img);
            if(white<0.2*black)
               BW(CC.PixelIdxList{i}) = 0; 
            end
        end
        maskedRGBImage = RGB;
        % Set background pixels where BW is false to zero.
        maskedRGBImage(repmat(~BW,[1 1 3])) = 0;
        R = maskedRGBImage(:,:,1);
        G = maskedRGBImage(:,:,2);
        B = maskedRGBImage(:,:,3);
        R(BW) = 255;
        G(BW) = 255;
        B(BW) = 0;
        maskedRGBImage(:,:,1)=R;
        maskedRGBImage(:,:,2)=G;
        maskedRGBImage(:,:,3)=B;
    otherwise
        warning('Unexpected band type. No elements created.')
end
obj = objBlobs(band);

%% Separazione delle chessboard
split = sum(BW);
split = find(split~=0);
if size(split,2) == 0
    return;
end
x = []; x(1,1) = split(1); k = 1;
for i=2:size(split,2)
    if(abs(split(i-1)-split(i))>1)
        x(k,2) = split(i-1);
        x(k+1,1) = split(i);
        k = k+1;
    end
end
x(k,2) = split(end); y = [];
for i = 1:size(x,1)
    bw_cut = BW(:,x(i,1):x(i,2),:);
    bw_cut = sum(bw_cut');
    idx = find(bw_cut~=0);
    y = [y; min(idx),max(idx)];
end

%% Calcolo della chessboard complementare & background
inv_BW = false(size(maskedRGBImage,1),size(maskedRGBImage,2));
blackwhite_BW = false(size(maskedRGBImage,1),size(maskedRGBImage,2));
colors = [{''}];
ii = [];
for i = 1:size(x,1)
    mask = false(size(maskedRGBImage,1),size(maskedRGBImage,2));
    mask(y(i,1):y(i,2),x(i,1):x(i,2))=true;
    color_square =  BW&mask;
    CC = bwconncomp(color_square);
    if(size(CC.PixelIdxList,2) == 1)
        BW_TMP = false(size(blackwhite_BW));
        BW_TMP(CC.PixelIdxList{1}) = 1;
        bbs = cell2mat(struct2cell(regionprops(BW_TMP,'BoundingBox')));
        if(bbs(3)*bbs(4) < 3000)
            uno = maskedRGBImage(:,:,1);
            due = maskedRGBImage(:,:,2);
            tre = maskedRGBImage(:,:,3);
            BW(CC.PixelIdxList{1}) = 0;
            uno(CC.PixelIdxList{1}) = 0;
            due(CC.PixelIdxList{1}) = 0;
            tre(CC.PixelIdxList{1}) = 0;        
            maskedRGBImage = cat(3,uno,due,tre);
            ii = [ii, i];
            colors(i) = {''};
            continue;
        end
    end
    if(size(CC.PixelIdxList,2) == 2) % Risolvere il problema quando abbiamo solo due componenti connesse all'interno di una sezione. Cio� quando abbiamo solo i due quadrati intermedi, per calcolare la maschera complementare dobbiamo simulare quello di sotpra e quello di sottoaltrimenti la regione sarebbe ritagliata troppo piccola.
        stats = cell2mat(struct2cell(regionprops(color_square,'BoundingBox'))); % Devo essere sicuro siano due blob isolati
        % TODO Se la distanza tra due centroidi di blob � piccola rimuovili
%         color_square(CC.PixelIdxList{1}) = 0;
%         idx = find(color_square == 1);
%         [idx,idy]=ind2sub(size(mask),idx);
%         j = boundary(idx,idy,0.98); % Parametro
%         mask = poly2mask(idy(j),idx(j), size(mask,1), size(mask,2));
%         contour_square = edge(imfill(mask,'holes'));
%         p_filt = fspecial('prewit');
%         x_sides = imfilter(contour_square,p_filt);
%         y_sides = imfilter(contour_square,p_filt');
%         xy_sides = x_sides & y_sides;
%         x_sides = bwareaopen(x_sides & ~xy_sides, 5);
%         y_sides = bwareaopen(y_sides & ~xy_sides, 5);
        if((stats(4) + stats(8))< 0.5*size(maskedRGBImage,1))
            mask = false(size(maskedRGBImage,1),size(maskedRGBImage,2));
            %en = abs(y(i,2)-y(i,1))/3;
            en = min(stats(4),stats(8)); % Modificato perch� a volte dei due blob quello sopra/sotto � pi� lungo dato che ha connesso anche la tessera adiacente, allungo solo se non vado oltre l'immagine
            if(round(y(i,1)-en)>0)
                y(i,1) = round(y(i,1)-en);
            end
            if(round(y(i,2)+en)<size(maskedRGBImage,1))
                y(i,2) = round(y(i,2)+en);
            end
            mask(y(i,1):y(i,2),x(i,1):x(i,2))=true;
            color_square =  BW&mask;
            side = mean(x(i,1),x(i,2))<size(maskedRGBImage,2)*0.5;
            color_square(y(i,1):y(i,1)+3,x(i,1)+10*(~side):x(i,2)-10*side) = 255;
            color_square(y(i,2)-3:y(i,2),x(i,1)+10*side:x(i,2)-10*(~side)) = 255;
        end
    end
    idx = find(color_square == 1);
    [idx,idy]=ind2sub(size(color_square),idx);
    j = boundary(idx,idy,0.1); % Parametro
    color_square = poly2mask(idy(j),idx(j), size(color_square,1), size(color_square,2));   
    mask = ~BW&mask&color_square;
    mask = imerode(mask,strel('line',7,0));
    image = rgb2gray(RGB);
    square = image.*uint8(mask);
    square =  square(y(i,1):y(i,2),x(i,1):x(i,2));
    square_bw = im2bw(square,0.4); % Parametro
    square_bw = bwareaopen(square_bw, 400); % Parametro
    blackwhite_BW(y(i,1):y(i,2),x(i,1):x(i,2)) = blackwhite_BW(y(i,1):y(i,2),x(i,1):x(i,2)) + square_bw;
    if(sum(sum(bwconvhull(square_bw)))>0.2*size(square_bw,1)*size(square_bw,2))
        colors(i) = {'White'};
    else
        colors(i) = {'Black'};
    end
    inv_BW = inv_BW + mask;
end
y(ii,:) = [];
x(ii,:) = [];
colors(:,ii) = [];
obj.bbox_x = x;
obj.bbox_y = y;
obj.color_mask = BW;

%% Fix dei triangoli del background
CC = bwconncomp(blackwhite_BW);
for i = 1:size(CC.PixelIdxList,2)
    BW_TMP = false(size(blackwhite_BW));
    BW_TMP(CC.PixelIdxList{i}) = 1;
    s = regionprops(BW_TMP,'BoundingBox');
    img = imcrop(BW_TMP,s.BoundingBox); 
    toFlip = isTriangle(img);
    white = sum(sum(img));
    black = size(img,1)*size(img,2); 
    %figure, subplot(211); imshow(img);
    if(white<0.6*black && black>500 && black<10000 && toFlip)
        img_fill = img |rot90(img,2);
        img_fill = imfill(img_fill,'holes');
        %subplot(212); imshow(img_fill);        
        blackwhite_BW(s.BoundingBox(2):s.BoundingBox(2)+s.BoundingBox(4),s.BoundingBox(1):s.BoundingBox(1)+s.BoundingBox(3)) = img_fill; 
    end
end
obj.black_white = imfill(blackwhite_BW,'holes');

%% Separazione delle Tessere
BW = bwareaopen(BW, 50); % Parametro
inv_BW = bwareaopen(inv_BW, 50); % Parametro
BW_edge = edge(BW);
assert(sum(sum(inv_BW&BW))==0);
inv_BW_edge = edge(inv_BW);
eedge = imdilate(BW_edge|inv_BW_edge,strel('square',3));
eedge = imerode(eedge,strel('square',2));
BW(eedge) = 0;
inv_BW(eedge) = 0;
bwD = BW | inv_BW;
bwD = bwareaopen(bwD, 50); % Parametro

for i = 1:size(x,1)
    chess = bwD(y(i,1):y(i,2),x(i,1):x(i,2));
    condition = true;
    figure
    while condition
        CC = bwconncomp(chess);
        bb = regionprops(CC,'BoundingBox'); bboxCorr = zeros(1,size(bb,1));
        for iter=1:size(bb,1)
            BW_TMP = false(size(chess));
            BW_TMP(CC.PixelIdxList{iter}) = 1;
            s = regionprops(BW_TMP,'BoundingBox');
            img = imcrop(BW_TMP,s.BoundingBox);
            toCompare = bwconvhull(img);
            bboxCorr(iter) = corr2(toCompare,imfill(img,'holes'));
        end
        if(find(bboxCorr<0.72))
            idx = find(bboxCorr<0.72);
            tmp_chess = zeros(size(chess));
            for k=1:size(idx)
                tmp_chess(CC.PixelIdxList{idx(k)}) = 1;
                chess(CC.PixelIdxList{idx(k)}) = 0;
                %tmp_chess = imerode(tmp_chess,strel('square',2));
                tmp_chess = bwdist(~tmp_chess,'chessboard');
                tmp_chess = tmp_chess>1;
                imshowpair(tmp_chess,chess,'falsecolor');
                chess = chess | tmp_chess;
                chess = bwareaopen(chess, 100);
            end      
        else
            condition = false;
        end        
    end
    bwD(y(i,1):y(i,2),x(i,1):x(i,2))=chess;
end

L = bwlabel(bwD);
RGB = label2rgb(L);
imshow(RGB);

%% Fix dei triangoli della machera complementare
CC = bwconncomp(bwD);
for i = 1:size(CC.PixelIdxList,2)
    BW_TMP = false(size(blackwhite_BW));
    BW_TMP(CC.PixelIdxList{i}) = 1;
    s = regionprops(BW_TMP,'BoundingBox');
    img = imcrop(BW_TMP,s.BoundingBox);
    toFlip = isTriangle(img);
    white = sum(sum(img));
    black = size(img,1)*size(img,2);    
%    figure, subplot(211); imshow(img);
    if(white<0.6*black && black>500 && black<10000 && toFlip)
        img_fill = img |rot90(img,2);
        img_fill = imfill(img_fill,'holes');
%        subplot(212); imshow(img_fill);
        img_fill = img_fill(4:end-3,4:end-3);
        bwD((s.BoundingBox(2)+3):(s.BoundingBox(2)+s.BoundingBox(4)-3),(s.BoundingBox(1)+3):(s.BoundingBox(1)+s.BoundingBox(3)-3)) = img_fill; 
    end
end
inv_BW = bwD-BW;
obj.inv_color_mask = inv_BW;


for i = 1:size(x,1)
    chess = bwD(y(i,1):y(i,2),x(i,1):x(i,2));
    chess = bwareaopen(chess, 50); % Parametro
    figure, imshow(chess);
    CC = bwconncomp(chess);
    s = regionprops(CC,'centroid');
    obj = putCenters(obj, x(i,1), y(i,1), s, i, size(chess,1));
    obj.chess(i).mask = false(size(bwD));
    obj.chess(i).mask(y(i,1):y(i,2),x(i,1):x(i,2)) = chess;
    obj.chess(i).background = colors(i);
%     centroids = cat(1, s.Centroid);
%     %Display original image and superimpose centroids.
%     hold on
%     plot(centroids(:,1), centroids(:,2), 'r*')    
end
close all;
figure,imshow(BW); hold on;
for i = 1:size(x,1)
    X = obj.chess(i).center_x;
    Y = obj.chess(i).center_y;
    ind = find(obj.chess(i).center_x==0);
    [id_row,id_col]=ind2sub(size(X),ind);
    X(id_row,:) = 0;
    Y(id_row,:) = 0;
    axis = 1:5; axis(id_row) = [];
    for k=1:size(X,2)
        value = X(:,k); value(X(:,k)==0) = [];
        if(size(value,2) <=1)
            continue
        end
        X(:,k) = pchip(axis,value,1:5);
        value = Y(:,k); value(Y(:,k)==0) = [];
        if(size(value,2) <=1)
            continue
        end
        Y(:,k) = pchip(axis,value,1:5); 
    end
    if(size(find(X==0),1)>0 || size(find(Y==0),1)>0)
        continue
    end
    % trovo la posizione giusta degli zeri in funzione di quelle che sono
    % le tessere mancanti
    for k=1:size(id_row,1)
        row_Xchess = obj.chess(i).center_x(id_row(k),:);
        row_Ychess = obj.chess(i).center_y(id_row(k),:);
        row_Xchess(row_Xchess ==0) = [];
        row_Ychess(row_Ychess ==0) = [];
        row_X = X(id_row(k),:);
        index = 1:size(row_Xchess,2);
        for j=1:size(row_Xchess,2)
            curr_val_x = row_Xchess(j);
            curr_val_y = row_Ychess(j);
            row = abs(row_X-curr_val_x);
            idx = find(row== min(row));
            X(id_row(k),idx) = curr_val_x;
            Y(id_row(k),idx) = curr_val_y;
            index(index == idx) = [];
        end
        X(id_row(k),index) = 0;
        Y(id_row(k),index) = 0;
    end  
    obj.chess(i).center_x = X;
    obj.chess(i).center_y = Y;
    ind = find(obj.chess(i).center_x==0);
    while(size(ind,1)>0)
        %       mod(4,2)   % Even  number
        %       ans = 0
        %       mod(5,2)   % even number 
        %       ans =1
        [id_row,id_col]=ind2sub(size(obj.chess(i).center_x),ind(1)); %Prendo sempre il primo
        X = obj.chess(i).center_x(:,id_col);
        Y = obj.chess(i).center_y(id_row,:);
        xv = 1:size(X,1);
        yv = 1:size(Y,2);    
        X_pchip = X; X_pchip(X==0)=[];
        Y_pchip = Y; Y_pchip(Y==0)=[];
        if(size(X_pchip,1) == 1) % Se abbiamo solo un punto switcha sulle righe/colle
            X = obj.chess(i).center_x(id_row,:);
            xv = 1:size(X,1);   
            X_pchip = X; X_pchip(X==0)=[];
        end
        if(size(Y_pchip,2) == 1)
            Y = obj.chess(i).center_y(:,id_col);
            yv = 1:size(Y,1);    
            Y_pchip = Y; Y_pchip(Y==0)=[];
        end
        xv_pchip = xv; xv_pchip(X==0)=[];
        yv_pchip = yv; yv_pchip(Y==0)=[];
        if(size(yv_pchip,2) == 1 || size(xv_pchip,2) == 1) % Se quanto fatto prima non risolve cancella la colonna
            continue;
        end
        val_x =  pchip(xv_pchip,X_pchip',xv);
        o = obj.chess(i).center_x;
        obj.chess(i).center_x(id_row,id_col) = val_x(id_row);
        o =  obj.chess(i).center_y;
        val_y = pchip(yv_pchip,Y_pchip,yv);
        obj.chess(i).center_y(id_row,id_col) = val_y(id_col); 
        ind = find(obj.chess(i).center_x==0);
        disp(['fix',num2str(i)]);
    end
    for k=1:size(obj.chess(i).center_x,1)
        for j=1:size(obj.chess(i).center_x,2)
            scatter(obj.chess(i).center_x(k,j),obj.chess(i).center_y(k,j));
        end
    end 
end


% 
% 
% CC = bwconncomp(bwD);
% s = regionprops(CC,'centroid');
% figure,imshow(BW); hold on;
% for i = 1:size(s,1)
%     scatter(s(i).Centroid(1),s(i).Centroid(2)) % Riferimento Y verso il basso
% end
% 
% 
% 
% centers = zeros(5,5,2,size(x,1));
% for i = 1:size(x,1)
%     color = BW(y(i,1):y(i,2),x(i,1):x(i,2));
%     color = imclose(color,strel('disk', 3));
%     d = size(color,1); j=1; 
%     for k=floor(d/10):2*floor(d/10):d
%         centers(j,:,2,i) = k;
%         vline = color(k,:);
%         idx = find(vline~=0);
%         if(size(idx,2) ~= 0)
%             xx = []; xx(1,1) = idx(1); kk = 1;
%             for jj=2:size(idx,2)
%                 if(abs(idx(jj-1)-idx(jj))>1)
%                     xx(kk,2) = idx(jj-1);
%                     xx(kk+1,1) = idx(jj);
%                     kk = kk+1;
%                 end
%             end
%             xx(kk,2) = idx(end);
%             for jj=1:size(xx,1)
%                 centers(j,jj,1,i) = 0.5*(xx(jj,2)-xx(jj,1))+xx(jj,1);
%             end
%         end
%         j = j+1;
%     end
% end
% 
% for i = 1:size(x,1)
%     color = BW(y(i,1):y(i,2),x(i,1):x(i,2));
%     figure,imshow(color); hold on;
%     for k=1:5
%         for j=1:5
%             scatter(centers(k,j,1,i),centers(k,j,2,i)) % Riferimento Y verso il basso
%         end
%     end
% end
% 
% 
% for i = 1:size(x,1)
%     j = size(centers,2);
%     for k=1:size(centers,2)
%         x_axis = [centers(1,k,1,i), centers(3,k,1,i), centers(5,k,1,i)];
%         y_axis = [centers(1,k,2,i), centers(3,k,2,i), centers(5,k,2,i)];
%         if(centers(1,k,1,i) == 0 || centers(3,k,1,i) == 0 || centers(5,k,1,i) ==0)
%             centers(1,k,1,i) = 0;
%             centers(3,k,1,i) = 0;
%             centers(5,k,1,i) = 0;
%         end
%         if(sum(x_axis)>0 && k<=j)
%               centers(2,j,1,i) = (x_axis(1)*centers(2,k,2,i) - x_axis(3)*centers(2,k,2,i) - x_axis(1)*y_axis(3) + x_axis(3)*y_axis(1))/(y_axis(1) - y_axis(3));
%               centers(4,j,1,i) = (x_axis(1)*centers(4,k,2,i) - x_axis(3)*centers(4,k,2,i) - x_axis(1)*y_axis(3) + x_axis(3)*y_axis(1))/(y_axis(1) - y_axis(3));
%         end
%         x_axis = [centers(2,k,1,i), centers(4,k,1,i)];
%         y_axis = [centers(2,k,2,i), centers(4,k,2,i)];
%         if(centers(2,k,1,i) == 0 || centers(4,k,1,i) == 0)
%             centers(2,k,1,i) = 0;
%             centers(4,k,1,i) = 0;
%         end
%         if(sum(x_axis)>0  && k<=j)
%             centers(1,j,1,i) = (x_axis(1)*centers(1,k,2,i) - x_axis(2)*centers(1,k,2,i) - x_axis(1)*y_axis(2) + x_axis(2)*y_axis(1))/(y_axis(1) - y_axis(2));
%             centers(3,j,1,i) = (x_axis(1)*centers(3,k,2,i) - x_axis(2)*centers(3,k,2,i) - x_axis(1)*y_axis(2) + x_axis(2)*y_axis(1))/(y_axis(1) - y_axis(2));
%             centers(5,j,1,i) = (x_axis(1)*centers(5,k,2,i) - x_axis(2)*centers(5,k,2,i) - x_axis(1)*y_axis(2) + x_axis(2)*y_axis(1))/(y_axis(1) - y_axis(2));
%         end
%         j = j-1;
%     end
%      for k=1:size(centers,1)
%          centers(k,:,1,i) = sort(centers(k,:,1,i));
%      end
% end
% 
% % for i = 1:size(centers,4) % Immagini
% %     color = BW(y(i,1):y(i,2),x(i,1):x(i,2));
% %     color = imclose(color,strel('disk', 3));
% %     for j = 1:size(centers,1) %Riga
% %         for k = 1:size(centers,2) %Colonna
% %             %[centers(j,k,1,i), centers(j,k,2,i)] = updatePos(centers(j,k,1,i), centers(j,k,2,i));
% %             if(centers(k,j,1,i) == 0)
% %                 first = k+1;
% %             end
% %             if(j==1 && k ==first)
% %                 neighbour = ['right', 'down'];
% %             elseif(j==5 && k == 5)
% %                 neighbour = ['left', 'up'];
% %             elseif(j==1 && k == 5)
% %                 neighbour = ['down', 'left'];
% %             elseif(j==5 && k == first)
% %                 neighbour = ['right', 'up'];   
% %             elseif(j==5 && k>first && k<5)
% %                 neighbour = ['right', 'left', 'up'];  
% %             elseif(j>1 && j<5 && k == first)
% %                 neighbour = ['right', 'down', 'up'];   
% %             elseif(j==1 && k>first && k<5)
% %                 neighbour = ['right', 'down', 'left'];
% %             elseif(j>1 && j<5 && k == 5 )
% %                 neighbour = ['down', 'left', 'up'];
% %             end
% %             
% %             if(color(centers(j,k,2,i),centers(j,k,1,i))) % Se il pixel � colorato fai una cosa altrimenti un altra
% %                 
% %             else
% %                  
% %             end
% %         end
% %     end
% % end
% 
% for i = 1:size(x,1)
%     color = BW(y(i,1):y(i,2),x(i,1):x(i,2));
%     figure,imshow(color); hold on;
%     for k=1:5
%         for j=1:5
%             scatter(centers(k,j,1,i),centers(k,j,2,i)) % Riferimento Y verso il basso
%         end
%     end
% end
% 
% 
% % Riorganizzare i punti a partire dalla configurazione approssimativa
% 
% % trovare nei centri l'oggetto reale, spostarsi a sinistra e a destra
% % muovendosi sempre sul pi� vicino rispetto a quello vero
% % aggiornare con la posizione assoluta
% 
% for i = 1:size(x,1)
%     for k=1:5
%         for j=1:5
%             if(centers(k,j,1,i) ~= 0)
%                 centers(k,j,1,i) = x(i,1) + centers(k,j,1,i);
%                 centers(k,j,2,i) = y(i,1) + centers(k,j,2,i);
%             else
%                 centers(k,j,2,i) = 0;
%             end        
%         end
%     end
% end
% 
% figure,imshow(BW); hold on;
% for i = 1:size(x,1)
%     for k=1:5
%         for j=1:5
%             scatter(centers(k,j,1,i),centers(k,j,2,i)) % Riferimento Y verso il basso
%         end
%     end
% end

obj.masked_rgb = im2uint8(maskedRGBImage);
disp('ciao');
close all;





