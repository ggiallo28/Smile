function [color_square, y_c, x_c, isDone, isErased] = checkbadthings(color_square, Threshold, y, x, transtions, obj_idx, chess_idx, Container, obj_chess)
% Usare direction per prelevare dei quadrati a sinistra e a destra per controllare l'altezza del blob corrente
% Erodere se il bordo di sopra è troppo piccolo, la fit avrà problemi
%% INIT
    idx = find(transtions(:,1)==obj_idx & transtions(:,2)==chess_idx, 1);
    isDone = ~isempty(idx); y_c = []; x_c = []; isErased = sum(sum(color_square)) == 0;
    if ~isDone
        return;
    end
    disp(['Analysing ', num2str(chess_idx),'...']);
    ymean = mean(abs(y(:,1)-y(:,2)));
%    ratio = abs(x(1)-x(2))/ymean;
%% GET LEFT/RIGHT SQUARES USING DIRECTIONS
    direction = transtions(idx,3);
   	ido = find(Container.raw_order(1,:)==obj_idx & Container.raw_order(2,:)==chess_idx) - direction;
    chess = obj_chess(Container.raw_order(1,ido)).raw_bw; cut_chess = false(size(chess));
    bbx = obj_chess(Container.raw_order(1,ido)).bbox_x(Container.raw_order(2,ido),:);
    bby = obj_chess(Container.raw_order(1,ido)).bbox_y(Container.raw_order(2,ido),:);
    cut_chess(bby(1):bby(2),bbx(1):bbx(2)) = chess(bby(1):bby(2),bbx(1):bbx(2));
    [lineTop, lineBot] = getTopBotLines(cut_chess);
    if direction == -1
        hindirection = abs(lineTop(bbx(1))-lineBot(bbx(1)));
    else
        hindirection = abs(lineTop(bbx(2))-lineBot(bbx(2)));
    end
    [idr,idc] = ind2sub(size(color_square),find(imfill(color_square,'holes')==1));
    rct = cv.minAreaRect([idr,idc]);
    if max(rct.size)<hindirection*1.05 && max(rct.size)>hindirection*0.95
         isDone = false; disp('tall');
         return;  
    end
%% 
    
    % Quando non devo farlo.
    % Sono due e sono buoni, sono 3 e sono buoni, è uno ed è alto giusto è
    % se è nell'intersezione?
   
    
%    isDone = ratio < 0.2 && (CC.NumObjects < 3 || cond) && CC.NumObjects ~=0;
    % Se ho esattamente due blob erodi solo se il ratio non è minore di una
    % certa soglia, quindi non è un quadrato
    
    CC = bwconncomp(color_square);
    if(CC.NumObjects == 2)
        Images = regionprops(CC,'Image','BoundingBox');
        c1 = Images(1).Image;
        c2 = Images(2).Image;
        c1 = prepareHull(c1);
        c2 = prepareHull(c2);
        ch1 = bwconvhull(c1);
        ch2 = bwconvhull(c2);
        w1 = sum(sum(c1));
        w2 = sum(sum(c2));
        wh1 = sum(sum(ch1)); 
        wh2 = sum(sum(ch2));
        ratio1 = w1/wh1;
        ratio2 = w2/wh2;
        isDone = isDone && (ratio1<0.7 || ratio2<0.7); %Provare ad erodere uno per volta
    end
    if(CC.NumObjects == 3)
        Images = regionprops(CC,'Image','BoundingBox','Centroid');
        c1 = Images(1).Image;
        c2 = Images(2).Image;
        c3 = Images(3).Image;
        c1 = prepareHull(c1);
        c2 = prepareHull(c2);
        c3 = prepareHull(c3);
        h1 = sqrt((Images(1).Centroid(1)-Images(2).Centroid(1)).^2 + (Images(1).Centroid(2)-Images(2).Centroid(2)).^2);
        h2 = sqrt((Images(3).Centroid(1)-Images(2).Centroid(1)).^2 + (Images(3).Centroid(2)-Images(2).Centroid(2)).^2);
        h3 = (h1+h2)/4;
        hh = h1+h2+h3;
        ch1 = bwconvhull(c1);
        ch2 = bwconvhull(c2);
        ch3 = bwconvhull(c3);
        w1 = sum(sum(c1));
        w2 = sum(sum(c2));
        w3 = sum(sum(c3));
        wh1 = sum(sum(ch1)); 
        wh2 = sum(sum(ch2));
        wh3 = sum(sum(ch3));
        ratio1 = w1/wh1;
        ratio2 = w2/wh2;
        ratio3 = w3/wh3;
        isDone = isDone && (ratio1<0.7 || ratio2<0.7 || ratio3<0.7 || ~(hh<hindirection*1.05 && hh>hindirection*0.95)); % Se ho 3 quadrati giusti l'artezza calcolata dovrà essere simile a quella media
    end
    
    if(isDone)
        backup_color_square = color_square;
        [~, ~, ~, ~, rct_before]  = getRotatedRectangle(backup_color_square);
        max_size = max(transtions(:,4)); side = sqrt(max_size);
        op_th = round(0.15*(transtions(idx,4)/side)); % 0.15, 0.1
        color_square = imerode(color_square,strel('rectangle',[round(op_th*1.5),op_th])); %0, 1.2
        color_square = bwareaopen(color_square,round(Threshold*1.3)); % 0, 1.15, 1.3
        color_square = imdilate(color_square,strel('rectangle',[round(op_th*1.5),op_th]));
        if sum(sum(color_square)) == 0
            disp('erased')
            isErased = true;
        else
            [~, ~, ~, ~, rct_after]  = getRotatedRectangle(color_square);
            if 0.95*rct_before.size(1) < rct_after.size(1) && 1.05*rct_before.size(1) > rct_after.size(1)
                color_square = backup_color_square;
                isDone = false;
                disp('restored');
            else
                split = medfilt2(sum(color_square),[1,10]);
                split = find(split~=0);
                [x_c, y_c] = compute_projections(color_square, split);
                disp('its bad');
            end
        end
    else
        disp('skipped');
    end
end