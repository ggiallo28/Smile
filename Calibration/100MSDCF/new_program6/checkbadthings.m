function [color_square, y_c, x_c, isDone, isErased] = checkbadthings(color_square, Threshold, y, x, transtions, obj_idx, chess_idx)
    ymean = mean(abs(y(:,1)-y(:,2)));
    ratio = abs(x(1)-x(2))/ymean;
    CC = bwconncomp(color_square);
    y_c = []; x_c = [];
    idx = find(transtions(:,1)==obj_idx & transtions(:,2)==chess_idx, 1);
    direction = transtions(idx,3);
    % Quando non devo farlo.
    % Sono due e sono buoni, sono 3 e sono buoni, è uno ed è alto giusto è
    % se è nell'intersezione?
    isDone = ~isempty(idx);
    isErased = sum(sum(color_square)) == 0;
%    isDone = ratio < 0.2 && (CC.NumObjects < 3 || cond) && CC.NumObjects ~=0;
    % Se ho esattamente due blob erodi solo se il ratio non è minore di una
    % certa soglia, quindi non è un quadrato
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
        isDone = isDone && (ratio1<0.7 || ratio2<0.7 || ratio3<0.7 || ~(hh<ymean*1.05 && hh>ymean*0.95)); % Se ho 3 quadrati giusti l'artezza calcolata dovrà essere simile a quella media
    end
    
    if(isDone)
        backup_color_square = color_square;
        [~, ~, ~, ~, rct_before]  = getRotatedRectangle(backup_color_square);
        max_size = max(transtions(:,4)); side = sqrt(max_size);
        op_th = round(0.15*(transtions(idx,4)/side));
        color_square = imerode(color_square,strel('rectangle',[op_th,round(op_th*1.2)])); %0, 1.2
        color_square = bwareaopen(color_square,round(Threshold*1.3)); % 0, 1.15, 1.3
        color_square = imdilate(color_square,strel('rectangle',[op_th,round(op_th*1.2)]));
        [~, ~, ~, ~, rct_after]  = getRotatedRectangle(backup_color_square);
        if 0.95*rct_before.size(1)< rct_after.size(1) && 1.05*rct_before.size(1) > rct_after.size(1)
            color_square = backup_color_square;
            isDone = false;
            disp('restored');
        elseif ~sum(sum(color_square)) == 0
            split = medfilt2(sum(color_square),[1,10]);
            split = find(split~=0);
            [x_c, y_c] = compute_projections(color_square, split);
            disp('its bad');
        else
            disp('erased')
            isErased = true;
        end
    else
        disp('skipped');
    end
end