function [color_square, y_c, isDone] = checkbadthings(color_square, Threshold, op_th, y, x, transtions, obj_idx, chess_idx, cond)
    ymean = mean(abs(y(:,1)-y(:,2)));
    ratio = abs(x(1)-x(2))/ymean;
    CC = bwconncomp(color_square);
    y_c = [];
    idx = find(transtions(:,1)==obj_idx & transtions(:,2)==chess_idx, 1);
    direction = transtions(idx,3);
    % Quando non devo farlo.
    % Sono due e sono buoni, sono 3 e sono buoni, è uno ed è alto giusto è
    % se è nell'intersezione?
    isDone = ~isempty(idx);
%    isDone = ratio < 0.2 && (CC.NumObjects < 3 || cond) && CC.NumObjects ~=0;
    % Se ho esattamente due blob erodi solo se il ratio non è minore di una
    % certa soglia, quindi non è un quadrato
    if(CC.NumObjects == 2)
        Images = regionprops(CC,'Image','BoundingBox');
        c1 = Images(1).Image;
        c2 = Images(2).Image;
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

    if(isDone)
        color_square = imopen(color_square,strel('rectangle',[2,op_th]));
        color_square = bwareaopen(color_square,Threshold);
        split = sum(color_square);
        split = find(split~=0);
        [~, y_c] = compute_projections(color_square, split);
        disp('its bad');
    else
        disp('skipped');
    end
end