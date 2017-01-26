function [inv_BW, blackwhite_BW, colors, obj] = compute_dualchessboard(x,y,BW,RGB,obj,dim, Container)
    inv_BW = false(dim);
    blackwhite_BW = false(dim);
    colors = [{''}];
    for i = 1:size(x,1)
        mask = false(dim);
        mask(y(i,1):y(i,2),x(i,1):x(i,2))=true;
        color_square =  BW&mask;
        [color_square, obj.chess(i)] = fitSquare(color_square,x(i,1),x(i,2),obj.chess(i), Container);
        mask = ~BW&mask&color_square;
        mask = imerode(mask,strel('line',7,0));
        image = rgb2gray(RGB);
        square = image.*uint8(mask);
        square =  square(y(i,1):y(i,2),x(i,1):x(i,2));
        square_bw = im2bw(square,0.4); % Parametro
        square_bw = bwareaopen(square_bw, 400); % Parametro
        if(sum(sum(square_bw))>0.5*sum(sum(mask)))%0.2*size(square_bw,1)*size(square_bw,2))
            colors(i) = {'White'};
            blackwhite_BW(y(i,1):y(i,2),x(i,1):x(i,2)) = blackwhite_BW(y(i,1):y(i,2),x(i,1):x(i,2)) + square_bw;
        else
            colors(i) = {'Black'};
        end
        inv_BW = inv_BW + mask;
    end
    inv_BW = bwareaopen(inv_BW,100);
end