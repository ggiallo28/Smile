function obj = putCenters(obj, xOffset, yOffset, centers, index, img_dim)
    arr= struct2cell(centers);
    x = zeros(1,size(centers,1));
    y = zeros(1,size(centers,1));
    offset = 0;
    figure,imshow(obj.color_mask); hold on;
    while(size(x,2)>0 && size(y,2) >0)
        for i=1:size(centers,1)
            v = arr{i};
            x(i) = v(1);
            y(i) = v(2);
        end 
        obj.chess(index).centroid = [xOffset+mean(x), yOffset+mean(y)];
        th = floor(img_dim/10);  
        v = [];
        for k=floor(img_dim/10):2*floor(img_dim/10):img_dim
            idx = find(abs(y-k)<th);
            v = [v;size(idx,2)]; % Conto quante tessere ci sono per ogni riga
        end
        j = 1;
        for k=floor(img_dim/10):2*floor(img_dim/10):img_dim
            idx = find(abs(y-k)<th);
            if size(idx,2)>mode(v)+offset
                idx = idx(1:mode(v)+offset);
            end
            for i=1:size(idx,2)
                obj.chess(index).center_x(j,i) = xOffset+x(idx(i));
                obj.chess(index).center_y(j,i) = yOffset+y(idx(i));
                hold on; scatter(obj.chess(index).center_x(j,i),obj.chess(index).center_y(j,i))
            end
            x(idx) = [];
            y(idx) = [];
            j = j+1;
        end
        offset = offset + 1;
    end
end

