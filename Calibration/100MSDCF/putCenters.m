function obj = putCenters(obj, xOffset, yOffset, centers, index, img_dim)
    arr= struct2cell(centers);
    obj.chess(index) = objChess();
    x = zeros(1,size(centers,1));
    y = zeros(1,size(centers,1));
    for i=1:size(centers,1)
        v = arr{i};
        x(i) = v(1);
        y(i) = v(2);
    end 
    th = floor(img_dim/10);
    j = 1;
    for k=floor(img_dim/10):2*floor(img_dim/10):img_dim
        idx = find(abs(y-k)<th);
        for i=1:size(idx,2)
            obj.chess(index).center_x(j,i) = xOffset+x(idx(i));
            obj.chess(index).center_y(j,i) = yOffset+y(idx(i));
        end
        j = j+1;
    end
end

