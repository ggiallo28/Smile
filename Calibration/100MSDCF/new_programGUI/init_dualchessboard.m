function obj = init_dualchessboard(x,y,BW,RGB,obj,dim)
    for i = 1:size(x,1)
        mask = false(dim);
        mask(y(i,1):y(i,2),x(i,1):x(i,2))=true;
        color_square =  BW&mask;
        color_square = bwconvhull(color_square);
        mask = ~BW&mask&color_square;
        mask = imerode(mask,strel('line',7,0));
        image = rgb2gray(RGB);
        square = image.*uint8(mask);
        square =  square(y(i,1):y(i,2),x(i,1):x(i,2));
        square_bw = im2bw(square,0.4); % Parametro
        square_bw = bwareaopen(square_bw, 400); % Parametro
        if(sum(sum(square_bw))>0.5*sum(sum(mask)))%0.2*size(square_bw,1)*size(square_bw,2))
            obj.chess(i).background = 'White';
        else
            obj.chess(i).background = 'Black';
        end
    end
end