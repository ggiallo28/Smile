function image = line2image(line,dim)
    image = false(dim);    
    coeffs = coeffvalues(line);
    y1 = 1;
    y2 = dim(1);
    x1 = (y1-coeffs(2))/(coeffs(1));
    x2 = (y2-coeffs(2))/(coeffs(1));
    if(x1<1 || x1>dim(2) || x2<1 || x2>dim(2))
        x1 = 1;
        x2 = dim(2);
        y1 = line(x1);
        y2 = line(x2);
    end
    image = logical(rgb2gray(insertShape(double(image),'Line',[x1 y1 x2 y2])));
%     
%     x = 1:0.001:dim(2);
%     y = floor(polyval(coeffs,x));
%     x = floor(x);
%     x(y<1 | y>dim(1)) = [];
%     y(y<1 | y>dim(1)) = [];
%     for j = 1:size(x,2)
%         image(y(j),x(j)) = 1;
%     end
end

