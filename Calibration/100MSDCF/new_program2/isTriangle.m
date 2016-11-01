function bool = isTriangle(image)
    triangle = zeros(size(image));
    x1 = 0;
    y1 = 0;
    x2 = size(image,2);
    y2 = size(image,1);
    for i=1:size(image,1)
        for j=1:size(image,1)
            y = (j*y1 - j*y2 + x1*y2 - x2*y1)/(x1 - x2);
            if(i>y)
                triangle(i,j) = true;
            end
        end
    end
    uno = triangle;
    due = fliplr(uno);
    tre = flipud(due);
    quattro = fliplr(tre);
    arr = abs([corr2(image,uno),corr2(image,due),corr2(image,tre),corr2(image,quattro)]);
    c = max(arr);
    bool = c>0.6;
end