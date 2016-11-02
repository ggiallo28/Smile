function image = line2image(line,dim)
    image = false(dim);    
    coeffs = coeffvalues(line);
    x = 1:0.001:dim(2);
    y = floor(polyval(coeffs,x));
    x = floor(x);
    x(y<1 | y>dim(1)) = [];
    y(y<1 | y>dim(1)) = [];
    for j = 1:size(x,2)
        image(y(j),x(j)) = 1;
    end
end

