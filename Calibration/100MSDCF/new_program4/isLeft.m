function tf = isLeft(line,x,y)
    coef = coeffvalues(line);
    xx = (y-coef(2))/coef(1);
    tf = xx>x;
end