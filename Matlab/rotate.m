function [X,Y] = rotate(x, y, centerX, centerY, alpha)
R = rotz(alpha);
dimX = size(x);
dimY = size(y);
assert(dimX(1) == dimY(1));
assert(dimX(2) == dimX(2));
N = max(dimX(1),dimX(2)); % E' indifferente se supera l'assert di prima.
X = size(dimX(1),dimX(2));
Y = size(dimX(1),dimX(2));
for i=1:N
    vect = [x(i)-centerX;y(i)-centerY;0];
    res = R*vect;
    X(i) = res(1)+centerX;
    Y(i) = res(2)+centerY;
end

