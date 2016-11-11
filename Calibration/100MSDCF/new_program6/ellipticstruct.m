% CIRCULARSTRUCT
%
% Function to construct a circular structuring element
% for morphological operations.
%
% function strel = circularstruct(radius)
%
% Note radius can be a floating point value though the resulting
% circle will be a discrete approximation
%
% Peter Kovesi   March 2000

function strel = ellipticstruct(a,b)

if a < 1 | b<1
  error('radius must be >= 1');
end

x_axis = ceil(2*a);  % Diameter of structuring element
y_axis = ceil(2*b);  % Diameter of structuring element

if mod(x_axis,2) == 0     % If diameter is a odd value
 x_axis = x_axis + 1;        % add 1 to generate a `centre pixel'
end

if mod(y_axis,2) == 0     % If diameter is a odd value
 y_axis = y_axis + 1;        % add 1 to generate a `centre pixel'
end

ra = fix(x_axis/2);
rb = fix(y_axis/2);
strel = zeros((rb*2+1),(ra*2+1));
for i =1:size(strel,1)
    for j =1:size(strel,2)
        xc = ra+1;
        yc = rb+1;
        el1 = (((j-xc)^2)/((ra)^2));
        el2 = (((i-yc)^2)/((rb)^2));
        if( el1 + el2 <=1 )
            strel(i,j) = 1;
        end
    end
end

