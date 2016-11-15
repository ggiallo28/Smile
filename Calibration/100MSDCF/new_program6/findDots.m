function points = findDots(Igray)
    Isin = im2single(Igray);
    % Provo la detect a differenti risoluzioni
    sigma = 2; 
    peakThreshold = 0.10; % threshold for corner metric

    G = fspecial('gaussian', round(sigma * 7)+1, sigma);
    Ig = imfilter(Isin, G, 'conv');

    derivFilter = [-1 0 1];

    % first derivatives
    Iy = imfilter(Ig, derivFilter', 'conv');
    Ix = imfilter(Ig, derivFilter, 'conv');

    % first derivative at 45 degrees
    I_45 = Ix * cos(pi/4) + Iy * sin(pi/4);
    I_n45 = Ix * cos(-pi/4) + Iy * sin(-pi/4);

    % second derivative
    Ixy = imfilter(Ix, derivFilter', 'conv');

    I_45_x = imfilter(I_45, derivFilter, 'conv');
    I_45_y = imfilter(I_45, derivFilter', 'conv');    

    I_45_45 = I_45_x * cos(-pi/4) + I_45_y * sin(-pi/4);

    % suppress the outer corners
    cxy = sigma^2 * abs(Ixy) - sigma * (abs(I_45) + abs(I_n45));
    cxy(cxy < 0) = 0;
    c45 = sigma^2 * abs(I_45_45) - sigma * (abs(Ix) + abs(Iy));
    c45(c45 < 0) = 0;


    Ix2 = Ix .^ 2;
    Iy2 = Iy .^ 2;
    Ixy = Ix .* Iy;

    G = fspecial('gaussian', 2, 15);

    Ix2 = imfilter(Ix2, G);
    Iy2 = imfilter(Iy2, G);
    Ixy = imfilter(Ixy, G);
    
    points = find_peaks(cxy, peakThreshold);
    scores0 = cxy(sub2ind(size(cxy), points(:, 2), points(:, 1)));
    if(size(points,1)>1)
        [~, i] = max(scores0);
        points = points(i,:);
    end
end

