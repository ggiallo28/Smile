function [obj] = computeChess(Container, bw, hist_size, output_minima_mid, band)
%% INIT
    warning off
    RGB = Container.I;
    Threshold = Container.Threshold;
%    op_th = Container.op_th;
%         y_mean = mean(abs(y(:,1)-y(:,2)));
%     % x_min = min(abs(x(:,1)-x(:,2)));
%     % Threshold = round(0.2*y_mean*0.5*x_min);
%     Threshold = round(5*y_mean);
%% LOGIC
    [maskedRGBImage, BW] = createColorMask(RGB, bw, output_minima_mid, Threshold, hist_size, band);
    obj = objBlobs(band);
%% Separazione delle chessboard
    BW = bwareaopen(BW,Threshold);
    split = medfilt2(sum(BW),[1,10]);
    split = find(split~=0);
    if size(split,2) == 0
        obj.isEmpty = true;
        return;
    end
    [x, y] = compute_projections(BW, split);
%% Checkerboard Filtering
    [x, y, BW, maskedRGBImage] = img_filtering(x, y, BW, maskedRGBImage, Container.size_square, Threshold);
    figure, imshow(im2double(Container.I).*repmat(BW,[1 1 3]));
    obj.bbox_x = x;
    obj.bbox_y = y;
%% RAW Elements for separation
    obj.raw_bw = BW;
    obj.raw_maskedRGBImage = maskedRGBImage;
    for i=1:size(x,1)
       obj.chess(i) = objChess();
       obj.chess(i).raw_centroid = [mean(x(i,:)), mean(y(i,:))];
    end
%% Inizializza
    obj = init_dualchessboard(x,y,BW,RGB,obj,size(BW));
    obj.isEmpty = false;




