function [x, y, BW, maskedRGBImage] = img_filtering(x, y, BW, maskedRGBImage, op_th)
    y_mean = mean(abs(y(:,1)-y(:,2)));
    % x_min = min(abs(x(:,1)-x(:,2)));
    % Threshold = round(0.2*y_mean*0.5*x_min);
    Threshold = round(5*y_mean);
    ii = [];
    for i = 1:size(x,1)
        mask = false(size(maskedRGBImage,1),size(maskedRGBImage,2));
        mask(y(i,1):y(i,2),x(i,1):x(i,2))=true;
        color_square =  BW&mask;
        CC = bwconncomp(color_square);
        if(size(CC.PixelIdxList,2) == 1)
            if(abs(y(i,1)-y(i,2)) < 0.5*y_mean)
                uno = maskedRGBImage(:,:,1);
                due = maskedRGBImage(:,:,2);
                tre = maskedRGBImage(:,:,3);
                BW(CC.PixelIdxList{1}) = 0;
                uno(CC.PixelIdxList{1}) = 0;
                due(CC.PixelIdxList{1}) = 0;
                tre(CC.PixelIdxList{1}) = 0;        
                maskedRGBImage = cat(3,uno,due,tre);
                ii = [ii, i];
                continue;
            end
        end
    end
    y(ii,:) = [];
    x(ii,:) = [];
end