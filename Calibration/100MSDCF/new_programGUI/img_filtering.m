function [x, y, BW, maskedRGBImage] = img_filtering(x, y, BW, maskedRGBImage, size_square, Threshold)
    ii = [];
    y_mean = mean(abs(y(:,1)-y(:,2)));
    y_mean = mean([y_mean, sqrt(size_square)*5]);
    
    for i = 1:size(x,1)
        mask = false(size(maskedRGBImage,1),size(maskedRGBImage,2));
        mask(y(i,1):y(i,2),x(i,1):x(i,2))=true;
        color_square =  BW&mask;
        color_square =  imclose(color_square,strel('Line',round(0.1*sqrt(size_square)),90));
        CC = bwconncomp(color_square);
        if(size(CC.PixelIdxList,2) == 1 || sum(sum(color_square)) < 2*Threshold)
            if(abs(y(i,1)-y(i,2)) < 0.51*y_mean || sum(sum(color_square)) < 2*Threshold)
                if( ~isempty(CC.PixelIdxList) )
                    uno = maskedRGBImage(:,:,1);
                    due = maskedRGBImage(:,:,2);
                    tre = maskedRGBImage(:,:,3);

                    BW(CC.PixelIdxList{1}) = 0;
                    uno(CC.PixelIdxList{1}) = 0;
                    due(CC.PixelIdxList{1}) = 0;
                    tre(CC.PixelIdxList{1}) = 0;        
                    maskedRGBImage = cat(3,uno,due,tre);
                end
                ii = [ii, i];
            end
        end
    end
    y(ii,:) = [];
    x(ii,:) = [];
end