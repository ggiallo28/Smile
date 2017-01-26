function [maskedRGBImage, BW] = createColorMask(RGB, bw, output_minima_mid, Threshold, hist_size, band)
    switch band
    %% RED
        case 'red'
            RGB = im2double(RGB);
            hsv_inImg = rgb2hsv(RGB);
            red_th_min = output_minima_mid(1)/hist_size;
            red_th_max = output_minima_mid(4)/hist_size;
            bw_red = (hsv_inImg(:,:,1)<=red_th_min | hsv_inImg(:,:,1)>=red_th_max) & bw;
            bw_red = bwareaopen(bw_red,Threshold);
%            figure, imshow(RGB.*repmat(bw_red,[1 1 3]));
            BW = bwareaopen(bw_red, 300); % Parametro
            BW2 =  imdilate(BW,strel('square',10));
            BW2 = bwareaopen(BW2, 2000); % Parametro
            BW = BW & BW2;
            BW = imerode(BW,strel('square',5));
            BW = imdilate(BW,strel('square',5));
            BW = bwareaopen(BW, 300); % Parametro
            CC = bwconncomp(BW);
            for i=1:size(CC.PixelIdxList,2)
                BW_TMP = false(size(BW));
                BW_TMP(CC.PixelIdxList{i}) = 1;
                s = regionprops(BW_TMP,'BoundingBox');
                img = imcrop(BW_TMP,s.BoundingBox);
                white = sum(sum(img));
                black = size(img,1)*size(img,2);
                %figure, imshow(img);
                if(white<0.3*black) % Portato a 0.3 prima era 0.2
                   BW(CC.PixelIdxList{i}) = 0; 
                end
            end
            maskedRGBImage = RGB;
            % Set background pixels where BW is false to zero.
            maskedRGBImage(repmat(~BW,[1 1 3])) = 0;
            R = maskedRGBImage(:,:,1);
            G = maskedRGBImage(:,:,2);
            B = maskedRGBImage(:,:,3);
            hold off
            R(BW) = 255;
            G(BW) = 0;
            B(BW) = 0;
            maskedRGBImage(:,:,1)=R;
            maskedRGBImage(:,:,2)=G;
            maskedRGBImage(:,:,3)=B;
            RGB = im2uint8(RGB);
    %% GREEN
        case 'green'
            RGB = im2double(RGB);
            hsv_inImg = rgb2hsv(RGB);
            green_th_min = output_minima_mid(2)/hist_size;
            green_th_max = output_minima_mid(3)/hist_size;
            bw_green = hsv_inImg(:,:,1)>=green_th_min & hsv_inImg(:,:,1)<=green_th_max & bw;
            bw_green = bwareaopen(bw_green,Threshold);
%            figure, imshow(RGB.*repmat(bw_green,[1 1 3]));
            BW = bwareaopen(bw_green, 300); % Parametro
            BW2 =  imdilate(BW,strel('square',10));
            BW2 = bwareaopen(BW2, 2000); % Parametro
            BW = BW & BW2;
            BW = imerode(BW,strel('square',5));
            BW = imdilate(BW,strel('square',5));
            BW = bwareaopen(BW, 300); % Parametro
            CC = bwconncomp(BW);
            for i=1:size(CC.PixelIdxList,2)
                BW_TMP = false(size(BW));
                BW_TMP(CC.PixelIdxList{i}) = 1;
                s = regionprops(BW_TMP,'BoundingBox');
                img = imcrop(BW_TMP,s.BoundingBox);
                white = sum(sum(img));
                black = size(img,1)*size(img,2);
                %figure, imshow(img);
                if(white<0.2*black)
                   BW(CC.PixelIdxList{i}) = 0; 
                end
            end
            maskedRGBImage = RGB;
            % Set background pixels where BW is false to zero.
            maskedRGBImage(repmat(~BW,[1 1 3])) = 0;
            R = maskedRGBImage(:,:,1);
            G = maskedRGBImage(:,:,2);
            B = maskedRGBImage(:,:,3);
            R(BW) = 0;
            G(BW) = 255;
            B(BW) = 0;
            maskedRGBImage(:,:,1)=R;
            maskedRGBImage(:,:,2)=G;
            maskedRGBImage(:,:,3)=B;
            RGB = im2uint8(RGB);
    %% BLUE
        case 'blue'
            RGB = im2double(RGB);
            hsv_inImg = rgb2hsv(RGB);
            blu_th_min = output_minima_mid(3)/hist_size;
            blu_th_max = output_minima_mid(4)/hist_size;
            bw_blu = hsv_inImg(:,:,1)>=blu_th_min & hsv_inImg(:,:,1)<=blu_th_max & bw;
            bw_blu = bwareaopen(bw_blu,Threshold);
            bw_blu = imclose(bw_blu,strel('square',10));
%            figure, imshow(RGB.*repmat(bw_blu,[1 1 3]));
            % Initialize output masked image based on input image.
            BW = bwareaopen(bw_blu, 300); % Parametro
            BW2 =  imdilate(BW,strel('square',10));
            BW2 = bwareaopen(BW2, 2000); % Parametro
            BW = BW & BW2;
            BW = imerode(BW,strel('square',5));
            BW = imdilate(BW,strel('square',5));
            BW = bwareaopen(BW, 300); % Parametro
            CC = bwconncomp(BW);
            for i=1:size(CC.PixelIdxList,2)
                BW_TMP = false(size(BW));
                BW_TMP(CC.PixelIdxList{i}) = 1;
                s = regionprops(BW_TMP,'BoundingBox');
                img = imcrop(BW_TMP,s.BoundingBox);
                white = sum(sum(img));
                black = size(img,1)*size(img,2);
                %figure, imshow(img);
                if(white<0.2*black)
                   BW(CC.PixelIdxList{i}) = 0; 
                end
            end
            maskedRGBImage = RGB;
            % Set background pixels where BW is false to zero.
            maskedRGBImage(repmat(~BW,[1 1 3])) = 0;
            R = maskedRGBImage(:,:,1);
            G = maskedRGBImage(:,:,2);
            B = maskedRGBImage(:,:,3);
            R(BW) = 0;
            G(BW) = 0;
            B(BW) = 255;
            maskedRGBImage(:,:,1)=R;
            maskedRGBImage(:,:,2)=G;
            maskedRGBImage(:,:,3)=B;
            RGB = im2uint8(RGB);
    %% YELLOW
        case 'yellow' 
            RGB = im2double(RGB);
            hsv_inImg = rgb2hsv(RGB);
            yellow_th_min = output_minima_mid(1)/hist_size;
            yellow_th_max = output_minima_mid(2)/hist_size;
            bw_yellow = hsv_inImg(:,:,1)>=yellow_th_min & hsv_inImg(:,:,1)<=yellow_th_max & bw;
            bw_yellow = bwareaopen(bw_yellow,Threshold);
%            figure, imshow(RGB.*repmat(bw_yellow,[1 1 3]));
            BW = bwareaopen(bw_yellow, 300); % Parametro
            BW2 =  imdilate(BW,strel('square',10));
            BW2 = bwareaopen(BW2, 2000); % Parametro
            BW = BW & BW2;
            BW = imerode(BW,strel('square',5));
            BW = imdilate(BW,strel('square',5));
            BW = bwareaopen(BW, 300); % Parametro
            CC = bwconncomp(BW);
            for i=1:size(CC.PixelIdxList,2)
                BW_TMP = false(size(BW));
                BW_TMP(CC.PixelIdxList{i}) = 1;
                s = regionprops(BW_TMP,'BoundingBox');
                img = imcrop(BW_TMP,s.BoundingBox);
                white = sum(sum(img));
                black = size(img,1)*size(img,2);
                %figure, imshow(img);
                if(white<0.2*black)
                   BW(CC.PixelIdxList{i}) = 0; 
                end
            end
            maskedRGBImage = RGB;
            % Set background pixels where BW is false to zero.
            maskedRGBImage(repmat(~BW,[1 1 3])) = 0;
            R = maskedRGBImage(:,:,1);
            G = maskedRGBImage(:,:,2);
            B = maskedRGBImage(:,:,3);
            R(BW) = 255;
            G(BW) = 255;
            B(BW) = 0;
            maskedRGBImage(:,:,1)=R;
            maskedRGBImage(:,:,2)=G;
            maskedRGBImage(:,:,3)=B;
            RGB = im2uint8(RGB);
        otherwise
            warning('Unexpected band type. No elements created.')
    end
end