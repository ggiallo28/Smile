function [bw, inImg] = subtract_background(I, I_BG, confidence, Threshold)
    lab_image_bg = rgb2lab(I_BG);
    lab_imgae_fg = rgb2lab(I);
    a_band = lab_image_bg(:,:,2);
    mu_a = mean2(a_band);
    si_a = std2(a_band); 
    ll1 = a_band(:); ll1(ll1<=mu_a) = [];
    ll2 = a_band(:); ll2(ll2>=mu_a) = [];
    % Controlliamo quanti pixel ci sono a sinistra e a destra della media per aggiustare la soglia
    p = size(ll1,1)/size(ll2,1);

    channel2Min = mu_a-confidence*p*si_a;
    channel2Max = mu_a+confidence*(2-p)*si_a;
    BW1 = (lab_imgae_fg(:,:,2) >= channel2Min ) & (lab_imgae_fg(:,:,2) <= channel2Max);
    maskedRGBImage = I; maskedRGBImage(repmat(bwareaopen(imopen(bwareaopen(~BW1,Threshold),strel('square',10)),round(Threshold/2)),[1 1 3])) = 0;
    % Ripetiamo quello che abbiamo fatto prima sulla banda b
    b_band = lab_image_bg(:,:,3); 
    mu_b = mean2(b_band);
    si_b = std2(b_band);
    ll1 = b_band(:); ll1(ll1<=mu_b) = [];
    ll2 = b_band(:); ll2(ll2>=mu_b) = [];
    % Controlliamo quanti pixel ci sono a sinistra e a destra della media per aggiustare la soglia
    p = size(ll1,1)/size(ll2,1);

    channel2Min = mu_b-confidence*p*si_b;
    channel2Max = mu_b+confidence*(2-p)*si_b;
    BW2 = (lab_imgae_fg(:,:,3) >= channel2Min ) & (lab_imgae_fg(:,:,3) <= channel2Max);
    maskedRGBImage(repmat(bwareaopen(imopen(bwareaopen(~BW2,Threshold),strel('square',10)),round(Threshold/2)),[1 1 3])) = 0;

    result = I-maskedRGBImage;
    bw = im2bw(result,0);
    stats = regionprops(bw,'BoundingBox','Area','Image');
    CC = bwconncomp(bw);
    for i=1:size(stats)
        if(stats(i).BoundingBox(3) > 2*stats(i).BoundingBox(4))
            if(stats(i).Area < 10000)
                bw(CC.PixelIdxList{i}) = 0;
            end
        end
    end
    figure, imshow(maskedRGBImage);
    inImg = (I-maskedRGBImage).*repmat(uint8(bw),[1 1 3]);
    figure, imshow(inImg);
end