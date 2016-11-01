im1 = imread('DSC00130.JPG');
[im1, BB] = imcrop(im1);
R = im2double(im1(:,:,1));
G = im2double(im1(:,:,2));
B = im2double(im1(:,:,3));
L = (0.2126 * R) + (0.7152 * G) + (0.0722 * B); % 
K = zeros(size(im1));
K(:,:,1) = R.*(2-L);
K(:,:,2) = G.*(2-L);
K(:,:,3) = B.*(2-L);
im1 = im2uint8(K);
im1 = imfilter(im1,[1 1 1]./3);

im2 = imcrop(imread('DSC00129.JPG'),BB);
R = im2double(im2(:,:,1));
G = im2double(im2(:,:,2));
B = im2double(im2(:,:,3));
L = (0.2126 * R) + (0.7152 * G) + (0.0722 * B); % 
K = zeros(size(im2));
K(:,:,1) = R.*(2-L);
K(:,:,2) = G.*(2-L);
K(:,:,3) = B.*(2-L);
im2 = im2uint8(K);
im2 = imfilter(im2,[1 1 1]./3);

lab2 = rgb2lab(im2);
lab1 = rgb2lab(im1);
mu = mean2(lab2(:,:,2));
si = std2(lab2(:,:,2));
ll = lab2(:,:,2); 
ll1 = ll(:); ll1(ll1<=mu) = [];
ll2 = ll(:); ll2(ll2>=mu) = [];
p = size(ll1,1)/size(ll2,1);
[r,c] = size(im1);
T = round(r*c/7000);
channel2Min = mu-2.58*p*si; %2.7 forse meglio
channel2Max = mu+2.58*(2-p)*si;
BW1 = (lab1(:,:,2) >= channel2Min ) & (lab1(:,:,2) <= channel2Max);
maskedRGBImage_im1 = im1;
maskedRGBImage_im1(repmat(bwareaopen(imopen(bwareaopen(~BW1,T),strel('square',10)),round(T/2)),[1 1 3])) = 0;

mu = mean2(lab2(:,:,3));
si = std2(lab2(:,:,3));
ll = lab2(:,:,3); 
ll1 = ll(:); ll1(ll1<=mu) = [];
ll2 = ll(:); ll2(ll2>=mu) = [];
p = size(ll1,1)/size(ll2,1);

channel2Min = mu-2.58*p*si;
channel2Max = mu+2.58*(2-p)*si;
BW2 = (lab1(:,:,3) >= channel2Min ) & (lab1(:,:,3) <= channel2Max);
maskedRGBImage_im1(repmat(bwareaopen(imopen(bwareaopen(~BW2,T),strel('square',10)),round(T/2)),[1 1 3])) = 0;

res = im1-maskedRGBImage_im1;
bw = im2bw(res,0);
stats = regionprops(bw,'BoundingBox','Area','Image');
CC = bwconncomp(bw);
for i=1:size(stats)
    if(stats(i).BoundingBox(3) > 2*stats(i).BoundingBox(4))
        if(stats(i).Area < 10000)
            bw(CC.PixelIdxList{i}) = 0;
        end
    end
end
figure, imshow(maskedRGBImage_im1);
figure, imshow((im1-maskedRGBImage_im1).*repmat(uint8(bw),[1 1 3]));
inImg = (im1-maskedRGBImage_im1).*repmat(uint8(bw),[1 1 3]);
R = inImg(:,:,1); R = R(bw);
G = inImg(:,:,2); G = G(bw);
B = inImg(:,:,3); B = B(bw);
RGB = cat(3,R,G,B);
hsv = rgb2hsv(RGB);
mpd = 15;
[output_peak, output_minima_low, output_minima_high, output_minima_mid, hist_size] = findlocalminima(hsv(:,:,1),mpd,5,0,1);
% Non sono in grado d distinguere tra il viola e il rosso, tra l'azzurro e il blu, quindi 4 cluster invece che 6
figure, imshow(imread('hsv.jpg'));
% RED
inImg = im2double(inImg);
hsv_inImg = rgb2hsv(inImg);
red_th_min = output_minima_mid(1)/hist_size;
red_th_max = output_minima_mid(4)/hist_size;
bw_red = (hsv_inImg(:,:,1)<=red_th_min | hsv_inImg(:,:,1)>=red_th_max) & bw;
bw_red = bwareaopen(bw_red,T);
bw_red = imclose(bw_red,strel('square',10));
figure, imshow(inImg.*repmat(bw_red,[1 1 3]));
% YELLOW
yellow_th_min = output_minima_mid(1)/hist_size;
yellow_th_max = output_minima_mid(2)/hist_size;
bw_yellow = hsv_inImg(:,:,1)>=yellow_th_min & hsv_inImg(:,:,1)<=yellow_th_max & bw;
bw_yellow = bwareaopen(bw_yellow,T);
bw_yellow = imclose(bw_yellow,strel('square',10));
figure, imshow(inImg.*repmat(bw_yellow,[1 1 3]));
% GREEN
green_th_min = output_minima_mid(2)/hist_size;
green_th_max = output_minima_mid(3)/hist_size;
bw_green = hsv_inImg(:,:,1)>=green_th_min & hsv_inImg(:,:,1)<=green_th_max & bw;
bw_green = bwareaopen(bw_green,T);
bw_green = imclose(bw_green,strel('square',10));
figure, imshow(inImg.*repmat(bw_green,[1 1 3]));
% BLU
blu_th_min = output_minima_mid(3)/hist_size;
blu_th_max = output_minima_mid(4)/hist_size;
bw_blu = hsv_inImg(:,:,1)>=blu_th_min & hsv_inImg(:,:,1)<=blu_th_max & bw;
bw_blu = bwareaopen(bw_blu,T);
bw_blu = imclose(bw_blu,strel('square',10));
figure, imshow(inImg.*repmat(bw_blu,[1 1 3]));








% noOfColors = 4;
% hsv = rgb2hsv(inImg);
% s_img = size(inImg);
% inputImg = double(hsv(:,:,1))+1;
% disp('K-Means Processing Started');
% opts = statset('Display','final');
% [idx, C] = kmeans(inputImg, noOfColors, 'Distance','cosine','Replicates',5, 'Options',opts);
%                           
% disp('K-Means Processing Completed');
% palette = round(C);
% 
% %Color Mapping
% idx = uint8(idx);
% outImg = zeros(s_img(1),s_img(2),2);
% temp = reshape(idx, [s_img(1) s_img(2)]);
% for i = 1 : 1 : s_img(1)
%     for j = 1 : 1 : s_img(2)
%         outImg(i,j,:) = palette(temp(i,j),:);
%     end
% end
% lumMat = lab(:,:,1).*bw;
% lumMat = bw.*(sum(sum(lumMat))./sum(sum(bw)));
% outImg = cat(3,lumMat,outImg);
% imshow(lab2rgb(outImg));
% 
% 
%  
% diffImage = imabsdiff(im2double(im2), im2double(im1));
% foregroundMask = zeros(size(diffImage));
% 
% threshold = 200.0/255;
% for j=1:size(diffImage,1)
%     for i=1:size(diffImage,2)
%         pix = diffImage(j,i,:);
%         dist = pix(1)*pix(1) + pix(2)*pix(2) + pix(3)*pix(3);
%         dist = sqrt(dist);
%         if dist>threshold
%             foregroundMask(j,i) = 255;
%         end
%     end
% end
% imshow(foregroundMask);
%         
% diff = lab1(:,:,2)-lab2(:,:,2);
% imshow(diff)
% 
% 
% diff = sqrt((im1(:,:,1)-im2(:,:,1)).^2+(im1(:,:,3)-im2(:,:,3)).^2+(im1(:,:,2)-im2(:,:,2)).^2);
% 
% im3 = im1;
% im3(:,:,1) = abs(im2(:,:,1)-im1(:,:,1));
% im3(:,:,2) = abs(im2(:,:,2)-im1(:,:,2));
% im3(:,:,3) = abs(im2(:,:,3)-im1(:,:,3));


% 
% im1 = im2double(imread('DSC00127.JPG'));
% im2 = im2double(imread('DSC00128.JPG'));
% im4 = im2double(imread('DSC00130.JPG'));
% [im4, BB] = imcrop(im4);
% im3 = imcrop(im2double(imread('DSC00129.JPG')),BB);
% im5 = imcrop(im2double(imread('DSC00131.JPG')),BB);
% im6 = imcrop(im2double(imread('DSC00132.JPG')),BB);
% im7 = imcrop(im2double(imread('DSC00133.JPG')),BB);
% im8 = imcrop(im2double(imread('DSC00134.JPG')),BB);
% im9 = imcrop(im2double(imread('DSC00135.JPG')),BB);
% im10 = imcrop(im2double(imread('DSC00136.JPG')),BB);
% im11 = imcrop(im2double(imread('DSC00137.JPG')),BB);
% im12 = imcrop(im2double(imread('DSC00138.JPG')),BB);
% im13 = imcrop(im2double(imread('DSC00139.JPG')),BB);
% im14 = imcrop(im2double(imread('DSC00140.JPG')),BB);
% im15 = imcrop(im2double(imread('DSC00141.JPG')),BB);
% im16 = imcrop(im2double(imread('DSC00142.JPG')),BB);
% 
% b(:,:,:,1) = im4;
% b(:,:,:,2) = im5;
% b(:,:,:,3) = im6;
% b(:,:,:,4) = im7;
% b(:,:,:,5) = im8;
% 
% imW = median(b,4);
% 
% b(:,:,:,1) = im3;
% b(:,:,:,2) = im9;
% b(:,:,:,3) = im10;
% b(:,:,:,4) = im11;
% b(:,:,:,5) = im12;
% b(:,:,:,6) = im13;
% b(:,:,:,7) = im14;
% b(:,:,:,8) = im15;
% b(:,:,:,9) = im16;
% 
% imO = median(b,4);
% 
% imshow([imO;imW]);
% 
% imr = imhist(
% imshow(imW-imO);


