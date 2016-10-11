close all; clear all; clc;
% A = [];
% for i = 2: 7
%     orig = imread(['_DSC007',num2str(i),'.JPG']);
%     I = imcrop(orig);
%     R = double(I(:,:,1));
%     G = double(I(:,:,2));
%     B = double(I(:,:,3));
%     lumMat = (0.2126 * R) + (0.7152 * G) + (0.0722 * B); % 
%     K = zeros(size(I));
%     K(:,:,1) = R./(255-lumMat);
%     K(:,:,2) = G./(255-lumMat);
%     K(:,:,3) = B./(255-lumMat);
%     K = im2uint8(K);
%     A = [A;K];
% end

orig = imread(['_DSC007',num2str(7),'.JPG']);
I = imcrop(orig);
% [edge_magnitude, edge_orientation] = coloredges(I);
% img1=rgb2ycbcr(I);
% dx1=edge(img1(:,:,1),'canny');
% med = medfilt2(edge_magnitude,[5,5]);
% res = imclose(dx1,strel('disk',5));
% res = median(median(med)).*im2double(res1);
% med = res + med;
% T = adaptthresh(edge_magnitude,0.3,'ForegroundPolarity','bright');
% BW = imbinarize(edge_magnitude,T);
% BW = medfilt2(BW,[9,9]);
% imshow([BW;med;T]);
R = double(I(:,:,1));
G = double(I(:,:,2));
B = double(I(:,:,3));
lumMat = (0.2126 * R) + (0.7152 * G) + (0.0722 * B); % 
K = zeros(size(I));
K(:,:,1) = R./(255-lumMat);
K(:,:,2) = G./(255-lumMat);
K(:,:,3) = B./(255-lumMat);
% K(:,:,1) = R./(R+G+B);
% K(:,:,2) = G./(R+G+B);
% K(:,:,3) = B./(R+G+B);
I = im2uint8(K);

% nColors = 7;
% kmeans: http://stackoverflow.com/questions/25691735/matlab-color-based-segmentation
%% Segmentazione
[BW_red,inv_BW_red, RedImage, centers_red] = createMask(I,'red');
[BW_green,inv_BW_green, GreenImage, centers_green] = createMask(I,'green');
[BW_blue,inv_BW_blue, BlueImage, centers_blue] = createMask(I,'blue');
[BW_yellow,inv_BW_yellow, YellowImage, centers_yellow] = createMask(I,'yellow');
inv_BW = uint8(inv_BW_red|inv_BW_green|inv_BW_blue|inv_BW_yellow);
col_BW = uint8(BW_red|BW_green|BW_blue|BW_yellow);
BW = 255*(inv_BW-col_BW);
fuse = RedImage+GreenImage+BlueImage+YellowImage+repmat(BW,[1 1 3]);
imshow(fuse); hold on;
for i = 1:size(centers_red,4)
    for k=1:5
        for j=1:5
            scatter(centers_red(k,j,1,i),centers_red(k,j,2,i)) % Riferimento Y verso il basso
        end
    end
end
for i = 1:size(centers_green,4)
    for k=1:5
        for j=1:5
            scatter(centers_green(k,j,1,i),centers_green(k,j,2,i)) % Riferimento Y verso il basso
        end
    end
end
for i = 1:size(centers_yellow,4)
    for k=1:5
        for j=1:5
            scatter(centers_yellow(k,j,1,i),centers_yellow(k,j,2,i)) % Riferimento Y verso il basso
        end
    end
end
for i = 1:size(centers_blue,4)
    for k=1:5
        for j=1:5
            scatter(centers_blue(k,j,1,i),centers_blue(k,j,2,i)) % Riferimento Y verso il basso
        end
    end
end
hold off;
%% Fissa il centro degli assi
BW = uint8(255*(inv_BW_red-BW_red));
Red = RedImage+repmat(BW,[1 1 3]);
BW = uint8(255*(inv_BW_yellow-BW_yellow));
Yellow = YellowImage+repmat(BW,[1 1 3]);
ry = Red+Yellow;

img = (circshift(Yellow,1)-circshift(Yellow,-1))+(circshift(Red,1)-circshift(Red,-1));
 
% imshowpair(Red,Yellow,'Montage');
% col = 15;
% for i=1:size(fuse,2)-col
%     imshow(fuse(:,i:col+i,:));
% end
 
% figure, imshow(fuse);
% [c r] = getpts(1);
% loc = int32([c r]);
% if size(loc)>1
%     loc = [loc(1,1) loc(1,2)];
% end

% I = im2double(I);
% r1 = (I(:,:,1)-I(:,:,2))./(I(:,:,1)+I(:,:,2)); %Red
% r5 = (I(:,:,3)-I(:,:,1))./(I(:,:,3)+I(:,:,1)); %Blue
% r2 = (I(:,:,1)-I(:,:,3))./(I(:,:,1)+I(:,:,3)); %-Red & Yellow
% r7 = r2-r1;                                    %Yellow
% r4 = (I(:,:,2)-I(:,:,1))./(I(:,:,2)+I(:,:,1)); %-Green & Blue
% r8 = r4-r5;                                    %Blue
% r3 = (I(:,:,2)-I(:,:,3))./(I(:,:,2)+I(:,:,3)); %Yellow
% r6 = (I(:,:,3)-I(:,:,2))./(I(:,:,3)+I(:,:,2));
 
% figure, imshow(r7)
 
% figure, imshow(I);
% Red = I(:,:,1);
% Green = I(:,:,2);
% Blue = I(:,:,3);
% corners = detectFASTFeatures(Red);
% imshow(I); hold on;
% plot(corners);
% corners = detectFASTFeatures(Green);
% plot(corners);
% corners = detectFASTFeatures(Blue);
% plot(corners);
 
% s = 30;
% Red_corners = corner(Red,'Harris');
% Green_corners = corner(Green,'Harris');
% Blue_corners = corner(Blue,'Harris');
% figure,imshow(I), hold on, scatter(Red_corners(:,1),Red_corners(:,2));
% scatter(Green_corners(:,1),Green_corners(:,2));
% scatter(Blue_corners(:,1),Blue_corners(:,2));
 
% Gray = rgb2gray(I);
% Gray_corners = corner(Gray,'Harris');
% scatter(Gray_corners(:,1),Gray_corners(:,2));
 
% rgbImage = imread('_DSC0072.JPG');
% % Extract individual color channels and cast to double to avoid clipping.
% R = double(rgbImage(:,:,1));
% G = double(rgbImage(:,:,2));
% B = double(rgbImage(:,:,3));
% blueRatio = uint8(((100 * B)./(1+R+G)) .* (256./(1+B+R+G)));
% % Display images
% fontSize = 28;
% subplot(2, 3, 1);
% imshow(rgbImage);
% title('Original RGB Image', 'FontSize', fontSize);
% subplot(2, 3, 2);
% imshow(uint8(R)); % Must be uint8 for display
% title('R Image', 'FontSize', fontSize);
% subplot(2, 3, 3);
% imshow(uint8(G)); % Must be uint8 for display
% title('G Image', 'FontSize', fontSize);
% subplot(2, 3, 4);
% imshow(uint8(B)); % Must be uint8 for display
% title('B Image', 'FontSize', fontSize);
% subplot(2, 3, 5);
% imshow(blueRatio);
% title('Blue Ratio Image', 'FontSize', fontSize);
% % Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % Give a name to the title bar.
% set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off') 
