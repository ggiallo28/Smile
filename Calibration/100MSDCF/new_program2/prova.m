clear all; close all;
orig = imread(['foto/DSC0010',num2str(9),'.JPG']);
%% Normalizzazione
figure, O = imcrop(orig);
R = im2double(O(:,:,1));
G = im2double(O(:,:,2));
B = im2double(O(:,:,3));
L = (0.2126 * R) + (0.7152 * G) + (0.0722 * B); % 
K = zeros(size(O));
K(:,:,1) = imadjust(R.*(2-L));
K(:,:,2) = imadjust(G.*(2-L));
K(:,:,3) = imadjust(B.*(2-L));
I = im2uint8(K);

R = im2double(I(:,:,1));
G = im2double(I(:,:,2));
B = im2double(I(:,:,3));
cmyk = rgb2cmyk(I);

figure, imshow(I);
red = R-G-B;                %ok
red(red<0) = 0;             %ok
red = im2bw(red+double(cmyk(:,:,2))./100,graythresh(red+double(cmyk(:,:,2))./100)); %ok
red = imopen(red,strel('square',10));   %ok
figure, imshowpair(I,red,'falsecolor');
title('red');
yellow = R+G-2*B;           %ok
T = graythresh(yellow);     %ok
yellow = im2bw(yellow,T);   %ok
yellow = yellow-red;        %ok
figure, imshowpair(I,yellow,'falsecolor');
title('yellow');
green = G-0.5*R-0.5*B;                      %ok
T = graythresh(green);                      %ok
green = im2bw(green,T)-yellow;              %ok
green(green<0) = 0;                         %ok
green = imopen(green,strel('square',5));   %ok
figure, imshowpair(I,green,'falsecolor');
title('green');

blu = B-R-G+double(cmyk(:,:,1))./100;               %ok
T = graythresh(blu);                                %ok
blu = im2bw(blu,T)-green;                           %ok 
light_blue = G+B-2*R;                               %ok
T = graythresh(light_blue);                         %ok
light_blue = im2bw(light_blue,T) & blu;             %ok
blu = ~light_blue & blu;                            %ok
light_blue = light_blue - blu;                      %ok
blu = imopen(blu,strel('square',5));                %ok  
light_blue = imopen(light_blue,strel('square',5));  %ok
figure, imshowpair(I,light_blue,'falsecolor');
title('light_blue');
figure, imshowpair(I,blu,'falsecolor');
title('blu');

mask = light_blue | blu | green | yellow | red;
mask = imfill(imdilate(mask,strel('square',10)),'holes');

violet = R.*(~mask)+B.*(~mask)-2*G.*(~mask);
violet(violet<0) = 0;
T = graythresh(violet);                      %ok
violet = im2bw(violet,T);                    %ok
figure, imshowpair(I,violet,'falsecolor');
title('violet');

blackwhite = R+G+B;

mul = 2;
raw_green = G>B & G>R;
stats = regionprops(raw_green,'Area');
T = mean([stats.Area]);
raw_green = bwareaopen(raw_green, mul*round(T));

raw_red = R>B & R>G;
stats = regionprops(raw_red,'Area');
T = mean([stats.Area]);
raw_red = bwareaopen(raw_red, mul*round(T));

raw_blu = B>R & B>G;
stats = regionprops(raw_blu,'Area');
T = mean([stats.Area]);
raw_blu = bwareaopen(raw_blu, mul*round(T));

raw_yellow = R>B & G>B;
stats = regionprops(raw_yellow,'Area');
T = mean([stats.Area]);
raw_yellow = bwareaopen(raw_yellow, mul*round(T));

raw_light_blu = B>R & G>R;
stats = regionprops(raw_light_blu,'Area');
T = mean([stats.Area]);
raw_light_blu = bwareaopen(raw_light_blu, mul*round(T));

raw_violet = R>G & R>B;
stats = regionprops(raw_violet,'Area');
T = mean([stats.Area]);
raw_violet = bwareaopen(raw_violet, mul*round(T));

raw_green = raw_green & ~raw_red & ~raw_blu;
raw_red = ~raw_green & raw_red & ~raw_blu & ~raw_light_blu & raw_violet;
raw_blu = ~raw_green & ~raw_red & raw_blu & ~raw_yellow & ~raw_light_blu & ~raw_violet;
raw_yellow = ~raw_green & ~raw_red & ~raw_blu & raw_yellow & ~raw_light_blu & ~raw_violet;
raw_light_blu = ~raw_green & ~raw_red & ~raw_blu & ~raw_yellow & raw_light_blu & ~raw_violet;
raw_violet = ~raw_green & ~raw_red & ~raw_blu & ~raw_yellow & ~raw_light_blu & raw_violet;




