color_markers = [ 255 0 0;
                    0 0 0; 
                  0 0 255; 
                255 255 0; 
                  0 255 0; 
              255 255 255; 
                0 255 255; 
                255 0 255]./255;
nClasses = size(color_markers,1);
color_labels = 0:nClasses-1;
distance = zeros([size(R),nClasses]);

% Perform classification

for count = 1:nClasses
    distance(:,:,count) = ( (R-color_markers(count,1)).^2 + (G-color_markers(count,2)).^2 +(B-color_markers(count,3)).^2).^0.5;
end

[value, label] = min(distance,[],3);
label = color_labels(label);

colors = [255 0 0;
            0 0 0; 
          0 0 255; 
        255 255 0; 
          0 255 0; 
      255 255 255; 
        0 255 255; 
        255 0 255]./255;
 
[r, c,p] = size(O);
y = zeros(size(O));
l = double(label)+1;
for m=1:r
	for n=1:c
		y(m,n,:) = colors(l(m,n),:);
	end
end

figure, imshow(y);
colorbar

bw = uint8(value<0.5);
bw = double(cat(3,bw,bw,bw));
imshow(y.*bw);

% scatter plot for the nearest neighbor classification
purple = [119/255 73/255 152/255];
l_blue = [0 1 1];
plot_labels = {'r','k','b','y','g','w',l_blue, purple};


R = im2double(I(:,:,1));
G = im2double(I(:,:,2));
B = im2double(I(:,:,3));
cmyk = rgb2cmyk(I);

yellow = im2bw(cymk(:,:,3),graythresh(cymk(:,:,3)));
yellow = yellow - red;

blu = B-R-G+im2double(cmyk(:,:,1));
blu(blu<0) = 0;
blu = im2bw(blu,graythresh(blu));


red = R-G-B;                %ok
red(red<0) = 0;             %ok
red = im2bw(red+double(cymk(:,:,2))./100,graythresh(red+double(cymk(:,:,2))./100)); %ok
red = imopen(red,strel('square',10));   %ok
yellow = R+G-2*B;           %ok
T = graythresh(yellow);     %ok
yellow = im2bw(yellow,T);   %ok
yellow = yellow-red;        %ok
green = G-0.5*R-0.5*B;                      %ok
T = graythresh(green);                      %ok
green = im2bw(green,T)-yellow;              %ok
green(green<0) = 0;                         %ok
green = imopen(green,strel('square',10));   %ok

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


mask = light_blue | blu | green | yellow | red;
mask = imfill(imdilate(mask,strel('square',10)),'holes');

violet = R.*(~mask)+B.*(~mask)-2*G.*(~mask);
violet(violet<0) = 0;
T = graythresh(violet);                      %ok
violet = im2bw(violet,T);                    %ok

blackwhite = R+G+B;





