%% Parametri
Container.num_square = 5;
Container.confidence = 2.8;
Container.img_dim  = size(Container.I);
%Container.Threshold = round(Container.img_dim(1)*Container.img_dim(2)/7000); % Soglia dimensione blob normalizzata alla dimensione dell'immagine
Container.fraction = 20;
Container.mpd = 30; % 30, 25
Container.windowSize = 6; % 6
Container.op_th = 15;
% if exist([path,name,'.mat'], 'file') == 2
%     load([path,name,'.mat'])
% end
%% Background Subtraction
[bw, inImg] = subtract_background(Container);
%% Histogram Peak Finding
Rp = inImg(:,:,1); Rp = Rp(bw);
Gp = inImg(:,:,2); Gp = Gp(bw);
Bp = inImg(:,:,3); Bp = Bp(bw);
RGB = cat(3,Rp,Gp,Bp);
hsv = rgb2hsv(RGB);
[output_peak, output_minima_mid, hist_size ] = findpeaksandminima(hsv(:,:,1),Container.windowSize,Container.mpd);
% [output_peak, output_minima_low, output_minima_high, output_minima_mid, hist_size] =...
%     findlocalminima(hsv(:,:,1),Container.mpd,Container.windowSize,0,1);
% Non sono in grado d distinguere tra il viola e il rosso, tra l'azzurro e il blu, quindi 4 cluster invece che 6
figure, imshow(imread('hsv.jpg'));
%% Stima la dimensione del quadrato più grande
bw_en = imdilate(bw,strel('square',20));
bw_en = bwareaopen(bw_en,floor(0.2*mean(cell2mat(struct2cell(regionprops(bw_en,'Area'))))));
bw = bw & bw_en;
Container = th_estimation(Container, bw);
% bisogna farlo per ogni riflesso siccome il massimo quadrato ha dimensione
% diversa
%% Segmentazione
obj_red = computeChess(Container, bw, hist_size, output_minima_mid, 'red');
obj_green = computeChess(Container, bw, hist_size, output_minima_mid, 'green');
obj_blue = computeChess(Container, bw, hist_size, output_minima_mid, 'blue');
obj_yellow = computeChess(Container, bw, hist_size, output_minima_mid, 'yellow');
[obj_chess, transtions, Container] = getColorBoundary(obj_red, obj_green, obj_blue, obj_yellow, checker_vector, Container);
obj_chess = completeComputeChess(obj_chess, transtions, Container);