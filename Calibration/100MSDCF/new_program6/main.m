close all; clear all; clc;
figure,imshow(imread('checkerboard.jpg'));
checker_vector = reshape([[0,0,0;255,0,255];[0,0,0;0,255,255];[0,0,0;255,255,0];[255,255,255;255,0,0];[255,255,255;0,255,0];[255,255,255;0,0,255]],[2,6,3]);
checker_center = [0.5*size(checker_vector,2),0.5*size(checker_vector,2)+1];
path = '../foto/ultime/';
name = ['DSC00',num2str(295)];
orig = imread([path,name,'.JPG']);
orig_bg = imread([path,'DSC00',num2str(296),'.JPG']);
%% Normalizzazione
Container = objContainer();
[Container.I, Container.I_BG, Container.O, Container.O_BG, Container.BB] = normalize_image(orig, orig_bg);
%% Parametri
Container.num_square = 5;
Container.confidence = 3;
Container.img_dim  = size(Container.I);
%Container.Threshold = round(Container.img_dim(1)*Container.img_dim(2)/7000); % Soglia dimensione blob normalizzata alla dimensione dell'immagine
Container.fraction = 20;
Container.mpd = 30; % 30
Container.windowSize = 6; % 6
Container.op_th = 15;
if exist([path,name,'.mat'], 'file') == 2
    load([path,name,'.mat'])
end
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
Container = th_estimation(Container, bw);
% bisogna farlo per ogni riflesso siccome il massimo quadrato ha dimensione
% diversa
%% Segmentazione
obj_red = computeChess(Container, bw, hist_size, output_minima_mid, 'red');
obj_green = computeChess(Container, bw, hist_size, output_minima_mid, 'green');
obj_blue = computeChess(Container, bw, hist_size, output_minima_mid, 'blue');
obj_yellow = computeChess(Container, bw, hist_size, output_minima_mid, 'yellow');
[obj_chess, transtions] = getColorBoundary(obj_red, obj_green, obj_blue, obj_yellow, checker_vector, Container);
obj_chess = completeComputeChess(obj_chess, transtions, Container);
%% Identifica colori usando Background, Error Check
Container.obj_chess = error_check(obj_chess(1), obj_chess(2), obj_chess(3),  obj_chess(4));
%% Plot Results
fuse = show_result(Container);
%% Separa i riflessi
Container = split_refletions(Container, fuse, checker_vector);
%% Creazione Maschere per separazione riflessi
Container = generate_mask(Container);%obj_chess, label, order, positions, types, size(fuse));
%% Calcolo convexhull dei riflessi: controllare se è necessario fare sta cosa
Container = generate_bwconvhull(Container);
%% Calcolo Convexhull tessere singole
for i=1:size(Container.obj_chess,1)
    if ( ~Container.obj_chess(i).isEmpty )
        for j=1:size(Container.obj_chess(i).chess,2)
            cut_x = Container.obj_chess(i).bbox_x(j,:);
            cut_y = Container.obj_chess(i).bbox_y(j,:); 
            Container.obj_chess(i).chess(j).ch_mask = false(size(fuse,1),size(fuse,2));
            mask = imdilate(bwconvhull(Container.obj_chess(i).chess(j).mask),strel('square',3));
            Container.obj_chess(i).chess(j).ch_mask(cut_y(1):cut_y(2),cut_x(1):cut_x(2)) =...
                mask(cut_y(1):cut_y(2),cut_x(1):cut_x(2));     
            figure, imshowpair(Container.obj_chess(i).chess(j).ch_mask,rgb2gray(Container.I),'falsecolor');
        end
    end
end
close all;
%% Full Mask
FullMask  = false(size(fuse,1),size(fuse,2));
if ( Container.isL2 )
    FullMask = FullMask | Container.maskL2I;
end
if ( Container.isL1 )
    FullMask = FullMask | Container.maskL1I;
end
if ( Container.isC )
    FullMask = FullMask | Container.maskCI;
end
if ( Container.isR1 )
    FullMask = FullMask | Container.maskR1I;
end
if ( Container.isR2 )
    FullMask = FullMask | Container.maskR2I;
end
Container.FullMask = FullMask;
%% Fissa il centro degli assi
[left_center_axis, right_center_axis, mid_center_axis] = generate_central_axis(Container);
%% Corner
pointsArray = calculate_corners(Container, left_center_axis, right_center_axis, mid_center_axis);
%% Divido i punti per ogni componente
Container = points_split(Container, pointsArray);
%% Padding delle matrici e allineamento
Container3 = points_allign(Container);
%% Mapping sulla superficie del Manifold
Container3 = generate_mapping(Container3, checker_vector);
show_mapping(Container3);
%% Calcolo angolo degli specchi
% Pleft = [];
% Pright = [];
% for l=1:size(obj_chess,1)
%     if ( ~obj_chess(l).isEmpty )
%         for i=1:size(obj_chess(l).chess(1).intersections_x,1)
%             for j=1:size(obj_chess(l).chess(1).intersections_x,2)
%                 vect_x = zeros(2,3);
%                 vect_y = zeros(2,3);
%                 idk = 1;
%                 for k =1:size(obj_chess(l).chess,2)           
%                     if (~strcmp(obj_chess(l).chess(k).type,types(2)))
%                         vect_x(1,idk) = obj_chess(l).chess(k).intersections_x(i,j); vect_x(2,idk) = k;
%                         vect_y(1,idk) = obj_chess(l).chess(k).intersections_y(i,j); vect_y(2,idk) = k;
%                         idk = idk +1;
%                     end                                        
%                 end
%                 vect_xleft = 0; vect_yleft = 0; vect_xmid = 0; vect_ymid = 0; vect_xright = 0; vect_yright = 0;
%                 if idk > 2
%                     for k=1:size(vect_x,2)
%                         idk = find(order(2,:)==vect_x(2,k) & order(1,:) == l);
%                         if ~isempty(idk)
%                             pos = label(3,idk);
%                             type = label(4,idk);
%                             assert(~strcmp(type,types(2)));
%                             if (strcmp(pos,positions(1)))
%                                 vect_xleft = vect_x(1,k);
%                                 vect_yleft = vect_y(1,k);
%                             end
%                             if (strcmp(pos,positions(2)))
%                                 vect_xmid = vect_x(1,k);
%                                 vect_ymid = vect_y(1,k);
%                             end
%                             if (strcmp(pos,positions(3)))
%                                 vect_xright = vect_x(1,k);
%                                 vect_yright = vect_y(1,k);
%                             end
%                         end
%                     end
%                     if(~(vect_xleft == 0 || vect_yleft == 0 || vect_xmid == 0 || vect_ymid == 0))
%                         Pleft = [Pleft;[vect_xmid,vect_ymid,vect_xleft,vect_yleft]];
%                     end
%                     if(~(vect_xright == 0 || vect_yright == 0 || vect_xmid == 0 || vect_ymid == 0))
%                         Pright = [Pright;[vect_xmid,vect_ymid,vect_xright,vect_yright]];
%                     end
%                 end
%             end
%         end
%     end
% end
% Pleft = unique(Pleft,'rows');
% Pright = unique(Pright,'rows');
% imshow(I); hold on;
% scatter([Pleft(:,1);Pleft(:,3)],[Pleft(:,2);Pleft(:,4)]);
% pause
% imshow(I); hold on;
% scatter([Pright(:,1);Pright(:,3)],[Pright(:,2);Pright(:,4)]);
%% Ora come associo i punti se levo il cilindro? Con il codice che già hanno ma usando la griglia.
%% TODO: Il prof ha detto di stimare la differenza di colore negli histogrammi del rosso al centro e del rosso a destra

%% Come stimo l'angolo degli specchi?

%% Salvataggio Setup
% save([path,name,'.mat'],'confidence','Threshold','mpd','windowSize','op_th');