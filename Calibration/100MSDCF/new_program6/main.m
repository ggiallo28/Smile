close all; clear all; clc;
figure,imshow(imread('checkerboard.jpg'));
checker_vector = reshape([[0,0,0;255,0,255];[0,0,0;0,255,255];[0,0,0;255,255,0];[255,255,255;255,0,0];[255,255,255;0,255,0];[255,255,255;0,0,255]],[2,6,3]);
checker_center = [0.5*size(checker_vector,2),0.5*size(checker_vector,2)+1];
% path = '../foto/ultime/';
% name = ['DSC00',num2str(545)];
% orig = imread([path,name,'.JPG']);
% orig_bg = imread([path,'DSC00',num2str(547),'.JPG']);
path = '../foto/Evaluation/4r_c/2/';
orig = imread([path,'fore.JPG']);
orig_bg = imread([path,'back.JPG']);
%% Segmentazione
%% Normalizzazione
Container = objContainer();
[Container.I, Container.I_BG, Container.O, Container.O_BG, Container.BB] = normalize_image(orig, orig_bg);
%% SEGMENTATION
segmentation;
%% Identifica colori usando Background, Error Check
Container.obj_chess = error_check(obj_chess(1), obj_chess(2), obj_chess(3),  obj_chess(4));
%% Plot Results
fuse = show_result(Container, false, path);
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
Container = points_allign(Container);
%% Mapping sulla superficie del Manifold
Container = generate_mapping(Container, checker_vector);
show_mapping(Container);
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