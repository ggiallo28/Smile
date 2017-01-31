close all; clear all; clc;
addpath('E:\ProgrammiHD\mexopencv-master');
figure,imshow(imread('checkerboard.jpg'));
checker_vector = reshape([[0,0,0;255,0,255];[0,0,0;0,255,255];[0,0,0;255,255,0];[255,255,255;255,0,0];[255,255,255;0,255,0];[255,255,255;0,0,255]],[2,6,3]);
checker_center = [0.5*size(checker_vector,2),0.5*size(checker_vector,2)+1];
% path = 'E:\Tesi\Ottica\Calibration\100MSDCF\foto\Barbara\5\';
% name = ['DSC00',num2str(563)];
% orig = imread([path,name,'.JPG']);
% orig_bg = imread([path,'DSC00',num2str(564),'.JPG']);
path = '../foto/Evaluation/4l_c/10/';
orig = imread([path,'fore.JPG']);
orig_bg = imread([path,'back.JPG']);
%% Normalizzazione
Container = objContainer(false);
if exist([path,'bb.mat'], 'file') == 2
    bb = load([path,'bb.mat']);
    bb = bb.bb;
else
    bb = [];
end
if exist([path,'toflip'], 'file')
    checker_vector = fliplr(checker_vector);
end
[Container.I, Container.I_BG, Container.O, Container.O_BG, Container.BB] = normalize_image(orig, orig_bg, bb);
%% Parametri
Container.num_square = 5;
if exist([path,'confidence'], 'file') == 2
    Container.confidence = 3.2;
else
    Container.confidence = 2.8;
end
Container.img_dim  = size(Container.I);
%Container.Threshold = round(Container.img_dim(1)*Container.img_dim(2)/7000); % Soglia dimensione blob normalizzata alla dimensione dell'immagine
Container.fraction = 20;
Container.mpd = 30; % 30, 25
Container.windowSize = 6; % 6
Container.op_th = 15;
%% SEGMENTATION
segmentation; % O(numero pixel)
obj_chess = completeComputeChess(obj_chess, transtions, Container);
%% Identifica colori usando Background, Error Check
Container.obj_chess = error_check(obj_chess(1), obj_chess(2), obj_chess(3),  obj_chess(4));
%% Plot Results
fuse = show_result(Container, false, path);
%% MODELLAZIONE
modeling;
% cyl = get_cylinder(Container);
%% Fissa il centro degli assi
[left_center_axis, right_center_axis, mid_center_axis] = generate_central_axis(Container);
% Corner detector numero_viste*[O(harris) |+ O(numero corner)]
%% COORNER IDENTIFICATION
%% Corner
pointsArray = calculate_corners(Container, left_center_axis, right_center_axis, mid_center_axis, path);
%% Divido i punti per ogni componente
Container = points_split(Container, pointsArray);
% Association numero_viste*O(numero corner)
%% Padding delle matrici e allineamento
Container = points_allign(Container);
%% Mapping sulla superficie del Manifold
Container = generate_mapping(Container, checker_vector);
show_mapping(Container);
points_match = save_mapping(Container);
figure, imshow(orig); hold on;
points_match_offset = points_match;
for i=1:size(points_match,1)
    row = points_match_offset(i,:);
    for k=1:size(row,2)
        EL = cell2mat(row(k));
        if ~isempty(EL)
            if EL(1) == 0 || EL(2) == 0
                row(k) = cell(1,1);
                disp('del');
            else
                row(k) = mat2cell([EL(1)+bb(1),EL(2)+bb(2)],1,2);
                scatter(EL(1)+bb(1),EL(2)+bb(2));
            end
        end
    end
    points_match_offset(i,:) = row;
end
%% Complessità totale
% O(numero pixel) + O(numero settori angolari) + numero_viste*[O(harris) |+ O(numero corner)] + numero_viste*O(numero corner) = 
% O(numero pixel) + O(numero settori angolari) + numero_viste*[O(harris) + O(numero corner)]

% %% MIRROR ANGLE
% idxCentral = 3;
% idxReflected = 2;
% M = []; fx = 6455.61820399344; P=[];
% figure, imshow(orig); hold on;
% for i=1:size(points_match_offset,1)
%     if ~isempty(points_match_offset{i,idxCentral}) && ~isempty(points_match_offset{i,idxReflected})
%        pi = points_match_offset{i,idxCentral};
%        if isempty(P) || isempty(find(P(:,1) == pi(1) & P(:,2) == pi(2), 1))
%             P = [P;pi];
%             scatter(pi(1),pi(2));
%             pip = points_match_offset{i,idxReflected};
%             scatter(pip(1),pip(2));
%             row = [(pi(2) - pip(2))*fx (-pi(1) + pip(1))*fx pi(1)*pip(2)-pi(2)*pip(1)]; 
%             M = [M; row];
%         end
%     end
% end
% LS = M'*M;
% [V,D] = eig(LS);
% D(D==0) = [];
% [err, idx] = min(D);
% u = V(:,idx);
% norm(M*u)
% % 
% % [U,S,V] = svd(M);
% % S(S==0)=[];
% % [err, idx] = min(S);
% % u2 = V(:,idx);
% % norm(M*u2)
% % 
% % [U,S,V] = svd(M,'econ');
% % u = V(:,size(M,2));
% % norm(M*u)
% 
% a = u(1); b = u(2); c = u(3);
% d = (fx*2.6*5)/785; P = []; % http://www.pyimagesearch.com/2015/01/19/find-distance-camera-objectmarker-using-python-opencv/
% for i=1:size(points_match_offset,1)
%     if ~isempty(points_match_offset{i,idxCentral}) && ~isempty(points_match_offset{i,idxReflected})
%         pi = points_match_offset{i,idxCentral};
%         pip = points_match_offset{i,idxReflected};
%         if isempty(P) || isempty(find(P(:,1) == pi(1) & P(:,2) == pi(2), 1))
%             P = [P;pi];
%             g11 = ((2*a*a-1)/(2*fx));
%             g12 = a*b/fx;
%             g13 = a*c; 
%             g14 = 0.5/fx;
%             g21 = g12;
%             g22 = ((2*b*b-1)/(2*fx));
%             g23 = b*c;
%             g24 = g14;
%             g31 = a*c/fx; 
%             g32 = b*c/fx; 
%             g33 = 0.5*(2*c*c-1);
%             g34 = 0.5;
%             
%             G = [g11*pi(1) + g12*pi(2) + g13, pip(1)*g14;
%                  g21*pi(1) + g22*pi(2) + g23, pip(2)*g24
%                  g31*pi(1) + g32*pi(2) + g33, g34];
%              
%             z = pinv(G)*d*u
%             scatter3(pi(1),pi(2),z(1)); hold on; scatter3(pip(1),pip(2),z(2));
%         end
%     end
% end
% %% Calcolo angolo degli specchi
% % Pleft = [];
% % Pright = [];
% % for l=1:size(Container.obj_chess,1)
% %     if ( ~Container.obj_chess(l).isEmpty )
% %         for i=1:size(Container.obj_chess(l).chess(1).intersections_x,1)
% %             for j=1:size(Container.obj_chess(l).chess(1).intersections_x,2)
% %                 vect_x = zeros(2,3);
% %                 vect_y = zeros(2,3);
% %                 idk = 1;
% %                 for k =1:size(Container.obj_chess(l).chess,2)           
% %                     if (~strcmp(Container.obj_chess(l).chess(k).type,Container.types(2)))
% %                         vect_x(1,idk) = Container.obj_chess(l).chess(k).intersections_x(i,j); vect_x(2,idk) = k;
% %                         vect_y(1,idk) = Container.obj_chess(l).chess(k).intersections_y(i,j); vect_y(2,idk) = k;
% %                         idk = idk +1;
% %                     end                                        
% %                 end
% %                 vect_xleft = 0; vect_yleft = 0; vect_xmid = 0; vect_ymid = 0; vect_xright = 0; vect_yright = 0;
% %                 if idk > 2
% %                     for k=1:size(vect_x,2)
% %                         idk = find(Container.order(2,:)==vect_x(2,k) & Container.order(1,:) == l);
% %                         if ~isempty(idk)
% %                             pos = Container.label(3,idk);
% %                             type = Container.label(4,idk);
% %                             assert(~strcmp(type,Container.types(2)));
% %                             if (strcmp(pos,Container.positions(1)))
% %                                 vect_xleft = vect_x(1,k);
% %                                 vect_yleft = vect_y(1,k);
% %                             end
% %                             if (strcmp(pos,Container.positions(2)))
% %                                 vect_xmid = vect_x(1,k);
% %                                 vect_ymid = vect_y(1,k);
% %                             end
% %                             if (strcmp(pos,Container.positions(3)))
% %                                 vect_xright = vect_x(1,k);
% %                                 vect_yright = vect_y(1,k);
% %                             end
% %                         end
% %                     end
% %                     if(~(vect_xleft == 0 || vect_yleft == 0 || vect_xmid == 0 || vect_ymid == 0))
% %                         Pleft = [Pleft;[vect_xmid,vect_ymid,vect_xleft,vect_yleft]];
% %                     end
% %                     if(~(vect_xright == 0 || vect_yright == 0 || vect_xmid == 0 || vect_ymid == 0))
% %                         Pright = [Pright;[vect_xmid,vect_ymid,vect_xright,vect_yright]];
% %                     end
% %                 end
% %             end
% %         end
% %     end
% % end
% % Pleft = unique(Pleft,'rows');
% % Pright = unique(Pright,'rows');
% % imshow(Container.I); hold on;
% % scatter([Pleft(:,1);Pleft(:,3)],[Pleft(:,2);Pleft(:,4)]);
% % pause
% % imshow(Container.I); hold on;
% % scatter([Pright(:,1);Pright(:,3)],[Pright(:,2);Pright(:,4)]);
% 
% %% Ora come associo i punti se levo il cilindro? Con il codice che già hanno ma usando la griglia.
% %% TODO: Il prof ha detto di stimare la differenza di colore negli histogrammi del rosso al centro e del rosso a destra
% 
% %% Come stimo l'angolo degli specchi?
% 
% %% Salvataggio Setup
% % save([path,name,'.mat'],'confidence','Threshold','mpd','windowSize','op_th');
% 
% % P = unique(P,'rows');
% % squareSize = 26; % in millimeters
% % boardSize = [5, 9];
% % worldPoints = generateCheckerboardPoints(boardSize, squareSize);
% % imshow(orig); hold on;
% % for i=1:size(P,1)/4
% %     range = (i-1)*4+1:(i-1)*4+4;
% %     [P(range,2), idx] = sort(P(range,2),'descend');
% %     idx = idx + (i-1)*4;
% %     P(idx,1) = P(idx,1);
% % end
% % % range = 1:4:size(P,1);
% % % P = [P(range,:); P(range+1,:); P(range+2,:); P(range+3,:)];
% % % for i=1:size(P,1)
% % %     scatter(P(i,1),P(i,2)); pause;
% % % end
% % % for i=1:size(worldPoints,1)
% % %     scatter(worldPoints(i,1),worldPoints(i,2)); hold on; pause
% % % end
% % [rotationMatrix, translationVector] = extrinsics(P, worldPoints, cameraParams);
