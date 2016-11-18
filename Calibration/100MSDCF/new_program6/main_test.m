close all; clear all; clc;
path = '../foto/Evaluation/';
checker_vector = reshape([[0,0,0;255,0,255];[0,0,0;0,255,255];[0,0,0;255,255,0];[255,255,255;255,0,0];[255,255,255;0,255,0];[255,255,255;0,0,255]],[2,6,3]);
checker_center = [0.5*size(checker_vector,2),0.5*size(checker_vector,2)+1];
folders = dir(path);
set(0,'DefaultFigureVisible','on');  % all subsequent figures "off"
for ff=1:size(folders,1)
    if ~(strcmp(folders(ff).name,'..') || strcmp(folders(ff).name,'.'))
        path_2 = [path,folders(ff).name,'/'];
        folders_2 = dir(path_2);
        for ff2=1:size(folders_2,1)
            if ~(strcmp(folders_2(ff2).name,'..') || strcmp(folders_2(ff2).name,'.'))
                path_3 = [path_2,folders_2(ff2).name,'/'];
                if exist([path_3,'skip'], 'file') == 2
                    continue;
                end
                 if exist([path_3,'Container.mat'], 'file') == 2
                    continue;
                end
                disp([path_3,'....']);
                orig = imread([path_3,'fore.JPG']);
                orig_bg = imread([path_3,'back.JPG']);
                %% Normalizzazione
                Container = objContainer();
                if exist([path_3,'bb.mat'], 'file') == 2
                     load([path_3,'bb.mat'])
                     Container.BB = bb;
                     Container.O = imcrop(orig,bb);
                else
                    figure, [Container.O, Container.BB] = imcrop(orig);
                    Container.BB = round(Container.BB);
                    Container.O = imcrop(orig,Container.BB);
                    bb =  Container.BB;
                    save([path_3,'bb.mat'],'bb');
                end
                R = im2double(Container.O(:,:,1));
                G = im2double(Container.O(:,:,2));
                B = im2double(Container.O(:,:,3));
                L = (0.2126 * R) + (0.7152 * G) + (0.0722 * B); 
                K = zeros(size(Container.O));
                K(:,:,1) = R.*(2-L);
                K(:,:,2) = G.*(2-L);
                K(:,:,3) = B.*(2-L);
                Container.I = im2uint8(K);
                Container.O_BG = imcrop(orig_bg,Container.BB);
                Rt = im2double(Container.O_BG(:,:,1));
                Gt = im2double(Container.O_BG(:,:,2));
                Bt = im2double(Container.O_BG(:,:,3));
                L = (0.2126 * Rt) + (0.7152 * Gt) + (0.0722 * Bt); 
                K = zeros(size(Container.O));
                K(:,:,1) = Rt.*(2-L);
                K(:,:,2) = Gt.*(2-L);
                K(:,:,3) = Bt.*(2-L);
                Container.I_BG = im2uint8(K);
                %% TEST SEGMENTAZIONE
                segmentation;
                logic_fuse = false(size(obj_chess(1).color_mask));
                for l=1:size(obj_chess,2)
                    logic_fuse = logic_fuse | obj_chess(l).color_mask;
                    logic_fuse = logic_fuse | obj_chess(l).inv_color_mask;
                end
                logic_fuse = imclose(logic_fuse,strel('square',15));
                full_image = false(size(orig,1),size(orig,2));
                Container.BB = round(Container.BB);
                full_image(Container.BB(2):Container.BB(2)+Container.BB(4),Container.BB(1):Container.BB(1)+Container.BB(3))= logic_fuse;
                imwrite(im2double(full_image),[path_3,'seg.jpg']);
                %% Identifica colori usando Background, Error Check
                Container.obj_chess = error_check(obj_chess(1), obj_chess(2), obj_chess(3),  obj_chess(4));
                %% Plot Results
                Container.fuse = show_result(Container, true, path_3);
                save([path_3,'Container.mat'],'Container');
                close all;
            end
        end
    end
end


for ff=1:size(folders,1)
    if ~(strcmp(folders(ff).name,'..') || strcmp(folders(ff).name,'.'))
        path_2 = [path,folders(ff).name,'/'];
        folders_2 = dir(path_2);
        for ff2=1:size(folders_2,1)
            if ~(strcmp(folders_2(ff2).name,'..') || strcmp(folders_2(ff2).name,'.'))
                path_3 = [path_2,folders_2(ff2).name,'/'];
                if ~exist([path_3,'manual_seg.JPG'], 'file')
                   continue;
                end
                if ~exist([path_3,'seg.JPG'], 'file')
                   continue;
                end
                disp([path_3,'....']);
                % auto_seg & manual_seg = bianco
                % ~auto_seg & manual_seg = verde
                % auto_seg & ~manual_seg = rosso
                % ~auto_seg & ~manual_seg = blu
                load([path_3,'bb.mat'])
                auto_seg = imcrop(imread([path_3,'seg.JPG'])==255,bb);
                manual_seg = imcrop(rgb2gray(imread([path_3,'manual_seg.JPG']))==255,bb);
                white = auto_seg & manual_seg;
                green = ~auto_seg & manual_seg;
                red = auto_seg & ~manual_seg;
                black = ~auto_seg & ~manual_seg;  
                full_image = double(cat(3,red|white,green|white,white|black));
                figure, imshow(full_image);
                imwrite(im2double(full_image),[path_3,'segmentation_comparison.jpg']);
%% TRUE POSITIVE, TRUE NEGATIVE, FALSE NEGATIVE, FALSE POSITIVE
                size_image = size(auto_seg,1)*size(auto_seg,2);
                true_positive = sum(sum(white));
                true_negative = sum(sum(black));
                false_positive = sum(sum(red));
                false_negative = sum(sum(green));
%% RATE
                precision = true_positive/(true_positive+false_positive);
                recall = true_positive/(true_positive+false_negative); 
                save_results_seg(true_positive, true_negative, false_positive, false_negative, precision, recall, path_3);
            end
        end
    end
end

for ff=1:size(folders,1)
    if ~(strcmp(folders(ff).name,'..') || strcmp(folders(ff).name,'.'))
        path_2 = [path,folders(ff).name,'/'];
        folders_2 = dir(path_2);
        for ff2=1:size(folders_2,1)
            if ~(strcmp(folders_2(ff2).name,'..') || strcmp(folders_2(ff2).name,'.'))
                path_3 = [path_2,folders_2(ff2).name,'/'];
                if ~exist([path_3,'Container.mat'], 'file')
                   continue;
                end
                load([path_3,'bb.mat'])
                load([path_3,'Container.mat']);
                disp([path_3,'....']);
                fuse = Container.fuse;
                modeling;
                [left_center_axis, right_center_axis, mid_center_axis] = generate_central_axis(Container);
                pointsArray = calculate_corners(Container, left_center_axis, right_center_axis, mid_center_axis);
                x_points = cell2mat(pointsArray.x_points); x_points = x_points(:)+bb(1);
                y_points = cell2mat(pointsArray.y_points); y_points = y_points(:)+bb(2);
                orig = imread([path_3,'fore.JPG']);
                image_gray = rgb2gray(orig);
                if ~exist([path_3,'manual_points_image.png'], 'file')
                   imwrite(im2double(image_gray),[path_3,'manual_points_image.png'],'png');  
                end
                for i=1:size(x_points)
                    image_gray = insertMarker(image_gray,[x_points(i), y_points(i)],'size',2);
                end
                imwrite(im2double(image_gray),[path_3,'points_image.png'],'png');
                save([path_3,'x_points.mat'],'x_points');
                save([path_3,'y_points.mat'],'y_points');
            end
        end
    end
end

for ff=1:size(folders,1)
    if ~(strcmp(folders(ff).name,'..') || strcmp(folders(ff).name,'.'))
        path_2 = [path,folders(ff).name,'/'];
        folders_2 = dir(path_2);
        for ff2=1:size(folders_2,1)
            if ~(strcmp(folders_2(ff2).name,'..') || strcmp(folders_2(ff2).name,'.'))
                path_3 = [path_2,folders_2(ff2).name,'/'];
                if ~exist([path_3,'manual_points_image.png'], 'file') || ...
                        ~exist([path_3,'points_image.png'], 'file') || ...
                            ~exist([path_3,'y_points.mat'], 'file') || ...
                                ~exist([path_3,'x_points.mat'], 'file')
                   continue;
                end
                load([path_3,'y_points.mat'])
                load([path_3,'x_points.mat']);
                manual_image = imread([path_3,'manual_points_image.png']);
                logical_manual_image = manual_image(:,:,1) == 255 & manual_image(:,:,2) == 0 & manual_image(:,:,3) == 0;
                props = regionprops(bwconncomp(logical_manual_image,8),'Centroid');
                props = struct2cell(props); props_x = []; props_y = [];
                for j=1:size(props,2)
                    tmp_props = props{j};
                    props_x = [props_x; tmp_props(1)];
                    props_y = [props_y; tmp_props(2)];
                end
                dist = 0;
                for i=1:size(x_points,1)
                    X = x_points(i); Y = y_points(i);
                    D = sqrt((props_x-X).^2 + (props_y-Y).^2);
                    [val, idx] = min(D);
                    if val <= 1.2
                        val = 0;
                    end
                    dist = dist + val;
                    props_x(idx) = [];
                    props_y(idx) = [];
                end
                dist = dist/size(x_points,1);
            end
        end
    end
end




