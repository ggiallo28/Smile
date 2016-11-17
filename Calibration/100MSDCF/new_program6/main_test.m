close all; clear all; clc;
path = '../foto/Evaluation/';
folders = dir(path);
for ff=1:size(folders,1)
    if ~(strcmp(folders(ff).name,'..') || strcmp(folders(ff).name,'.'))
        path_2 = [path,folders(ff).name,'/'];
        folders_2 = dir(path_2);
        for ff2=1:size(folders_2,1)
            if ~(strcmp(folders_2(ff2).name,'..') || strcmp(folders_2(ff2).name,'.'))
                path_3 = [path_2,folders_2(ff2).name,'/'];
                if exist([path_3,'MySavedPlot.png'], 'file') == 2
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
                fuse = show_result(Container, true, path_3);
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
                % auto_seg & manual_seg = bianco
                % ~auto_seg & manual_seg = verde
                % auto_seg & ~manual_seg = rosso
                % ~auto_seg & ~manual_seg = nero
                load([path_3,'bb.mat'])
                auto_seg = imcrop(imread([path_3,'seg.JPG'])==255,bb);
                manual_seg = imcrop(rgb2gray(imread([path_3,'manual_seg.JPG']))==255,bb);
                white = auto_seg & manual_seg;
                green = ~auto_seg & manual_seg;
                red = auto_seg & ~manual_seg;
                black = ~auto_seg & ~manual_seg;  
                full_image = double(cat(3,red|white,green|white,white));
                imshow(full_image);
%% TRUE POSITIVE, TRUE NEGATIVE, FALSE NEGATIVE, FALSE POSITIVE
                size_image = size(auto_seg,1)*size(auto_seg,2);
                true_positive = sum(sum(white));
                true_negative = sum(sum(black));
                false_positive = sum(sum(red));
                false_negative = sum(sum(green));
%% RATE
                fpr = false_positive/(false_positive+true_negative);
                precision = true_positive/(true_positive+false_positive);
                recall = true_positive/(true_positive+false_negative); 
            end
        end
    end
end
