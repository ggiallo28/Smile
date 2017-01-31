%% Separa i riflessi
close all;
Container = split_refletions(Container, fuse, checker_vector); % O(numero viste)
%% Creazione Maschere per separazione riflessi
Container = generate_mask(Container);%obj_chess, label, order, positions, types, size(fuse));
%% Calcolo convexhull dei riflessi: controllare se è necessario fare sta cosa
Container = generate_bwconvhull(Container);
%% Calcolo Convexhull tessere singole
if Container.isGUI
%    cla(Container.app.UIAxes9);
end
% Modelling O(numero settori angolari)
for i=1:size(Container.obj_chess,1)
    if ( ~Container.obj_chess(i).isEmpty )
        for j=1:size(Container.obj_chess(i).chess,2)
            cut_x = Container.obj_chess(i).bbox_x(j,:);
            cut_y = Container.obj_chess(i).bbox_y(j,:); 
            Container.obj_chess(i).chess(j).ch_mask = false(size(fuse,1),size(fuse,2));
            mask = imdilate(bwconvhull(Container.obj_chess(i).chess(j).mask),strel('square',3));
            Container.obj_chess(i).chess(j).ch_mask(cut_y(1):cut_y(2),cut_x(1):cut_x(2)) =...
                mask(cut_y(1):cut_y(2),cut_x(1):cut_x(2));  
            if Container.isGUI
%               imshowpair(Container.obj_chess(i).chess(j).ch_mask,rgb2gray(Container.I),'falsecolor','Parent',Container.app.UIAxes9);
%               drawnow;
            else    
               figure, imshowpair(Container.obj_chess(i).chess(j).ch_mask,rgb2gray(Container.I),'falsecolor');
            end
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

