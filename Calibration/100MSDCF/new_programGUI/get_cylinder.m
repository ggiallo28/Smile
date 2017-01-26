function pointsArray = get_cylinder(Container)
%% INIT
    I = Container.I;
    O = Container.O;
    obj_chess = Container.obj_chess;
    FullMask = Container.FullMask;
    maskCI = Container.maskCI;
    isC = Container.isC;
    maskL1I = Container.maskL1I;
    isL1 = Container.isL1;
    maskL2I = Container.maskL2I;
    isL2 = Container.isL2;
    maskR1I = Container.maskR1I;
    isR1 = Container.isR1;
    maskR2I = Container.maskR2I;
    isR2 = Container.isR2;
    label = Container.label;
    maskC = Container.maskC;
    maskL1 = Container.maskL1;
    maskL2 = Container.maskL2;
    maskR1 = Container.maskR1;
    maskR2 = Container.maskR2;
    order = Container.order;
%% LOGIC
    cymk = rgb2cmyk(I);
    Cyano = im2double(cymk(:,:,1));
    Magenta = im2double(cymk(:,:,2));
    Yellow = im2double(cymk(:,:,3));
    TheColors = Cyano+Magenta+Yellow;
    Key = im2double(cymk(:,:,4));

    KEY = (1-Key);
    ORIG = rgb2gray(im2double(O));

    color_filtered = abs(imfilter(im2double(TheColors),fspecial('prewit')));
    key_filtered = absimfilter(KEY,fspecial('prewit'));
    orig_filtered = absimfilter(ORIG,fspecial('prewit'));
    img_filtered = (key_filtered+orig_filtered+color_filtered)./3;
    bw = im2bw(img_filtered,graythresh(img_filtered.*FullMask)); % Se si cambia la soglia qua non funge più un cazzo
    bw = imclose(bw,strel('square',5));
    bwV = bwareaopen(bw, 100); % Parametro
    pointsArray = struct();
    horizontalCell = cell(5,20);
    verticalCell = cell(5,20);
    for i=1:5
        switch (i)
            case 1
                if ( ~isC )
                    continue;
                end
                % Center Real Reflection
                idl = find(strcmp(label(3,:),'Center') & strcmp(label(4,:),'Real'));
                grid_pos = 'Center';
                grid_typ = 'Real';
                maskI = maskCI;
                mask = maskC;
                dilateLeft = 12;
                dilateRight = 12;
                toTakeLeft = false;
                toTakeRight = false;
            case 2
                if ( ~isL1 )
                    continue;
                end
                % Left Primary Reflection
                idl = find(strcmp(label(3,:),'Left') & strcmp(label(4,:),'Primary'));
                grid_pos = 'Left';
                grid_typ = 'Primary';
                maskI = maskL1I;
                mask = maskL1;
                dilateRight = 12;
                toTakeRight = false;
                if ( isL2 )
                    dilateLeft = 15;
                    toTakeLeft = false;
                else
                    dilateLeft = 12;
                    toTakeLeft = false;
                end 
            case 3
                if ( ~isL2 )
                    continue;
                end
                % Left Secondary Reflection
                idl = find(strcmp(label(3,:),'Left') & strcmp(label(4,:),'Secondary'));
                grid_pos = 'Left';
                grid_typ = 'Secondary';
                maskI = maskL2I;
                mask = maskL2;
                dilateLeft = 12;
                dilateRight = 15;
                toTakeRight  = false;
                toTakeLeft = false;
            case 4
                if ( ~isR1 )
                    continue;
                end
                % Right Primary Reflection
                idl = find(strcmp(label(3,:),'Right') & strcmp(label(4,:),'Primary'));
                grid_pos = 'Right';
                grid_typ = 'Primary';
                maskI = maskR1I;
                mask = maskR1;
                dilateLeft = 12;
                toTakeLeft = false;
                if ( isR2 )
                   dilateRight = 15; 
                   toTakeRight  = false;
                else
                   dilateRight = 12; 
                   toTakeRight  = false;  
                end
            case 5
                if ( ~isR2 )
                    continue;
                end
                % Right Secondary Reflection
                idl = find(strcmp(label(3,:),'Right') & strcmp(label(4,:),'Secondary'));
                grid_pos = 'Right';
                grid_typ = 'Secondary';
                maskI = maskR2I;
                mask = maskR2;
                dilateLeft = 15; % Dobbiamo dilatare di una quantità diversa a causa delle occlusioni tra checkerboard che introducono contorni inesistenti
                dilateRight = 12;
                toTakeLeft = false;
                toTakeRight = false;
        end   
        id_obj_chess = order(1,idl);
        id_chess = order(2,idl);
        gapGRIDLeft = line2image(obj_chess(id_obj_chess(1)).chess(id_chess(1)).v_lines{1},size(mask));
        gapGRIDRight = line2image(obj_chess(id_obj_chess(end)).chess(id_chess(end)).v_lines{end},size(mask));
        [GRIDhh, GRIDh] = getHImage(grid_pos,grid_typ,obj_chess, order, label, size(maskI), maskI);
        GRIDh = filledgegaps(GRIDh.*imerode(maskI,strel('rectangle',[30,2])),30);
        condition = true;
        percent = 0.2;
        while condition
            GRIDh = bwareaopen(GRIDh,round(percent*sqrt(Container.size_small))); % Parametro
            % Togliere sta cosa e controllare di avere 4 linee, scegliere
            % quelle che sono fittate con più punti e che hanno coefficente
            % angolare più piccolo poichè orizontali, oppure semplicemente se
            % so di più erodi maggiormente sopra e sotto
            CC_GRIDh = bwconncomp(GRIDh,8);
            condition = CC_GRIDh.NumObjects ~=4; % Legato al numero di quadrati, se so 5 ho 4 linee in mezzo.
            percent = percent + 0.01;
        end
        figure, imshow(GRIDh); hold on
        for h_count=1:CC_GRIDh.NumObjects
            [circley,circlex] = ind2sub(size(mask),CC_GRIDh.PixelIdxList{h_count}); 
             horizontalCell{i,h_count} = createLine(circlex,circley);
             plot(horizontalCell{i,h_count});
        end
        
        hold off

        GRIDvv = getVImage(bwV, maskI, GRIDh);
        GRID = imerode(maskI,strel('rectangle',[30,12])) & ~mask & ~imdilate(gapGRIDLeft,strel('disk',dilateLeft)) & ~imdilate(gapGRIDRight,strel('disk',dilateRight));% | GRIDvv;
        %thinedImage = bwmorph(GRID,'thin',inf);
        % ordina i blob in base alla distanza dalla retta destra e sinistra,
        % misura la distanza media, usa questa per definire di quanto allargare
        % le rette che tagliano il bordo più esterno
        GRIDv = abs(imfilter(im2double(GRID),[-1 0 1]));

        GRIDh = imdilate(GRIDh,strel('disk',10)).*maskI;
        stats = cell2mat(struct2cell(regionprops(GRIDh,'Area')));
        GRIDh = imopen(bwareaopen(GRIDh,round(0.1*(mean(stats)-std(stats)))),strel('disk',10));

        GRIDv = GRIDv - GRIDh;
        GRIDv(GRIDv<0) = 0;
        stats = cell2mat(struct2cell(regionprops(bwconncomp(GRIDv,8),'Area')));
        GRIDv = bwareaopen(GRIDv,floor(0.8*sqrt(Container.size_big))); % Aggiustare questa soglia

        figure, imshow(GRIDv);

        props_center = regionprops(GRIDv,'Centroid'); % [X; Y]
        props_center = reshape(cell2mat(struct2cell(props_center)),2,size(props_center,1));
        props = [1:size(props_center,2);props_center];
        [props(2,:), props(1,:)] = sort(props(2,:));
        props(3,:) = props(3, props(1,:));
        CC = bwconncomp(GRIDv,8); usu_center = [zeros(1,size(props_center,2));props_center]; usu = CC.NumObjects;
        figure, imshow(I); hold on; image_line = false(size(GRIDv)); idx_line = 1;
        for j=1:size(idl,2)
            for l=1:size(obj_chess(id_obj_chess(j)).chess(id_chess(j)).v_lines_centroid,2)
                line = obj_chess(id_obj_chess(j)).chess(id_chess(j)).v_lines_centroid{l};
                lines = [];
                if ~isempty(line)
                    num_taken = 0;
                    for k=1:CC.NumObjects
                        if(usu_center(1,k)==0)
                            if(isLeft(line,usu_center(2,k),usu_center(3,k)))
                                usu_center(1,k) = 1;
                                usu = usu - 1;
                                lines = [lines; CC.PixelIdxList{k}];
                            end
                        end
                    end
                    if (l==1 && j==1 && ~toTakeLeft)
                        continue;
                    end
                    if(~isempty(lines))
                        [liney,linex] = ind2sub(size(GRIDv),lines); 
                        [linefitresult, ~] = createLineInv(liney, linex, size(GRIDv));
                        image_line = image_line | line2image(linefitresult,size(GRIDv));
                        plot(linefitresult);
                        verticalCell{i,idx_line} = linefitresult;
                        idx_line = idx_line + 1;
                    end
                end
            end
        end
        lines = [];
        if usu > 0 && toTakeRight && ~isempty(line)
            for k=1:CC.NumObjects
                if(usu_center(1,k)==0)
                    if(~isLeft(line,usu_center(2,k),usu_center(3,k)))
                        usu_center(1,k) = 1;
                        usu = usu - 1;
                        lines = [lines; CC.PixelIdxList{k}];
                    end
                end
            end
            if(~isempty(lines))
                [liney,linex] = ind2sub(size(GRIDv),lines); 
                [linefitresult, ~] = createLineInv(liney, linex, size(GRIDv));
                image_line = image_line | line2image(linefitresult,size(GRIDv));
                plot(linefitresult);
                verticalCell{i,idx_line} = linefitresult;
                idx_line = idx_line + 1;
            end
        end
    end
end