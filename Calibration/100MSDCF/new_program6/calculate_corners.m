function pointsArray = calculate_corners(Container, left_center_axis, right_center_axis, mid_center_axis)
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

    color_filtered = abs(imfilter(im2double(TheColors+left_center_axis+right_center_axis+mid_center_axis),fspecial('prewit')));
    key_filtered = absimfilter(KEY,fspecial('prewit'));
    orig_filtered = absimfilter(ORIG,fspecial('prewit'));
    img_filtered = (key_filtered+orig_filtered+color_filtered)./3;
    bw = im2bw(img_filtered,graythresh(img_filtered.*FullMask)); % Se si cambia la soglia qua non funge pi� un cazzo
    bw = imclose(bw,strel('square',5));
    bwV = bwareaopen(bw, 100); % Parametro
    pointsArray = struct();
    horizontalCell = cell(5,20);
    verticalCell = cell(5,20);
    for i=1:5 % Non sei indipendente dal numero di riflessi
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
                toTakeLeft = true;
                toTakeRight = true;
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
                toTakeRight = true;
                if ( isL2 )
                    dilateLeft = 15;
                    toTakeLeft = false;
                else
                    dilateLeft = 12;
                    toTakeLeft = true;
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
                toTakeLeft = true;
                if ( isR2 )
                   dilateRight = 15; 
                   toTakeRight  = false;
                else
                   dilateRight = 12; 
                   toTakeRight  = true;  
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
                dilateLeft = 15; % Dobbiamo dilatare di una quantit� diversa a causa delle occlusioni tra checkerboard che introducono contorni inesistenti
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
            % quelle che sono fittate con pi� punti e che hanno coefficente
            % angolare pi� piccolo poich� orizontali, oppure semplicemente se
            % so di pi� erodi maggiormente sopra e sotto
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
        % le rette che tagliano il bordo pi� esterno
        GRIDv = abs(imfilter(im2double(GRID),[-1 0 1]));

        GRIDh = imdilate(GRIDh,strel('disk',10)).*maskI;
        stats = cell2mat(struct2cell(regionprops(GRIDh,'Area')));
        GRIDh = imopen(bwareaopen(GRIDh,round(0.1*(mean(stats)-std(stats)))),strel('disk',10));

        GRIDv = GRIDv - GRIDh;
        GRIDv(GRIDv<0) = 0;
        stats = cell2mat(struct2cell(regionprops(bwconncomp(GRIDv,8),'Area')));
        GRIDv = bwareaopen(GRIDv,50); % Aggiustare questa soglia

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
    % Provo con la funzione di matlab, se questa non rintraccia alcuni punti
    % allora provo con il mio metodo
        [imagePoints,~] = detectCheckerboardPoints(I.*repmat(uint8(maskI),1,1,3));
        imshow(I); hold on; plot(imagePoints(:,1), imagePoints(:,2), 'ro')
    % Orizzontale
        R = im2double(I(:,:,1)); G = im2double(I(:,:,2)); B = im2double(I(:,:,3));
        GRIDv2 = imdilate(image_line,strel('disk',10));
        MINRGB = min(R,G); MINRGB = min(MINRGB,B);
        MAXRGB = max(R,G); MAXRGB = max(MAXRGB,B);
        gray = 0.5*MINRGB+0.5*MAXRGB;
        gray = imadjust(im2double(rgb2gray(I)));
        sep_mask = GRIDv2.*imdilate(GRIDh,strel('disk',15));
%% Aggiunta
        CC = bwconncomp(sep_mask,8);
        size_cc = CC.NumObjects; condition = true;
        dil = 1; sep_mask_tmp = sep_mask;
        while condition
            sep_mask_tmp = imdilate(sep_mask,strel('square',dil));
            CC = bwconncomp(sep_mask_tmp,8);
            if CC.NumObjects < size_cc
                break;
            end
            dil = dil + 2;
        end
        sep_mask = imdilate(sep_mask,strel('square',floor(dil*0.2)));
%% Aggiunta
        sep_squares = gray.*sep_mask;
        sep_squares = sep_squares.*bwareaopen(sep_squares>0,10);
        CC = bwconncomp(sep_squares>0,8);
        mask_sep_squares = sep_squares>0; ii = []; points_pos_image = false(size(mask_sep_squares));
        for points_i=1:size(imagePoints,1)
            points_pos = logical(rgb2gray(insertMarker(zeros(size(mask_sep_squares)),imagePoints(points_i,:),'size',1)));
            if(sum(sum(points_pos.*mask_sep_squares))==0)
                ii = [ii,points_i];
            end 
            points_pos_image = (points_pos_image | points_pos).*mask_sep_squares;
        end
        imagePoints(ii,:) = [];
        for l=1:CC.NumObjects
            points_pos = false(size(mask_sep_squares));
            points_pos(CC.PixelIdxList{l}) = 1;
            if(sum(sum(points_pos_image.*points_pos)) > 0)
                sep_squares(CC.PixelIdxList{l}) = 0;
            end
        end
        CC = bwconncomp(sep_squares>0,8);
    %     for l=1:CC.NumObjects
    %         if(isFound(imagePoints,CC.PixelIdxList{l},size(sep_squares)))
    %             sep_squares(CC.PixelIdxList{l}) = 0;
    %         end
    %     end
    %     CC = bwconncomp(sep_squares>0,8);
        sss_image = zeros(size(sep_squares,1),size(sep_squares,2),3); P = [];
%       deb = zeros(size(sep_squares));
        for l=1:CC.NumObjects
            [cutR,cutC] = ind2sub(size(GRIDv),CC.PixelIdxList{l});
            tmp = sep_squares(min(cutR):max(cutR),min(cutC):max(cutC));
%            deb(min(cutR):max(cutR),min(cutC):max(cutC)) = tmp;
            tt = 0.3*size(tmp,1)*size(tmp,2);
            condition = true;
            th = 0.9; store = cell(2,1);
            dot = findDots(tmp);
            %figure, imshow(tmp), hold on, scatter(dot(:,1),dot(:,2))
            P = [P;min(cutR)+dot(:,2),min(cutC)+dot(:,1)];      
            if(isempty(dot))
                figure
                while(condition && th>0)
                    bw_tmp = im2bw(tmp,th);
                    th = th-0.005;
                    CC_tmp = bwconncomp(bw_tmp,4);
%                     imshow(bw_tmp);
%                     pause(0.1)
%                     th
                    if(CC_tmp.NumObjects ==1 && size(CC_tmp.PixelIdxList{1},1)>tt)
                        condition = false;
                        [blob1y,blob1x] = ind2sub(size(tmp),store{1});
                        [blob2y,blob2x] = ind2sub(size(tmp),store{2});
                        D = zeros(5,size(blob1y,1)*(size(blob2y,1)-1)); idist = 1;
                        for k=1:size(blob1y,1)
                            for j=1:size(blob2y,1)
                                D(1,idist) = blob1x(k);
                                D(2,idist) = blob1y(k);
                                D(3,idist) = blob2x(j);
                                D(4,idist) = blob2y(j);
                                D(5,idist) = pdist([blob1x(k) blob1y(k); blob2x(j) blob2y(j)],'euclidean');
                                idist = idist+1;
                            end
                        end
                        if(isempty(D))
                             condition = true;
                        end
                    end
                    if(CC_tmp.NumObjects >= 2) 
%                        CC_tmp
                       store(1) = CC_tmp.PixelIdxList(1);
                       store(2) = CC_tmp.PixelIdxList(2);          
                    end 
                end
                %% INVECE DIFARE STA MANFRINA USARE corner(bw_tmp,'Harris'
                if(isempty(D))
                    YY = size(bw_tmp,1)*0.5;
                    XX = size(bw_tmp,2)*0.5;
                else
                    [D(5,:), idd] = sort(D(5,:));
                    D(1,:) = D(1,idd);
                    D(2,:) = D(2,idd);
                    D(3,:) = D(3,idd);
                    D(4,:) = D(4,idd);
                    if size(D,2)>=5
                        XX = round(mean([D(1,1:5),D(3,1:5)]));
                        YY = round(mean([D(2,1:5),D(4,1:5)]));
                    else
                        XX = round(mean([D(1,:),D(3,:)]));
                        YY = round(mean([D(2,:),D(4,:)]));
                    end
                        %figure, imshow(bw_tmp), hold on, scatter(XX,YY);
                        sss_image(min(cutR):max(cutR),min(cutC):max(cutC),1) = bw_tmp;
                        sss_image(min(cutR)+YY,min(cutC)+XX,2) = 1;
                        P = [P;min(cutR)+YY,min(cutC)+XX];
                end
            end
        end
        if ~isempty(imagePoints)
            P = [P;[imagePoints(:,2),imagePoints(:,1)]];
        end
        figure, imshow(I), hold on, scatter(P(:,2),P(:,1));
        [y_rj_matrix, x_cj_matrix] = orderPoints(P, horizontalCell(i,:), verticalCell(i,:), maskI);
        pointsArray.x_points{i} = x_cj_matrix;
        pointsArray.y_points{i} = y_rj_matrix;
        figure, imshow(I), hold on, scatter(x_cj_matrix(:),y_rj_matrix(:));
%         figure, imshow(I), hold on;
%         for k=1:size(x_cj_matrix,1)
%             for j=1:size(x_cj_matrix,2)
%                 scatter(x_cj_matrix(k,j),y_rj_matrix(k,j))
%                 pause
%             end
%         end
    end
end