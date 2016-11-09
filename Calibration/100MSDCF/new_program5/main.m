close all; clear all; clc;
figure,imshow(imread('checkerboard.jpg'));
checker_vector = reshape([[0,0,0;255,0,255];[0,0,0;0,255,255];[0,0,0;255,255,0];[255,255,255;255,0,0];[255,255,255;0,255,0];[255,255,255;0,0,255]],[2,6,3]);
checker_center = [0.5*size(checker_vector,2),0.5*size(checker_vector,2)+1];
path = 'foto/';
name = ['DSC0014',num2str(6)];
orig = imread([path,name,'.JPG']);
orig_bg = imread([path,'DSC0014',num2str(9),'.JPG']);
%% Normalizzazione
[I, I_BG, O, O_BG, BB] = normalize_image(orig, orig_bg);
%% Parametri
confidence = 2.8;
[r,c] = size(I); Threshold = round(r*c/7000); % Soglia dimensione blob normalizzata alla dimensione dell'immagine
mpd = 15;
windowSize = 6;
op_th = 15;
if exist([path,name,'.mat'], 'file') == 2
    load([path,name,'.mat'])
end
%% Background Subtraction
[bw, inImg] = subtract_background(I, I_BG, confidence, Threshold);
%% Histogram Peak Finding
Rp = inImg(:,:,1); Rp = Rp(bw);
Gp = inImg(:,:,2); Gp = Gp(bw);
Bp = inImg(:,:,3); Bp = Bp(bw);
RGB = cat(3,Rp,Gp,Bp);
hsv = rgb2hsv(RGB);
[output_peak, output_minima_low, output_minima_high, output_minima_mid, hist_size] =...
    findlocalminima(hsv(:,:,1),mpd,windowSize,0,1);
% Non sono in grado d distinguere tra il viola e il rosso, tra l'azzurro e il blu, quindi 4 cluster invece che 6
figure, imshow(imread('hsv.jpg'));
%% Segmentazione RGB, hist_size, output_minima_mid, Threshold, band
obj_red = createMask(I, bw, hist_size, output_minima_mid, Threshold, 'red', op_th);
obj_green = createMask(I, bw, hist_size, output_minima_mid, Threshold, 'green', op_th);
obj_blue = createMask(I, bw, hist_size, output_minima_mid, Threshold, 'blue', op_th);
obj_yellow = createMask(I, bw, hist_size, output_minima_mid, Threshold, 'yellow', op_th);
% Se il contorno di sopra o sotto è troppo piccolo usare bbox
% Per eliminare merde laterali puoi calcolare la maskera con la fit square,
% erodere e risegmentare
% dilatare un pochino il bbox perchè nella fit square può succedere che si
% vada a finire fuori
%% Identifica colori usando Background, Error Check
obj_chess = error_check(obj_red, obj_green, obj_blue, obj_yellow);
%% Plot Results
fuse = show_result(obj_chess, I);
%% Separa i riflessi
[order, label, positions, types, obj_chess] = split_refletions(obj_chess, fuse, checker_vector);
%% Creazione Maschere per separazione riflessi
[maskC, isC, maskL1, isL1, maskL2, isL2, maskR1, isR1, maskR2, isR2] =...
    generate_mask(obj_chess, label, order, positions, types, size(fuse));
%% Calcolo convexhull dei riflessi: controllare se è necessario fare sta cosa
[maskCI, maskL1I, maskL2I, maskR1I, maskR2I] =...
    generate_bwconvhull(maskC, isC, maskL1, isL1, maskL2, isL2, maskR1, isR1, maskR2, isR2);
%% Calcolo Convexhull tessere singole
for i=1:size(obj_chess,1)
    if ( ~obj_chess(i).isEmpty )
        for j=1:size(obj_chess(i).chess,2)
            cut_x = obj_chess(i).bbox_x(j,:);
            cut_y = obj_chess(i).bbox_y(j,:); 
            obj_chess(i).chess(j).ch_mask = false(size(fuse,1),size(fuse,2));
            mask = imdilate(bwconvhull(obj_chess(i).chess(j).mask),strel('square',3));
            obj_chess(i).chess(j).ch_mask(cut_y(1):cut_y(2),cut_x(1):cut_x(2)) =...
                mask(cut_y(1):cut_y(2),cut_x(1):cut_x(2));     
            figure, imshowpair(obj_chess(i).chess(j).ch_mask,rgb2gray(I),'falsecolor');
        end
    end
end
close all;
%% Full Mask
FullMask  = false(size(fuse,1),size(fuse,2));
if ( isL2 )
    FullMask = FullMask | maskL2I;
end
if ( isL1 )
    FullMask = FullMask | maskL1I;
end
if ( isC )
    FullMask = FullMask | maskCI;
end
if ( isR1 )
    FullMask = FullMask | maskR1I;
end
if ( isR2 )
    FullMask = FullMask | maskR2I;
end
%% Fissa il centro degli assi
[left_center_axis, right_center_axis, mid_center_axis] =...
    generate_central_axis(I, obj_chess, FullMask, maskCI, isC, maskL1I, isL1, maskL2I, isL2, maskR1I, isR1, maskR2I, isR2);
%% Corner
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
bw = im2bw(img_filtered,graythresh(img_filtered.*FullMask)); % Se si cambia la soglia qua non funge più un cazzo
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
    GRIDh_filled = filledgegaps(GRIDh.*imerode(maskI,strel('rectangle',[30,2])),30);
    CC_GRIDh = bwconncomp(GRIDh_filled,8);
    for h_count=1:CC_GRIDh.NumObjects
        [circley,circlex] = ind2sub(size(mask),CC_GRIDh.PixelIdxList{h_count}); 
         horizontalCell{i,h_count} = createLineInv(circley,circlex,size(mask));
    end

    GRIDvv = getVImage(bwV, maskI, GRIDh);
    GRID = imerode(maskI,strel('rectangle',[30,10])) & ~mask & ~imdilate(gapGRIDLeft,strel('disk',dilateLeft)) & ~imdilate(gapGRIDRight,strel('disk',dilateRight));% | GRIDvv;
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
    if usu > 0 && toTakeRight
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
%    [imagePoints,boardSize] = detectCheckerboardPoints(I.*repmat(uint8(maskI),1,1,3));
% Orizzontale 
    R = im2double(I(:,:,1)); G = im2double(I(:,:,2)); B = im2double(I(:,:,3));
    GRIDv2 = imdilate(image_line,strel('disk',10));
    MINRGB = min(R,G); MINRGB = min(MINRGB,B);
    MAXRGB = max(R,G); MAXRGB = max(MAXRGB,B);
    gray = 0.5*MINRGB+0.5*MAXRGB;
    sep_mask = GRIDv2.*imdilate(GRIDh,strel('disk',10));
    sep_squares = gray.*sep_mask;
    sep_squares = sep_squares.*bwareaopen(sep_squares>0,10);
    CC = bwconncomp(sep_squares>0,8);
%     for l=1:CC.NumObjects
%         if(isFound(imagePoints,CC.PixelIdxList{l},size(sep_squares)))
%             sep_squares(CC.PixelIdxList{l}) = 0;
%         end
%     end
%     CC = bwconncomp(sep_squares>0,8);
    sss_image = zeros(size(sep_squares,1),size(sep_squares,2),3); P = [];
    for l=1:CC.NumObjects
        [cutR,cutC] = ind2sub(size(GRIDv),CC.PixelIdxList{l});
        tmp = sep_squares(min(cutR):max(cutR),min(cutC):max(cutC));
        tt = 0.3*size(tmp,1)*size(tmp,2);
        condition = true;
        th = 0.9; store = cell(2,1);
        dot = findDots(tmp);
        %figure, imshow(tmp), hold on, scatter(dot(:,1),dot(:,2))
        P = [P;min(cutR)+dot(:,2),min(cutC)+dot(:,1)];      
        if(isempty(dot))
            while(condition && th>0)
                bw_tmp = im2bw(tmp,th);
                th = th-0.01;
                CC_tmp = bwconncomp(bw_tmp,4);        
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
                if(CC_tmp.NumObjects == 2)                  
                   store(1) = CC_tmp.PixelIdxList(1);
                   store(2) = CC_tmp.PixelIdxList(2);          
                end 
            end
                [D(5,:), idd] = sort(D(5,:));
                D(1,:) = D(1,idd);
                D(2,:) = D(2,idd);
                D(3,:) = D(3,idd);
                D(4,:) = D(4,idd);
                XX = round(mean([D(1,1:5),D(3,1:5)]));
                YY = round(mean([D(2,1:5),D(4,1:5)]));
                %figure, imshow(bw_tmp), hold on, scatter(XX,YY);
                sss_image(min(cutR):max(cutR),min(cutC):max(cutC),1) = bw_tmp;
                sss_image(min(cutR)+YY,min(cutC)+XX,2) = 1;
                P = [P;min(cutR)+YY,min(cutC)+XX];        
        end
    end
    figure, imshow(I), hold on, scatter(P(:,2),P(:,1));
    [y_rj_matrix, x_cj_matrix] = orderPoints(P, horizontalCell(i,:), verticalCell(i,:), maskI);
    pointsArray.x_points{i} = x_cj_matrix;
    pointsArray.y_points{i} = y_rj_matrix;
    figure, imshow(I), hold on, scatter(x_cj_matrix(:),y_rj_matrix(:));
end
%% Divido i punti per ogni componente
% Scansioniamo le tessere, facciamo intersezione con le maschere di ogni
% singola tessera, i punti che ricadono dell'intersezione vengono assegnati
% alla tessera associata alla specifica maschera.
% Per fare questo usiamo il vettore delle label e order creato
% precedentemente per poter facilemente recuperare gli indici di riga e
% colonna associati ad ogni singolo elemento.
for k=1:5
    switch(k)
        case 1
            if ( ~isC )
                continue;
            end
            mask = maskCI;
            idx = find(strcmp(label(3,:),positions(2)));
        case 2
            if ( ~isL1 )
                continue;
            end
            mask = maskL1I;
            idx = find(strcmp(label(3,:),positions(1)) & strcmp(label(4,:),types(1)));
        case 3
            if ( ~isL2 )
                continue;
            end
            mask = maskL2I;
            idx = find(strcmp(label(3,:),positions(1)) & strcmp(label(4,:),types(2)));
        case 4
            if ( ~isR1 )
                continue;
            end
            mask = maskR1I;
            idx = find(strcmp(label(3,:),positions(3)) & strcmp(label(4,:),types(1)));
        case 5
            if ( ~isR2 )
                continue;
            end
            mask = maskR2I;
            idx = find(strcmp(label(3,:),positions(3)) & strcmp(label(4,:),types(2)));
    end
    for i=1:size(idx,2)
        idx_chess_vector = order(1,idx(i));
        idx_color_chess = order(2,idx(i));
        x_cj_matrix = pointsArray.x_points{k};
        y_rj_matrix = pointsArray.y_points{k};
        ch_mask = obj_chess(idx_chess_vector).chess(idx_color_chess).ch_mask;
        w = find(sum(ch_mask) ~=0);
        w_min_x = min(w); w_max_y = max(w);
        check_matrix = false(size(x_cj_matrix));
        for j=1:size(x_cj_matrix,1)
            for q=1:size(x_cj_matrix,2)
                if(x_cj_matrix(j,q)>=(w_min_x-10) && x_cj_matrix(j,q)<=(w_max_y+10)) % 10 pixel di tolleranza: Parametro, speriamo vada bene
                    check_matrix(j,q) = true;
                end
            end
        end
        num_cols = mode(sum(check_matrix,2));
        sum_cols = sum(check_matrix);
        check_matrix = false(size(check_matrix)); % Fix, dobbiamo sempre prendere tutta una colonna, scegliamo il valore più frequente e completiamo
        for j=1:num_cols
            id = find(sum_cols ==max(sum_cols), 1 ); % prendi solo il primo
            check_matrix(:,id) = 1;
            sum_cols(id) = 0;
        end    
        y_rj_matrix_copy = y_rj_matrix.*check_matrix;
        x_cj_matrix_copy = x_cj_matrix.*check_matrix;
        id = find(y_rj_matrix_copy(1,:) == 0);
        y_rj_matrix_copy(:,id) = [];
        x_cj_matrix_copy(:,id) = [];
        
        obj_chess(idx_chess_vector).chess(idx_color_chess).intersections_x = x_cj_matrix_copy;
        obj_chess(idx_chess_vector).chess(idx_color_chess).intersections_y = y_rj_matrix_copy;   
        imshow([I;repmat(255.*ch_mask,1,1,3)]); hold on;
        scatter(x_cj_matrix_copy(:),y_rj_matrix_copy(:));
    end    
end
%% Padding delle matrici e allineamento
% Per ogni combinazione colore/sfondo controlliamo quant'e il massimo
% numero di colonne, aggiungiamo zeri alle matrici che hanno un numero
% inferiore di colonne
for l=1:size(obj_chess,1)
    if ( ~obj_chess(l).isEmpty )
        v = zeros(1,size(obj_chess(l).chess,2));
        for q=1:size(obj_chess(l).chess,2)
           v(q) = size(obj_chess(l).chess(q).intersections_x,2); 
        end
        max_v = max(v);
        for q=1:size(obj_chess(l).chess,2)
            num_col = size(obj_chess(l).chess(q).intersections_x,2);
            num_row = size(obj_chess(l).chess(q).intersections_x,1);
            if(num_col<max_v)
                obj_chess(l).chess(q).intersections_x = ...
                    [obj_chess(l).chess(q).intersections_x, zeros(num_row,max_v-num_col)];
                obj_chess(l).chess(q).intersections_y = ...
                    [obj_chess(l).chess(q).intersections_y, zeros(num_row,max_v-num_col)];
            end
            % Allineo a sinistra o a destra in funzione del fatto che ho un elemento a sinistra o a destra.
            % Ciò mi è utile per tenere tutti i punti compatti verso il centro
            % di ciò che è visibile sulla checherboard.
            % Per capire se ho un elemento a sinistra controllo che a sinistra
            % ci sia una griglia dello stesso tipo e nella stessa posizione
            % della corrente.
            % Per capire se ho un elemento a destra faccio la medesima cosa
            % guardando a destra.
            % Se abbiamo una griglia sinistra e una griglia destra sto vedendo
            % tutta la griglia corrente, quindi non ha importanza in che
            % direzione shifto.
            % Una volta che gli elementi sono stati allineati in questo modo è
            % sufficiente flippare le matrici associate ai riflessi primari,
            % unico caso in cui abbiamo un'inversione.
            % Gli indici di riga/colonna rappresentano le associazioni, abbiamo
            % degli zeri dove l'associazione non esiste poichè il punto non
            % risulta visibile.
            % LABEL: Colore Foreground, Colore Background, Posizione, Tipo
            % ORDER: Indice nel vettore obj_chess, Indice Chess, Centroide
            idl = find(order(1,:)==l);
            idq = find(order(2,:)==q);
            idx_order = intersect(idl,idq);
            if(hasLeft(label,idx_order))
              obj_chess(l).chess(q).intersections_x = ...
                  allignleftdouble(obj_chess(l).chess(q).intersections_x);
              obj_chess(l).chess(q).intersections_y = ...
                  allignleftdouble(obj_chess(l).chess(q).intersections_y);
            end
            if(hasRight(label,idx_order))
               obj_chess(l).chess(q).intersections_x = ...
                   allignrightdouble(obj_chess(l).chess(q).intersections_x);
               obj_chess(l).chess(q).intersections_y = ...
                   allignrightdouble(obj_chess(l).chess(q).intersections_y);
            end
            % positions = [{'Left'} {'Center'} {'Right'}];
            % types = [{'Primary'} {'Secondary'} {'Real'}];
            if(strcmp(label(4,idx_order),types(1)))
               obj_chess(l).chess(q).intersections_x = ...
                   fliplr(obj_chess(l).chess(q).intersections_x);
               obj_chess(l).chess(q).intersections_y = ...
                   fliplr(obj_chess(l).chess(q).intersections_y);
            end
        end
    end
end

%% Mapping sulla superficie del Manifold
axis_line = reshape([192,192,192;192,192,192],2,1,3);
checher_vector_with_axis = uint8([checker_vector(:,1:3,:),axis_line,checker_vector(:,4:6,:)]); 
figure, imshow(checher_vector_with_axis);

% La checker è suddivisa verticalmente in 4 quadrati, ciò significa 5 punti compresi quelli di intersezione
% Ogni quadrante rappresenta 45°, quindi ad ogni punto corrisponde un angolo di 9°
for l=1:size(obj_chess,1)
    if ( ~obj_chess(l).isEmpty )
        for i=1:size(obj_chess(l).chess(1).intersections_x,1)
            for j=1:size(obj_chess(l).chess(1).intersections_x,2)
                imshow(I); hold on;
                vect_x = [];
                vect_y = [];
                for k =1:size(obj_chess(l).chess,2)
                    axis_distance = axisdistance(obj_chess(l).name,obj_chess(l).chess(k).background,checher_vector_with_axis);
                    adjusted_distance = axis_distance-1*sign(axis_distance); % Se positivo togliamo 1, se negativo aggiungiamo 1, le tessere gialle e rosse stanno a distanza zero.
                    if(sign(axis_distance) >0)
                        adjusted_offset = (j-1)*11.25;
                    else
                        adjusted_offset = (j-size(obj_chess(l).chess(k).intersections_x,2))*11.25;
                    end
                    angle = adjusted_distance*45 + adjusted_offset
                    h = i*2.6
                    vect_x = [vect_x, obj_chess(l).chess(k).intersections_x(i,j)];
                    vect_y = [vect_y, obj_chess(l).chess(k).intersections_y(i,j)];
                end
                scatter(vect_x(:),vect_y(:)); hold off
                pause
            end
        end
    end
end
%% Calcolo angolo degli specchi
Pleft = [];
Pright = [];
for l=1:size(obj_chess,1)
    if ( ~obj_chess(l).isEmpty )
        for i=1:size(obj_chess(l).chess(1).intersections_x,1)
            for j=1:size(obj_chess(l).chess(1).intersections_x,2)
                vect_x = zeros(2,3);
                vect_y = zeros(2,3);
                idk = 1;
                for k =1:size(obj_chess(l).chess,2)           
                    if (~strcmp(obj_chess(l).chess(k).type,types(2)))
                        vect_x(1,idk) = obj_chess(l).chess(k).intersections_x(i,j); vect_x(2,idk) = k;
                        vect_y(1,idk) = obj_chess(l).chess(k).intersections_y(i,j); vect_y(2,idk) = k;
                        idk = idk +1;
                    end                                        
                end
                vect_xleft = 0; vect_yleft = 0; vect_xmid = 0; vect_ymid = 0; vect_xright = 0; vect_yright = 0;
                if idk > 2
                    for k=1:size(vect_x,2)
                        idk = find(order(2,:)==vect_x(2,k) & order(1,:) == l);
                        if ~isempty(idk)
                            pos = label(3,idk);
                            type = label(4,idk);
                            assert(~strcmp(type,types(2)));
                            if (strcmp(pos,positions(1)))
                                vect_xleft = vect_x(1,k);
                                vect_yleft = vect_y(1,k);
                            end
                            if (strcmp(pos,positions(2)))
                                vect_xmid = vect_x(1,k);
                                vect_ymid = vect_y(1,k);
                            end
                            if (strcmp(pos,positions(3)))
                                vect_xright = vect_x(1,k);
                                vect_yright = vect_y(1,k);
                            end
                        end
                    end
                    if(~(vect_xleft == 0 || vect_yleft == 0 || vect_xmid == 0 || vect_ymid == 0))
                        Pleft = [Pleft;[vect_xmid,vect_ymid,vect_xleft,vect_yleft]];
                    end
                    if(~(vect_xright == 0 || vect_yright == 0 || vect_xmid == 0 || vect_ymid == 0))
                        Pright = [Pright;[vect_xmid,vect_ymid,vect_xright,vect_yright]];
                    end
                end
            end
        end
    end
end
Pleft = unique(Pleft,'rows');
Pright = unique(Pright,'rows');
imshow(I); hold on;
scatter([Pleft(:,1);Pleft(:,3)],[Pleft(:,2);Pleft(:,4)]);
pause
imshow(I); hold on;
scatter([Pright(:,1);Pright(:,3)],[Pright(:,2);Pright(:,4)]);
%% Ora come associo i punti se levo il cilindro? Con il codice che già hanno ma usando la griglia.
%% TODO: Il prof ha detto di stimare la differenza di colore negli histogrammi del rosso al centro e del rosso a destra

%% Come stimo l'angolo degli specchi?

%% Salvataggio Setup
save([path,name,'.mat'],'confidence','Threshold','mpd','windowSize','op_th');