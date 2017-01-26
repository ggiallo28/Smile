function obj_chess = completeComputeChess(obj_chess, transtions, Container)
%% INIT
    warning off
    RGB = Container.I;
    Threshold = Container.Threshold;
    op_th = Container.op_th;
%% LOGIC
    for kk=1:size(obj_chess,2)
        obj = obj_chess(kk);
        if ~obj.isEmpty
            disp(['Computing ',obj.name,' ...']);
            x = obj.bbox_x;
            y = obj.bbox_y;
            BW = obj.raw_bw;
            maskedRGBImage = obj.raw_maskedRGBImage;
            if Container.isGUI
                axes_gui = Container.app.UIAxes8;
                imshow(obj.true_color_mask,'Parent',axes_gui);
            end
    %% Remove bad things
            ii = [];
            for i = 1:size(x,1)
                mask = false(size(maskedRGBImage,1),size(maskedRGBImage,2));
                mask(y(i,1):y(i,2),x(i,1):x(i,2))=true;
                color_square =  BW&mask;
                [color_square, yy, xx, isDone, isErased] = checkbadthings(color_square, Threshold,y, x(i,:), transtions, kk, i, Container, obj_chess);
                if ( isDone & ~isErased)
                    BW(y(i,1):y(i,2),x(i,1):x(i,2)) = color_square(y(i,1):y(i,2),x(i,1):x(i,2));
                    y(i,:) = yy;
                    x(i,:) = xx;
                end
                if ( isErased )
                    ii = [ii,i];
                end
            end
            y(ii,:) = [];
            x(ii,:) = [];
            obj.chess(ii) = [];
            obj.bbox_x = x;
            obj.bbox_y = y;
    %% Adjust 2 squares
            [BW, y] = twosquares_adjust(x, y, BW, RGB);
            obj.bbox_y = y;
            %% Calcolo della chessboard complementare & background
            [inv_BW, blackwhite_BW, colors, obj] = compute_dualchessboard(x, y, BW, RGB, obj, size(BW), Container);
            % Fix dei triangoli del background DA RIMUOVERE
            % CC = bwconncomp(blackwhite_BW);
            % for i = 1:size(CC.PixelIdxList,2)
            %     BW_TMP = false(size(blackwhite_BW));
            %     BW_TMP(CC.PixelIdxList{i}) = 1;
            %     s = regionprops(BW_TMP,'BoundingBox');
            %     img = imcrop(BW_TMP,s.BoundingBox); 
            %     toFlip = isTriangle(img);
            %     white = sum(sum(img));
            %     black = size(img,1)*size(img,2); 
            %     %figure, subplot(211); imshow(img);
            %     if(white<0.6*black && black>Threshold/3 && black<Threshold*14 && toFlip)
            %         img_fill = img |rot90(img,2);
            %         img_fill = imfill(img_fill,'holes');
            %         %subplot(212); imshow(img_fill);        
            %         blackwhite_BW(s.BoundingBox(2):s.BoundingBox(2)+s.BoundingBox(4),s.BoundingBox(1):s.BoundingBox(1)+s.BoundingBox(3)) = img_fill; 
            %         imshow(img_fill);
            %     end
            % end
            obj.black_white = imfill(blackwhite_BW,'holes');

            %% Separazione delle Tessere
            bwD = mask_split(BW, inv_BW, x, y, Container);

            % Fix dei triangoli della machera complementare
            % CC = bwconncomp(bwD);
            % for i = 1:size(CC.PixelIdxList,2)
            %     BW_TMP = false(size(blackwhite_BW));
            %     BW_TMP(CC.PixelIdxList{i}) = 1;
            %     s = regionprops(BW_TMP,'BoundingBox');
            %     img = imcrop(BW_TMP,s.BoundingBox);
            %     toFlip = isTriangle(img);
            %     white = sum(sum(img));
            %     black = size(img,1)*size(img,2);    
            % %   figure, subplot(211); imshow(img);
            %     if(white<0.6*black && black>500 && black<10000 && toFlip)
            %         img_fill = img |rot90(img,2);
            %         img_fill = imfill(img_fill,'holes');
            % %       subplot(212); imshow(img_fill);
            %         img_fill = img_fill(4:end-3,4:end-3);
            %         bwD((s.BoundingBox(2)+3):(s.BoundingBox(2)+s.BoundingBox(4)-3),(s.BoundingBox(1)+3):(s.BoundingBox(1)+s.BoundingBox(3)-3)) = img_fill; 
            %     end
            % end
            inv_BW = bwD & ~BW;
            obj.inv_color_mask = inv_BW;
            obj.color_mask = bwD & ~inv_BW;
            obj.color_mask = imopen(obj.color_mask,strel('square',4));

            %% Save Result
            for i = 1:size(x,1)
                chess = bwD(y(i,1):y(i,2),x(i,1):x(i,2));
                chess = bwareaopen(chess, 50); % Parametro
                height_arr = cell2mat(struct2cell(regionprops(chess,'MajorAxisLength')));
                height = mean(height_arr);
                toRemove = height_arr<height*0.7;
                CC = bwconncomp(chess);
                if(sum(toRemove)>0)
                    idx = find(toRemove == 1);
                    for j=1:size(idx,2)
                        chess(CC.PixelIdxList{idx(j)})=0;
                    end
                    CC = bwconncomp(chess);
                end
                s = regionprops(CC,'centroid');
                obj = putCenters(obj, x(i,1), y(i,1), s, i, Container);
                obj.chess(i).mask = false(size(bwD));
                obj.chess(i).mask(y(i,1):y(i,2),x(i,1):x(i,2)) = chess;
                obj.chess(i).background = colors(i);
            end

            %% Fix positions
            if(strcmp(obj.name,'red'))
                disp('fuck');
            end
            obj = fix_positions(x, obj, BW, false);
            %% Se su una stessa riga ho un centroide in più rimuovilo
            for i=1:size(obj.chess,2)
                X = obj.chess(i).center_x;
                Y = obj.chess(i).center_y;
                obj.chess(i).v_lines_centroid = cell(1,size(X,2));
                for k=1:size(X,2)
                    if(isempty(find(X(:,k)==0, 1)))
                        [fitresult, ~] = createLineInv(Y(:,k),X(:,k),size(obj.color_mask));
                        if(abs(fitresult.p1) == Inf || abs(fitresult.p2) == Inf)
                           [fitresult, ~] = createLine(X(:,k),Y(:,k));
                        end
                        obj.chess(i).v_lines_centroid{k} = fitresult;
                    else
                        obj.chess(i).v_lines_centroid{k} = [];
                    end
                end   
            end       
            maskedRGBImage(repmat(~obj.color_mask ,[1 1 3]))= 0;
            obj.masked_rgb = im2uint8(maskedRGBImage);
            obj.isEmpty = false;
            close all;
            obj_chess(kk)=obj;
        end
    end
end