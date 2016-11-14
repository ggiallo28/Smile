function [BW, y] = twosquares_adjust(x, y, BW, RGB)
% Funzione che risolvere gestisce la situazione in cui abbiamo solo due componenti connesse all'interno di una sezione, cioè quando abbiamo solo i due quadrati intermedi.
% per calcolare la maschera complementare dobbiamo simulare quello di sotpra e quello di sotto altrimenti la regione sarebbe ritagliata troppo piccola.
    for i=1:size(x,1)
        mask = false(size(BW,1),size(BW,2));
        mask(y(i,1):y(i,2),x(i,1):x(i,2))=true;
        color_square =  BW&mask;
        CC = bwconncomp(color_square); 
        if(size(CC.PixelIdxList,2) == 2) 
            stats = cell2mat(struct2cell(regionprops(color_square,'BoundingBox'))); % Devo essere sicuro siano due blob isolati
            if((stats(4) + stats(8))< 0.5*size(BW,1))
                mask = false(size(BW,1),size(BW,2));
                en = 0.9*min(stats(4),stats(8)); 
                if(round(y(i,1)-en)>0)
                    y(i,1) = round(y(i,1)-en);
                else
                    y(i,1) = 2;
                end
                if(round(y(i,2)+en)<size(BW,1))
                    y(i,2) = round(y(i,2)+en);
                else
                    y(i,2) = size(BW,1)-2;
                end
                mask(y(i,1):y(i,2),x(i,1):x(i,2))=true;
                props = regionprops(color_square,'Centroid','MinorAxisLength');
                len = max([props(1).MinorAxisLength,props(2).MinorAxisLength]);
                [fitresult, ~] = createLine([props(1).Centroid(2),props(2).Centroid(2)],[props(1).Centroid(1),props(2).Centroid(1)]);
                Xtop = fitresult(y(i,1));
                Xbot = fitresult(y(i,2));
                color_square =  BW&mask;  
                % Se abbiamo due square il rettangolo sopra e sotto manca, lo simuliamo replicando i contorni inferiore e superiore
                color_square_edge = imfilter(color_square,[-1 0 1]') | imfilter(color_square,[1 0 -1]');
                color_square_edge = bwareaopen(color_square_edge,round(0.9*len));
                CC_color_square_edge = bwconncomp(color_square_edge,8);
                centroid_color_square_edge = reshape(cell2mat(struct2cell(regionprops(color_square_edge,'Centroid'))),2,CC_color_square_edge.NumObjects);
                edge_to_top = centroid_color_square_edge(2,:)-y(i,1);
                edge_to_bot = y(i,2)-centroid_color_square_edge(2,:);
                [~, idx_top] = sort(edge_to_top);
                [~, idx_bot] = sort(edge_to_bot); 
                top_edge = CC_color_square_edge.PixelIdxList{idx_top(1)};
                [y_top,x_top] = ind2sub(size(color_square),top_edge);
                x_top = x_top + floor(Xtop - mean(x_top));
                y_top = y_top - min(y_top)+y(i,1);
                top_edge = sub2ind(size(color_square),y_top,x_top);
                color_square(top_edge) = 1;          

                bot_edge = CC_color_square_edge.PixelIdxList{idx_bot(1)};
                [y_bot,x_bot] = ind2sub(size(color_square),bot_edge);            
                x_bot = x_bot + floor(Xbot - mean(x_bot));
                y_bot = y_bot - max(y_bot) + y(i,2);
                top_edge = sub2ind(size(color_square),y_bot,x_bot);
                color_square(top_edge) = 1;
                BW(y(i,1):y(i,2),x(i,1):x(i,2)) = color_square(y(i,1):y(i,2),x(i,1):x(i,2));
            end
        end
    end
end