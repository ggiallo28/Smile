function Container = split_refletions(Container, fuse, checker_vector)
    obj_chess = Container.obj_chess;
    
    enlarged = imdilate(fuse,strel('disk',20));
    CC = bwconncomp(rgb2gray(enlarged)~=0);
    assert(CC.NumObjects == 3, 'Ci sono troppi riflessi');
    maskCenter = false(size(rgb2gray(enlarged))); maskCenter(CC.PixelIdxList{2}) = 1;
    maskLeft = false(size(rgb2gray(enlarged))); maskLeft(CC.PixelIdxList{1}) = 1;
    maskRight = false(size(rgb2gray(enlarged))); maskRight(CC.PixelIdxList{3}) = 1;
    positions = [{'Left'} {'Center'} {'Right'}];
    types = [{'Primary'} {'Secondary'} {'Real'}];
    order = []; k =1;
    % Order Vector of chess
    for l=1:size(obj_chess,1)
        if ( ~obj_chess(l).isEmpty )
            for i = 1:size(obj_chess(l).chess,2)
                order = [order,[l;i;obj_chess(l).chess(i).centroid(1);obj_chess(l).chess(i).centroid(2)]]; k = k+1;
            end
        end
    end
    % Find Objects Position
    for l=1:size(obj_chess,1)
        left = 0;
        right = 0;
        arr = [];
        if ( ~obj_chess(l).isEmpty )
            for i = 1:size(obj_chess(l).chess,2)
                v(1) = sum(sum(obj_chess(l).chess(i).mask&maskLeft));
                v(2) = sum(sum(obj_chess(l).chess(i).mask&maskCenter));
                v(3) = sum(sum(obj_chess(l).chess(i).mask&maskRight));
                idx = find(v==max(v));
                obj_chess(l).chess(i).position = positions(idx);
                if idx==2
                    obj_chess(l).chess(i).type =types(3);
                end
            end
        end
    end
    % Adds label
    label = cell(4,k-1); k =1;
    for l=1:size(obj_chess,1)
        if ( ~obj_chess(l).isEmpty )
            for i = 1:size(obj_chess(l).chess,2)
                label(1,k) = cellstr(obj_chess(l).name);
                label(2,k) = cellstr(obj_chess(l).chess(i).background);
                label(3,k) = cellstr(obj_chess(l).chess(i).position);
                k = k+1;
            end
        end
    end
    % Ordering of chess
    minimum = order(3,1); % Inizializzazione
    for l=1:size(order,2)-1
        for i=l:size(order,2)
            if(order(3,l)>order(3,i))
                tmp = order(:,l);
                tmp2 = label(:,l);
                order(:,l) = order(:,i);
                label(:,l) = label(:,i);
                order(:,i) = tmp;
                label(:,i) = tmp2;
            end
        end
    end
    % Find Reflection type using colors
    for l=1:size(obj_chess,1)
        left = 0;
        right = 0;
        arr = [];
        if ( ~obj_chess(l).isEmpty )
            for i = 1:size(obj_chess(l).chess,2)
                v(1) = sum(sum(obj_chess(l).chess(i).mask&maskLeft));
                v(2) = sum(sum(obj_chess(l).chess(i).mask&maskCenter));
                v(3) = sum(sum(obj_chess(l).chess(i).mask&maskRight));
                idx = find(v==max(v));
                if(idx==1)
                   left = left +1; 
                end
                if(idx==3)
                   right = right +1; 
                end
                arr = [arr,obj_chess(l).chess(i).position];
            end
            %% TODO Sistemare sta cosa non va bene fare assunzioni sul numero di riflessi che vedi
            if(left >=2 && right>=2)
                idx = find(strcmp(arr,positions(1)));
                obj_chess(l).chess(idx(1)).type = types(2);
                obj_chess(l).chess(idx(2)).type = types(1); 
                idx = find(strcmp(arr,positions(3)));
                obj_chess(l).chess(idx(1)).type = types(1);
                obj_chess(l).chess(idx(2)).type = types(2);
            end
            if(left >= 2 && right == 1)
               idx = find(strcmp(arr,positions(1)));
               obj_chess(l).chess(idx(1)).type = types(2);
               obj_chess(l).chess(idx(2)).type = types(1); 
               idx = find(strcmp(arr,positions(3)));
               obj_chess(l).chess(idx(1)).type = types(1);
            end
            if(right >= 2 && left == 1)
               idx = find(strcmp(arr,positions(3)));
               obj_chess(l).chess(idx(1)).type = types(1);
               obj_chess(l).chess(idx(2)).type = types(2);  
               idx = find(strcmp(arr,positions(1)));
               obj_chess(l).chess(idx(1)).type = types(1);
            end
            empty_type = false;
            for i = 1:size(obj_chess(l).chess,2)
                if(isempty(obj_chess(l).chess(i).type))
                    empty_type = true;
                end
            end
            % Se sono rimaste tessere senza nessuna labels associata
            if(empty_type)
                for i = 1:size(obj_chess(l).chess,2)
                    % Se non si tratta di una tessera centrale, non ha senso analzizarla di sicuro non è ne primaria ne secondaria
                    if(~strcmp(obj_chess(l).chess(i).position,positions(2)))
                        curr_color = name2code(obj_chess(l).name);
                        % Cerchiamo la tessera corrente nel vettore di ordine
                        idx = find(order(3,:) == obj_chess(l).chess(i).centroid(1));
                        % Cerchiamo la posizione della tessera centrale
                        labs = find(strcmp(label(3,:),positions(2)));
                        % Se idx è precedente della prima tessera centrale:
                        %   > è la prima tessera? Allora guardiamo a destra
                        %   > non è la prima tessera? Allora guardiamo a sinistra
                        % Se idx è successivo all'ultima tessera centrale:
                        %   > è la prima l'ultima tessera? Allora guardiamo a sinistra
                        %   > non è l'ultima tessera? Allora guardiamo a destra
                        % Usiamo l'ordine dei colori per determinare l'ordinamento, i casi con due colori a sinistra 
                        % o a destra sono banali e sono stati gestiti precedentemente.
                        if(idx < labs(1))
                            if(idx~=1)
                                color = name2code(obj_chess(order(1,idx-1)).name);                   
                                % Sfruttiamo il checker vector, se i colori sono invertiti rispetto all'ordine originale allora
                                % abbiamo un riflesso primario, altrimenti abbiamo un riflesso secondario. 
                                check = checker_vector(2,:,:);
                                prev = find(check(:,:,1) == color(1) & check(:,:,2) == color(2) & check(:,:,3) == color(3));
                                curr = find(check(:,:,1) == curr_color(1) & check(:,:,2) == curr_color(2) & check(:,:,3) == curr_color(3));
                                if(prev==1 && curr == size(check,2))
                                    obj_chess(l).chess(i).type = types(1);
                                elseif (prev==size(check,2) && curr ==1)
                                    obj_chess(l).chess(i).type = types(2);
                                elseif( prev < curr )
                                    obj_chess(l).chess(i).type = types(2);
                                else
                                    obj_chess(l).chess(i).type = types(1);
                                end
                            else
                                color = name2code(obj_chess(order(1,idx+1)).name);
                                check = checker_vector(2,:,:);
                                next = find(check(:,:,1) == color(1) & check(:,:,2) == color(2) & check(:,:,3) == color(3));
                                curr = find(check(:,:,1) == curr_color(1) & check(:,:,2) == curr_color(2) & check(:,:,3) == curr_color(3));                   
                                if(next == size(check,2) && curr == 1)
                                    obj_chess(l).chess(i).type = types(1);
                                elseif (next==1 && curr == size(check,2))
                                    obj_chess(l).chess(i).type = types(2);
                                elseif(curr < next)
                                     obj_chess(l).chess(i).type = types(2);
                                else
                                     obj_chess(l).chess(i).type = types(1);
                                end
                            end
                        elseif(idx > labs(end))
                            if (idx~=size(order,2))
                                color = name2code(obj_chess(order(1,idx+1)).name);
                                check = checker_vector(2,:,:);
                                next = find(check(:,:,1) == color(1) & check(:,:,2) == color(2) & check(:,:,3) == color(3));
                                curr = find(check(:,:,1) == curr_color(1) & check(:,:,2) == curr_color(2) & check(:,:,3) == curr_color(3));
                                if(next == size(check,2) && curr == 1)
                                    obj_chess(l).chess(i).type = types(1);
                                elseif (next==1 && curr == size(check,2))
                                    obj_chess(l).chess(i).type = types(2);
                                elseif(curr < next)
                                     obj_chess(l).chess(i).type = types(2);
                                else
                                     obj_chess(l).chess(i).type = types(1);
                                end
                            else
                                color = name2code(obj_chess(order(1,idx-1)).name);
                                check = checker_vector(2,:,:);
                                prev = find(check(:,:,1) == color(1) & check(:,:,2) == color(2) & check(:,:,3) == color(3));
                                curr = find(check(:,:,1) == curr_color(1) & check(:,:,2) == curr_color(2) & check(:,:,3) == curr_color(3));
                                if(prev==1 && curr == size(check,2))
                                    obj_chess(l).chess(i).type = types(1);
                                elseif (prev==size(check,2) && curr ==1)
                                    obj_chess(l).chess(i).type = types(2);
                                elseif(prev < curr)
                                     obj_chess(l).chess(i).type = types(2);
                                else
                                     obj_chess(l).chess(i).type = types(1);
                                end 
                            end
                        end
                    end 
                end
            end
        end
    end
    imshow(maskCenter|maskLeft|maskRight)
    figure, imshow(enlarged);
    for j=1:size(order,2)
        idx_chess_vector = order(1,j);
        idx_color_chess = order(2,j);
        label(4,j) = obj_chess(idx_chess_vector).chess(idx_color_chess).type;
    end
    Container.order = order;
    Container.label = label;
    Container.positions = positions;
    Container.types = types; 
    Container.obj_chess = obj_chess;
end