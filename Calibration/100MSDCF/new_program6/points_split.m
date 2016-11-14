function Container = points_split(Container, pointsArray)
    % Scansioniamo le tessere, facciamo intersezione con le maschere di ogni
    % singola tessera, i punti che ricadono dell'intersezione vengono assegnati
    % alla tessera associata alla specifica maschera.
    % Per fare questo usiamo il vettore delle label e order creato
    % precedentemente per poter facilemente recuperare gli indici di riga e
    % colonna associati ad ogni singolo elemento.
%% INIT
    I = Container.I;
    obj_chess = Container.obj_chess;
    isC = Container.isC;
    isL1 = Container.isL1;
    isL2 = Container.isL2;
    isR1 = Container.isR1;
    isR2 = Container.isR2;
    label = Container.label;
    positions =  Container.positions;
    order =  Container.order;
    types = Container.types;
% LOGIC
    for k=1:5
        switch(k)
            case 1
                if ( ~isC )
                    continue;
                end
                idx = find(strcmp(label(3,:),positions(2)));
            case 2
                if ( ~isL1 )
                    continue;
                end
                idx = find(strcmp(label(3,:),positions(1)) & strcmp(label(4,:),types(1)));
            case 3
                if ( ~isL2 )
                    continue;
                end
                idx = find(strcmp(label(3,:),positions(1)) & strcmp(label(4,:),types(2)));
            case 4
                if ( ~isR1 )
                    continue;
                end
                idx = find(strcmp(label(3,:),positions(3)) & strcmp(label(4,:),types(1)));
            case 5
                if ( ~isR2 )
                    continue;
                end
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
    %         for ii=1:size(x_cj_matrix_copy,1)
    %             for jj=1:size(x_cj_matrix_copy,2)
    %                 imshow(I); hold on; scatter(x_cj_matrix_copy(ii,jj),y_rj_matrix_copy(ii,jj));
    %                 pause
    %             end
    %         end
            imshow([I;repmat(255.*ch_mask,1,1,3)]); hold on;
            scatter(x_cj_matrix_copy(:),y_rj_matrix_copy(:));
            pause(1)
        end    
    end
    Container.obj_chess = obj_chess;
end