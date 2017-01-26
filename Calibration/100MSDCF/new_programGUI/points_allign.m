function Container = points_allign(Container)
    % Per ogni combinazione colore/sfondo controlliamo quant'e il massimo
    % numero di colonne, aggiungiamo zeri alle matrici che hanno un numero
    % inferiore di colonne
    %% INIT
    obj_chess = Container.obj_chess;
    label = Container.label;
    order =  Container.order;
    types = Container.types;
    for l=1:size(obj_chess,1)
        if ( ~obj_chess(l).isEmpty )
            max_v = 5; % numero di punti per checkerboard     
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
                if(hasLeft(label,idx_order, Container))
                  obj_chess(l).chess(q).intersections_x = ...
                      allignleftdouble(obj_chess(l).chess(q).intersections_x);
                  obj_chess(l).chess(q).intersections_y = ...
                      allignleftdouble(obj_chess(l).chess(q).intersections_y);
                end
                if(hasRight(label,idx_order, Container))
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
                   obj_chess(l).chess(q).isFlipped = true;
                end
                obj_chess(l).chess(q).h_matrix = zeros(size(obj_chess(l).chess(q).intersections_y,1),max_v);
                obj_chess(l).chess(q).angle_matrix = zeros(size(obj_chess(l).chess(q).intersections_y,1),max_v);
            end
        end
    end
    Container.obj_chess = obj_chess;
end