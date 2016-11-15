function Container = generate_mapping(Container, checker_vector)
    obj_chess = Container.obj_chess;
    angle_sector = 60;
    angle_square = 15;
    height_square = 2.6;
    axis_line = reshape([192,192,192;192,192,192],2,1,3);
    checher_vector_with_axis = uint8([checker_vector(:,1:3,:),axis_line,checker_vector(:,4:6,:)]); 
    % La checker è suddivisa verticalmente in 4 quadrati, ciò significa 5 punti compresi quelli di intersezione
    % Ogni quadrante rappresenta 60°, quindi ad ogni punto corrisponde un angolo di 15°
    for l=1:size(obj_chess,1)
        if ( ~obj_chess(l).isEmpty )
            disp(['Calculating ',obj_chess(l).name,'...']);
            axis_distance = axisdistance(obj_chess(l).name,checher_vector_with_axis);
            adjusted_distance = axis_distance-1*sign(axis_distance); % Se positivo togliamo 1, se negativo aggiungiamo 1, le tessere gialle e rosse stanno a distanza zero.        
            for i=1:size(obj_chess(l).chess(1).intersections_x,1)
                for j=1:size(obj_chess(l).chess(1).intersections_x,2)
                    if(sign(axis_distance) > 0)
                        adjusted_offset = (j-1)*angle_square;
                    else
                        adjusted_offset = (j-size(obj_chess(l).chess(1).intersections_x,2))*angle_square;
                    end
                    for k=1:size(obj_chess(l).chess,2)
                        obj_chess(l).chess(k).angle_matrix(i,j) =  adjusted_distance*angle_sector + adjusted_offset;
                        obj_chess(l).chess(k).h_matrix(i,j) = i*height_square;
                    end
                end
            end
        end
    end
    Container.obj_chess = obj_chess;
end