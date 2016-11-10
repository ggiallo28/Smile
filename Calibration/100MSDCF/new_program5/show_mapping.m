function show_mapping(obj_chess, I)
    for l=1:size(obj_chess,1)
        if ( ~obj_chess(l).isEmpty )
            for i=1:size(obj_chess(l).chess(1).intersections_x,1)
                for j=1:size(obj_chess(l).chess(1).intersections_x,2)
                    imshow(I); hold on;
                    vect_x = []; 
                    vect_y = []; 
                    for k =1:size(obj_chess(l).chess,2)
                        vect_x = [vect_x, obj_chess(l).chess(k).intersections_x(i,j)];
                        vect_y = [vect_y, obj_chess(l).chess(k).intersections_y(i,j)];
                    end
                    obj_chess(l).chess(1).h_matrix(i,j)
                    obj_chess(l).chess(1).angle_matrix(i,j)
                    scatter(vect_x(:),vect_y(:)); hold off
                    pause
                end
            end
        end
    end
end