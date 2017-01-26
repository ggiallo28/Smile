function bool_value = hasRight(label,id, Container)
    if(id==size(label,2))
        bool_value = false;
        return;
    end
    order = Container.order;
    order_left = order(:,id+1); %[id_in_obj_chess; id_chess;...]
    chess_int_matrix = Container.obj_chess(order_left(1)).chess(order_left(2)).intersections_x;
    label_left = label(:,id+1);
    label_curr = label(:,id);
    if( strcmp(label_left(3),label_curr(3)) && strcmp(label_left(4),label_curr(4)) )
        if (isempty(chess_int_matrix) || sum(sum(chess_int_matrix))==0)
            bool_value = false;
        else
            bool_value = true;
        end
    else
       bool_value = false;
    end
end
