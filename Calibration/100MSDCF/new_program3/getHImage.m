function image = getHImage(pos,type, obj_chess, order, label, dim)
    idx_pos = find(strcmp(label(3,:),pos));
    idx_type = find(strcmp(label(4,:),type));
    idx = intersect(idx_type,idx_pos);
    obj_chess_idx = order(1,idx);
    chess_idx = order(2,idx);
    image = false(dim);
    for i=1:size(idx,2)
        tmp_image = false(dim);
        for j=1:size(obj_chess(obj_chess_idx(i)).chess(chess_idx(i)).h_lines,2)
            tmp_image = tmp_image | line2image(obj_chess(obj_chess_idx(i)).chess(chess_idx(i)).h_lines{j},dim);
        end
        image = image | tmp_image & obj_chess(obj_chess_idx(i)).chess(chess_idx(i)).ch_mask;
    end   
end

