function [image1, image2] = getHImage(pos,type, obj_chess, order, label, dim, maskI)
    idx_pos = find(strcmp(label(3,:),pos));
    idx_type = find(strcmp(label(4,:),type));
    idx = intersect(idx_type,idx_pos);
    obj_chess_idx = order(1,idx);
    chess_idx = order(2,idx);
    image1 = false(dim);
    image2 = false(dim);
    for i=1:size(idx,2)
        tmp_image = false(dim);
        for j=1:size(obj_chess(obj_chess_idx(i)).chess(chess_idx(i)).h_lines,2)
            tmp_image = tmp_image | line2image(obj_chess(obj_chess_idx(i)).chess(chess_idx(i)).h_lines{j},dim);
        end
        image1 = image1 | tmp_image & imdilate(obj_chess(obj_chess_idx(i)).chess(chess_idx(i)).ch_mask,strel('rectangle',[20,2]));
        image2 = image2 | image1 | tmp_image & ~maskI;
    end   
end

