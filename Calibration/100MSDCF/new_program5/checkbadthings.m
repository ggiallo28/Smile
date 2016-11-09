function [color_square, y1, y2, isDone] = checkbadthings(color_square,mask,Threshold, op_th)
    mask_bbox = regionprops(mask,'BoundingBox');
    bw_cut = medfilt2(sum(color_square'),[1,2]);
    idx = find(bw_cut~=0);
    y1 = min(idx);
    y2 = max(idx);
    h_color_mask = y2-y1;
    isDone = false;
    if(0.9*mask_bbox.BoundingBox(4)>h_color_mask)
        color_square = imopen(color_square,strel('rectangle',[2,op_th]));
        color_square = bwareaopen(color_square,Threshold);
        bw_cut = sum(color_square');
        idx = find(bw_cut~=0);
        y1 = min(idx);
        y2 = max(idx);
        isDone = true;
    end
end