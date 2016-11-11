function verticalEdgeImage = getVImage( bw, mask, maskH )
%     mask = imdilate(mask,strel('rectangle',[2 12]));
%     mask = imerode(mask,strel('rectangle',[14 2]));
    bw_edge = filledgegaps(bw.*mask,1);
    CC = bwconncomp(~bw_edge,4);
    for i=1:size(CC.PixelIdxList,2)
        if(size(CC.PixelIdxList{i},1)<100)
            bw_edge(CC.PixelIdxList{i}) = 1;
        end
    end
    bw_edge = filledgegaps(bw_edge,1);
    horizontal_image = imdilate(maskH,strel('disk',10));
    verticalEdgeImage = imfilter(bw_edge, [-1 0 1]);
    verticalEdgeImage = bwareaopen(verticalEdgeImage, 20); % Parametro
    verticalEdgeImage = (verticalEdgeImage - horizontal_image)>0;
end

