function [fitresult_top, fitresult_bot] = getTopBotLines(color_square)
    tmp = bwconvhull(color_square);
    color_square_bb = regionprops(tmp,'BoundingBox');
    color_square_bb = color_square_bb.BoundingBox;
    color_square_edge_v = abs(imfilter(imfill(color_square,'holes'),[-1 0 1])) + abs(imfilter(imfill(color_square,'holes'),[1 0 -1]));
    color_square_edge_h = abs(imfilter(imfill(color_square,'holes'),[1 0 -1]'))+abs(imfilter(imfill(color_square,'holes'),[-1 0 1]'));
%% HORIZONTAL
    color_square_edge_h = color_square_edge_h-color_square_edge_v;
    color_square_edge_h(color_square_edge_h<0) = 0;
    color_square_edge_h = imclose(color_square_edge_h,strel('disk',2));
    area = bwconncomp(color_square_edge_h); acc = 0;
    for count=1:size(area.PixelIdxList,2)
        acc = acc + size(cell2mat(area.PixelIdxList(count)),1);
    end
    filt = round((acc/size(area.PixelIdxList,2))/4);    
    color_square_edge_h = bwareaopen(color_square_edge_h, filt); % Parametro
    CC_edge_h = regionprops(color_square_edge_h,'Centroid','Image');
    CC_h = bwconncomp(color_square_edge_h);

    yT = color_square_bb(2);
    yB = color_square_bb(2)+color_square_bb(4);
    vect_top = zeros(2,size(CC_edge_h,1));
    vect_bot = zeros(2,size(CC_edge_h,1));
    size_chess = abs(yT-yB);
    size_square = size_chess./5;
    
    for count=1:size(CC_edge_h,1)   
        yC = CC_edge_h(count).Centroid(2);
        vect_top(1,count) = count;
        vect_top(2,count) = abs(yC-yT);
        vect_bot(1,count) = count;
        vect_bot(2,count) = abs(yC-yB);
    end
    
    [vect_top(2,:),vect_top(1,:)] = sort(vect_top(2,:));
    [vect_bot(2,:),vect_bot(1,:)] = sort(vect_bot(2,:));
    
    iitop = find(vect_top(2,:)<0.7*size_square);
    iibot = find(vect_bot(2,:)<0.7*size_square);
    tmptop = CC_h.PixelIdxList(vect_top(1,iitop));
    tmpbot = CC_h.PixelIdxList(vect_bot(1,iibot));
    cctop = []; ccbot=[];   
    for count=1:size(tmptop,2)
        cctop = [cctop; cell2mat(tmptop(count))];
    end
    for count=1:size(tmpbot,2)
        ccbot = [ccbot; cell2mat(tmpbot(count))];
    end
    [idx_top,idy_top]=ind2sub(size(color_square),cctop);
    [fitresult_top, ~] = createLine(idy_top,idx_top);
    [fitresult_top_tmp, isDone] = adjustLine(fitresult_top,idy_top,idx_top);
    if(isDone)
        fitresult_top = fitresult_top_tmp;
    end
    %imshow(color_square); hold on;
%    plot(fitresult_top); l1 = legend(); set(l1,'visible','off');  
    [idx_bot,idy_bot]=ind2sub(size(color_square),ccbot);
    [fitresult_bot, ~] = createLine(idy_bot,idx_bot);
    [fitresult_bot_tmp, isDone] = adjustLine(fitresult_bot,idy_bot,idx_bot);
    if(isDone)
        fitresult_bot = fitresult_bot_tmp;
    end
    %plot(fitresult_bot); l1 = legend(); set(l1,'visible','off');
end