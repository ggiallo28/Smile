function color_square = fitSquare( color_square )
    color_square_bb = regionprops(color_square,'BoundingBox');
    color_square_bb = color_square_bb.BoundingBox; imshow(color_square); hold on
    plot([color_square_bb(1) color_square_bb(1)+color_square_bb(3)  color_square_bb(1)+color_square_bb(3)...
    color_square_bb(1) color_square_bb(1)],[color_square_bb(2)  color_square_bb(2) color_square_bb(2)+color_square_bb(4)...
    color_square_bb(2)+color_square_bb(4) color_square_bb(2)]);
    color_square_edge_v = abs(imfilter(imfill(color_square,'holes'),[-1 0 1])) + abs(imfilter(imfill(color_square,'holes'),[1 0 -1]));
    color_square_edge_h = abs(imfilter(imfill(color_square,'holes'),[1 0 -1]'))+abs(imfilter(imfill(color_square,'holes'),[-1 0 1]'));

%% VERTICAL
    color_square_edge_v = color_square_edge_v - color_square_edge_h;
    color_square_edge_v(color_square_edge_v<0) = 0;
    color_square_edge_v = imclose(color_square_edge_v,strel('disk',2));
    color_square_edge_v = bwareaopen(color_square_edge_v, 100); % Parametro
    CC_edge_v = regionprops(color_square_edge_v,'Centroid','Image');
    CC_v = bwconncomp(color_square_edge_v);

    xL = color_square_bb(1);
    xR = color_square_bb(1)+color_square_bb(3);
    vect_left = zeros(4,size(CC_edge_v,1));
    vect_right = zeros(4,size(CC_edge_v,1));
    
    for count=1:size(CC_edge_v,1)   
        xC = CC_edge_v(count).Centroid(1);
        vect_left(1,count) = count;
        vect_left(2,count) = abs(xC-xL);
        vect_left(3,count) = CC_edge_v(count).Centroid(2);
        vect_right(1,count) = count;
        vect_right(2,count) = abs(xC-xR);   
        vect_right(3,count) = CC_edge_v(count).Centroid(2);
    end
    [vect_left(2,:),vect_left(1,:)] = sort(vect_left(2,:));
    [vect_right(2,:),vect_right(1,:)] = sort(vect_right(2,:));
    vect_right(3,:) = vect_right(3,vect_right(1,:));   
    vect_right(4,:) = abs(diff([vect_right(3,1),vect_right(3,:)]));
    vect_left(4,:) = abs(diff([vect_left(3,1),vect_left(3,:)]));
    
    ccleft = [cell2mat(CC_v.PixelIdxList(vect_left(1,1))); cell2mat(CC_v.PixelIdxList(vect_left(1,2)))];
    ccright = [cell2mat(CC_v.PixelIdxList(vect_right(1,1))); cell2mat(CC_v.PixelIdxList(vect_right(1,2)))];
    yStd = std([vect_left(4,2), vect_right(4,2)]);
    yDist = max([vect_left(4,1), vect_left(4,2), vect_right(4,1), vect_right(4,2)])+2*yStd;
    if(vect_left(4,3) < yDist)
        ccleft = [ccleft;cell2mat(CC_v.PixelIdxList(vect_left(1,3)))];
    end
    if(vect_right(4,3) < yDist)
        ccright = [ccright;cell2mat(CC_v.PixelIdxList(vect_right(1,3)))];
    end    
    [idx_left,idy_left]=ind2sub(size(color_square),ccleft);
    [fitresult_left, ~] = createLine(idy_left,idx_left);
    plot(fitresult_left);
    [idx_right,idy_right]=ind2sub(size(color_square),ccright);
    [fitresult_right, ~] = createLine(idy_right,idx_right);
    plot(fitresult_right);
%% HORIZONTAL
    color_square_edge_h = color_square_edge_h-color_square_edge_v;
    color_square_edge_h(color_square_edge_h<0) = 0;
    color_square_edge_h = imclose(color_square_edge_h,strel('disk',2));
    color_square_edge_h = bwareaopen(color_square_edge_h, 100); % Parametro
    CC_edge_h = regionprops(color_square_edge_h,'Centroid','Image');
    CC_h = bwconncomp(color_square_edge_h);

    yT = color_square_bb(2);
    yB = color_square_bb(2)+color_square_bb(4);
    vect_top = zeros(2,size(CC_edge_h,1));
    vect_bot = zeros(2,size(CC_edge_h,1));
    
    for count=1:size(CC_edge_h,1)   
        yC = CC_edge_h(count).Centroid(2);
        vect_top(1,count) = count;
        vect_top(2,count) = abs(yC-yT);
        vect_bot(1,count) = count;
        vect_bot(2,count) = abs(yC-yB);
    end
    
    [vect_top(2,:),vect_top(1,:)] = sort(vect_top(2,:));
    [vect_bot(2,:),vect_bot(1,:)] = sort(vect_bot(2,:));
    size_chess = abs(yT-yB);
    size_square = size_chess./5;
    iitop = find(vect_top(2,:)<0.85*size_square);
    iibot = find(vect_bot(2,:)<0.85*size_square);
    tmptop = CC_h.PixelIdxList(vect_top(1,iitop));
    tmpbot = CC_h.PixelIdxList(vect_bot(1,iibot));
    cctop = []; ccbot=[];   
    for count=1:size(tmptop,2)
        cctop = [cctop; cell2mat(tmptop(count))];
    end
    for count=1:size(tmptop,2)
        ccbot = [ccbot; cell2mat(tmpbot(count))];
    end
    [idx_top,idy_top]=ind2sub(size(color_square),cctop);
    [fitresult_top, ~] = createLine(idy_top,idx_top);
    plot(fitresult_top);
    [idx_bot,idy_bot]=ind2sub(size(color_square),ccbot);
    [fitresult_bot, ~] = createLine(idy_bot,idx_bot);
    plot(fitresult_bot);
    
    xx = 1:0.001:size(color_square,2);
    [~,ii] = min(abs(fitresult_bot(xx(:)) - fitresult_left(xx(:))));
    ybotleft = fitresult_bot(xx(ii));
    xbotleft = xx(ii);
    scatter(xbotleft,ybotleft)
    
    [~,ii] = min(abs(fitresult_bot(xx(:)) - fitresult_right(xx(:))));
    ybotright = fitresult_bot(xx(ii));
    xbotright = xx(ii);
    scatter(xbotright,ybotright)
    
    [~,ii] = min(abs(fitresult_top(xx(:)) - fitresult_left(xx(:))));
    ytopleft = fitresult_bot(xx(ii));
    xtopleft = xx(ii);
    scatter(xtopleft,ytopleft)
    
    [~,ii] = min(abs(fitresult_top(xx(:)) - fitresult_right(xx(:))));
    ytopright = fitresult_bot(xx(ii));
    xtopright = xx(ii);
    scatter(ytopright,xtopright)
    
    

        
    idx = find(color_square == 1);
    [idx,idy]=ind2sub(size(color_square),idx);
    j = boundary(idx,idy,0.1); % Parametro
    color_square = poly2mask(idy(j),idx(j), size(color_square,1), size(color_square,2));  
end
