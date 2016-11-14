function [color_square, chess]  = fitSquare( color_square, x1, x2, chess)
%     idx = find(color_square == 1);
%     [idx,idy]=ind2sub(size(color_square),idx);
%     j = boundary(idx,idy,0.1); % Parametro
%     tmp = poly2mask(idy(j),idx(j), size(color_square,1), size(color_square,2));  
% Se il contorno di sopra o sotto � troppo piccolo usare bbox
% Per eliminare merde laterali puoi calcolare la maskera con la fit square,
% erodere e risegmentare
% dilatare un pochino il bbox perch� nella fit square pu� succedere che si
% vada a finire fuori
% proiettare a sinistra e a destra tutti i contorni prima di interpolare

    tmp = bwconvhull(color_square);
    color_square_bb = regionprops(tmp,'BoundingBox');
    color_square_bb = color_square_bb.BoundingBox; imshow(color_square); hold on
    plot([color_square_bb(1) color_square_bb(1)+color_square_bb(3)  color_square_bb(1)+color_square_bb(3)...
    color_square_bb(1) color_square_bb(1)],[color_square_bb(2)  color_square_bb(2) color_square_bb(2)+color_square_bb(4)...
    color_square_bb(2)+color_square_bb(4) color_square_bb(2)]); l1 = legend(); set(l1,'visible','off');
    color_square_edge_v = abs(imfilter(imfill(color_square,'holes'),[-1 0 1])) + abs(imfilter(imfill(color_square,'holes'),[1 0 -1]));
    color_square_edge_h = abs(imfilter(imfill(color_square,'holes'),[1 0 -1]'))+abs(imfilter(imfill(color_square,'holes'),[-1 0 1]'));
    [LTp, RTp, LBp, RBp]  = getRotatedRectangle(color_square);
    [rct_line_left, ~] = createLine([LBp(2) LTp(2)],[LBp(1) LTp(1)]); % Al contrario perch� ho necessit� di sostituire la y
    [rct_line_right, ~] = createLine([RBp(2) RTp(2)],[RBp(1) RTp(1)]);
    [rct_line_left_test, ~] = createLine([LBp(1) LTp(1)],[LBp(2) LTp(2)]);
    [rct_line_right_test, ~] = createLine([RBp(1) RTp(1)],[RBp(2) RTp(2)]);
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
    plot(fitresult_top); l1 = legend(); set(l1,'visible','off');  
    [idx_bot,idy_bot]=ind2sub(size(color_square),ccbot);
    [fitresult_bot, ~] = createLine(idy_bot,idx_bot);
    [fitresult_bot_tmp, isDone] = adjustLine(fitresult_bot,idy_bot,idx_bot);
    if(isDone)
        fitresult_bot = fitresult_bot_tmp;
    end
    plot(fitresult_bot); l1 = legend(); set(l1,'visible','off');
%% HORIZONTAL 2
    yLB = fitresult_bot(x1);
    yRB = fitresult_bot(x2);
    yLT = fitresult_top(x1);
    yRT = fitresult_top(x2);
    yLL = linspace(yLB,yLT,6);
    yRR = linspace(yRB,yRT,6);
    [line1, ~] = createLine([x1 x2],[yLL(2) yRR(2)]);
    plot(line1); l1 = legend(); set(l1,'visible','off');  
    [line2, ~] = createLine([x1 x2],[yLL(3) yRR(3)]);
    plot(line2); l1 = legend(); set(l1,'visible','off');  
    [line3, ~] = createLine([x1 x2],[yLL(4) yRR(4)]);
    plot(line3); l1 = legend(); set(l1,'visible','off');  
    [line4, ~] = createLine([x1 x2], [yLL(5) yRR(5)]);
    plot(line4); l1 = legend(); set(l1,'visible','off');  
    chess.h_lines = cell(1,6);
    chess.h_lines{1} = fitresult_bot;
    chess.h_lines{2} = line1;
    chess.h_lines{3} = line2;
    chess.h_lines{4} = line3;
    chess.h_lines{5} = line4;
    chess.h_lines{6} = fitresult_top;
    
    image = line2image(line1,size(color_square)) | line2image(line2,size(color_square)) |...
        line2image(line3,size(color_square)) | line2image(line4,size(color_square));
%% VERTICAL
    image = imdilate(image,strel('disk',10));
    color_square_edge_v = color_square_edge_v - color_square_edge_h-image;
    color_square_edge_v(color_square_edge_v<0) = 0;
    color_square_edge_v = imclose(color_square_edge_v,strel('disk',2));
    area = bwconncomp(color_square_edge_v); acc = 0;
    for count=1:size(area.PixelIdxList,2)
        acc = acc + size(cell2mat(area.PixelIdxList(count)),1);
    end
    filt = round((acc/size(area.PixelIdxList,2))/4);    
    color_square_edge_v = bwareaopen(color_square_edge_v, filt); % Parametro
    CC_edge_v = regionprops(color_square_edge_v,'Centroid','Image');
    CC_v = bwconncomp(color_square_edge_v);
% 
%     xL = color_square_bb(1);
%     xR = color_square_bb(1)+color_square_bb(3);
    vect_left = zeros(5,size(CC_edge_v,1));
    vect_right = zeros(5,size(CC_edge_v,1));
    
    for count=1:size(CC_edge_v,1)   
        xC = CC_edge_v(count).Centroid(1);
        yC = CC_edge_v(count).Centroid(2);
        xL = rct_line_left(yC);
        xR = rct_line_right(yC);
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
    vect_right(5,vect_right(1,:)==vect_left(1,1) | vect_right(1,:)==vect_left(1,2)) = 1;
    vect_left(5,1:2) = 1;   
    ccright = [cell2mat(CC_v.PixelIdxList(vect_right(1,1))); cell2mat(CC_v.PixelIdxList(vect_right(1,2)))];
    vect_left(5,vect_left(1,:)==vect_right(1,1) | vect_left(1,:)==vect_right(1,2)) = 1;
    vect_right(5,1:2) = 1;  
    
%     yStd = std([vect_left(4,2), vect_right(4,2)]);
%     yDist = max([vect_left(4,1), vect_left(4,2), vect_right(4,1), vect_right(4,2)])+2*yStd;
%     if(vect_left(4,3) < yDist && vect_left(5,3) ~=1) % Se non l'ho gi� usato
%         ccleft = [ccleft;cell2mat(CC_v.PixelIdxList(vect_left(1,3)))];
%     end
%     if(vect_right(4,3) < yDist && vect_right(5,3) ~=1)
%         ccright = [ccright;cell2mat(CC_v.PixelIdxList(vect_right(1,3)))];
%     end  

    [idx_left,idy_left] = ind2sub(size(color_square),ccleft);
    
    [fitresult_left, ~] = createLineInv(idx_left,idy_left,size(color_square));
    plot(fitresult_left);
    l1 = legend();
    set(l1,'visible','off');
    [idx_right,idy_right]=ind2sub(size(color_square),ccright);
    [fitresult_right, ~] = createLineInv(idx_right,idy_right,size(color_square));
    plot(fitresult_right);
    chess.v_lines = cell(1,2);
    chess.v_lines{1} = fitresult_left;
    chess.v_lines{2} = fitresult_right;
    l1 = legend();
    set(l1,'visible','off');
%% PLOT   
    xx = -40:0.001:size(color_square,2)+40;
    [~,ii] = min(abs(fitresult_bot(xx(:)) - fitresult_left(xx(:))));
    ybotleft = fitresult_bot(xx(ii));
    xbotleft = xx(ii);
    scatter(xbotleft,ybotleft)
    if(xbotleft<1)
        xbotleft = 1; 
    end
    if(xbotleft>size(color_square,2))
        xbotleft = size(color_square,2); 
    end
    if(ybotleft<1)
        ybotleft = 1; 
    end
    if(ybotleft>size(color_square,1))
        ybotleft = size(color_square,1); 
    end   
    [~,ii] = min(abs(fitresult_bot(xx(:)) - fitresult_right(xx(:))));
    ybotright = fitresult_bot(xx(ii));
    xbotright = xx(ii);
    scatter(xbotright,ybotright)
    if(xbotright<1)
        xbotright = 1; 
    end
    if(xbotright>size(color_square,2))
        xbotright = size(color_square,2); 
    end
    if(ybotright<1)
        ybotright = 1; 
    end
    if(ybotright>size(color_square,1))
        ybotright = size(color_square,1); 
    end   
    [~,ii] = min(abs(fitresult_top(xx(:)) - fitresult_left(xx(:))));
    ytopleft = fitresult_top(xx(ii));
    xtopleft = xx(ii);
    scatter(xtopleft,ytopleft)
    if(xtopleft<1)
        xtopleft = 1; 
    end
    if(xtopleft>size(color_square,2))
        xtopleft = size(color_square,2); 
    end
    if(ytopleft<1)
        ytopleft = 1; 
    end
    if(ytopleft>size(color_square,1))
        ytopleft = size(color_square,1); 
    end   
    [~,ii] = min(abs(fitresult_top(xx(:)) - fitresult_right(xx(:))));
    ytopright = fitresult_top(xx(ii));
    xtopright = xx(ii);
    scatter(xtopright,ytopright)
    if(xtopright<1)
        xtopright = 1; 
    end
    if(xtopright>size(color_square,2))
        xtopright = size(color_square,2); 
    end
    if(ytopright<1)
        ytopright = 1; 
    end
    if(ytopright>size(color_square,1))
        ytopright = size(color_square,1); 
    end   
    
    color_square(round(ybotleft),round(xbotleft)) = 1;
    color_square(round(ybotright),round(xbotright)) = 1;
    color_square(round(ytopleft),round(xtopleft)) = 1;
    color_square(round(ytopright),round(xtopright)) = 1;
     
    idx = find(color_square == 1);
    [idx,idy]=ind2sub(size(color_square),idx);
    j = boundary(idx,idy,0.1); % Parametro
    color_square = poly2mask(idy(j),idx(j), size(color_square,1), size(color_square,2));  
end

