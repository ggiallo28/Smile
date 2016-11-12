function [left_center_axis, right_center_axis, mid_center_axis] = generate_central_axis(Container);
%% INIT
    I = Container.I;
    obj_chess = Container.obj_chess;
    FullMask = Container.FullMask;
    maskCI = Container.maskCI;
    isC = Container.isC;
    maskL1I = Container.maskL1I;
    isL1 = Container.isL1;
    maskL2I = Container.maskL2I;
    isL2 = Container.isL2;
    maskR1I = Container.maskR1I;
    isR1 = Container.isR1;
    maskR2I = Container.maskR2I;
    isR2 = Container.isR2;
%% LOGIC

    BlurredMask = imgaussfilt(double(FullMask),10);
    R = im2double(I(:,:,1)); G = im2double(I(:,:,2)); B = im2double(I(:,:,3));
    MINRGB = min(R,G); MINRGB = min(MINRGB,B);
    MINMASK = false(size(MINRGB));
    for i=1:size(obj_chess,1)
        for j=1:size(obj_chess(i).chess,2) 
                if(~strcmp(obj_chess(i).chess(j).background,'Black'))
                    tmp = obj_chess(i).chess(j).mask & obj_chess(i).color_mask;
                    MINMASK = MINMASK | tmp;
                end
        end
    end
    MINRGB = imgaussfilt((MINRGB+MINMASK).*BlurredMask,2);
    T = graythresh(MINRGB);
    bw = im2bw(MINRGB,T);
    area = cell2mat(struct2cell(regionprops(bwconncomp(bw,8),'Area')));
    RIGHTMASK = bwareaopen(bw, round(0.13*mean(area))); % parametro
    RIGHTMASK = imopen(RIGHTMASK,strel('square',10));
    for i=1:3
        switch(i)
            case 1
                if (~isC)
                    continue;
                end
                mask = maskCI;
            case 2
                if (~isL1 && ~isL2)
                    continue;
                end
                mask = false(size(MINRGB));
                if ( isL1 )
                    mask = mask | maskL1I;
                end
                if ( isL2 )
                    mask = mask | maskL2I;
                end
            case 3
                if (~isR1 && ~isR2)
                    continue;
                end
                mask = false(size(MINRGB));
                if ( isR1 )
                    mask = mask | maskR1I;
                end
                if ( isR2 )
                    mask = mask | maskR2I;
                end
        end
         RIGHTMASK = RIGHTMASK | bwconvhull(RIGHTMASK.*mask);
    end
    LEFTMASK = FullMask-RIGHTMASK;
    LEFTMASK(LEFTMASK<0) = 0;
    LEFTMASK = imopen(LEFTMASK,strel('rectangle',[10,30]));
    LEFTMASK = imclose(LEFTMASK,strel('square',20));
    area = cell2mat(struct2cell(regionprops(bwconncomp(LEFTMASK,8),'Area')));
    LEFTMASK = bwareaopen(LEFTMASK, round(0.13*mean(area))); % parametro
    left = imdilate(edge(LEFTMASK),strel('disk',3));
    right = imdilate(edge(RIGHTMASK),strel('disk',3));
    center_axis = left & right;

    imshow(I); hold on;
    v_legend = [];
    mask = false(size(MINRGB));
    if ( isL1 ) 
        mask = mask | maskL1I;
    end
    if ( isL2 ) 
        mask = mask | maskL2I;
    end
    left_center_axis = center_axis.*mask;
    idx = find(left_center_axis == 1);
    if ( ~isempty(idx) )
        [idy,idx] = ind2sub(size(RIGHTMASK),idx);
        [left_fitresult, left_gof] = createLineInv(idy,idx,size(FullMask));
        plot(left_fitresult,'b'); 
        v_legend = [v_legend, 'left axis'];
    end

    mask = false(size(MINRGB));
    if ( isR1 ) 
        mask = mask | maskR1I;
    end
    if ( isR2 ) 
        mask = mask | maskR2I;
    end
    right_center_axis = center_axis.*mask;
    idx = find(right_center_axis == 1);
    if ( ~isempty(idx) )
        [idy,idx] = ind2sub(size(RIGHTMASK),idx);
        [right_fitresult, right_gof] = createLineInv(idy,idx,size(FullMask));
        plot(right_fitresult,'r'); 
        v_legend = [v_legend, 'right axis'];
    end

    if ( isC )
        mid_center_axis = center_axis.*(maskCI);
        idx = find(mid_center_axis == 1);
        [idy,idx] = ind2sub(size(RIGHTMASK),idx);
        [mid_fitresult, mid_gof] = createLineInv(idy,idx,size(FullMask));
        plot(mid_fitresult, 'y');
        v_legend = [v_legend, 'center axis'];
    end
    legend(v_legend);
end