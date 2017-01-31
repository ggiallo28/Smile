function [left_center_axis, right_center_axis, mid_center_axis] = generate_central_axis_two(Container)
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
    WHITEMASK = false(size(FullMask));
    BLACKMASK = false(size(FullMask));
    for i=1:size(Container.obj_chess,1)
        if ( ~Container.obj_chess(i).isEmpty )
            for j=1:size(Container.obj_chess(i).chess,2)
                switch (Container.obj_chess(i).chess(j).background{1})
                    case 'White'
                        WHITEMASK = WHITEMASK | Container.obj_chess(i).chess(j).ch_mask;
                    case 'Black'
                        BLACKMASK = BLACKMASK | Container.obj_chess(i).chess(j).ch_mask;
                end
            end
        end
    end
    WHITEMASK = imclose(WHITEMASK,strel('square',20));
    BLACKMASK = imclose(BLACKMASK,strel('square',20));
    left = imdilate(edge(WHITEMASK),strel('disk',7));
    right = imdilate(edge(BLACKMASK),strel('disk',7));
    center_axis = left & right;
    if (~Container.isGUI)
        fig = figure;
        imshow(I); hold on;
        v_legend = [];
    end
    mask = false(size(FullMask));
    if ( isL1 ) 
        mask = imerode(maskL1I,strel('square',10));
    end
    left1_center_axis = center_axis.*mask;
    idx = find(left1_center_axis == 1);
    if ( ~isempty(idx) && ~Container.isGUI)
        [idy,idx] = ind2sub(size(FullMask),idx);
        [left1_fitresult, ~] = createLineInv(idy,idx,size(FullMask));
        plot(left1_fitresult,'b'); 
        v_legend = [v_legend, 'left axis'];
    end
    if ( isL2 ) 
        mask = imerode(maskL2I,strel('square',10));
    end
    left2_center_axis = center_axis.*mask;
    idx = find(left2_center_axis == 1);
    if ( ~isempty(idx) && ~Container.isGUI)
        [idy,idx] = ind2sub(size(FullMask),idx);
        [left2_fitresult, ~] = createLineInv(idy,idx,size(FullMask));
        plot(left2_fitresult,'b'); 
        v_legend = [v_legend, 'left axis'];
    end
    left_center_axis = left2_center_axis | left1_center_axis;

    if ( isR1 ) 
        mask = imerode(maskR1I,strel('square',10));
    end
    right1_center_axis = center_axis.*mask;
    idx = find(right1_center_axis == 1);
    if ( ~isempty(idx) && ~Container.isGUI)
        [idy,idx] = ind2sub(size(FullMask),idx);
        [right1_fitresult, ~] = createLineInv(idy,idx,size(FullMask));
        plot(right1_fitresult,'r'); 
        v_legend = [v_legend, 'right axis'];
    end
    if ( isR2 ) 
        mask = imerode(maskR2I,strel('square',10));
    end
    right2_center_axis = center_axis.*mask;
    idx = find(right2_center_axis == 1);
    if ( ~isempty(idx) && ~Container.isGUI)
        [idy,idx] = ind2sub(size(FullMask),idx);
        [right2_fitresult, ~] = createLineInv(idy,idx,size(FullMask));
        plot(right2_fitresult,'r'); 
        v_legend = [v_legend, 'right axis'];
    end
    right_center_axis = right2_center_axis | right1_center_axis;

    if ( isC )
        mid_center_axis = center_axis.*imerode(maskCI,strel('square',10));
        idx = find(mid_center_axis == 1);
        if ( ~isempty(idx) && ~Container.isGUI)
            [idy,idx] = ind2sub(size(FullMask),idx);
            [mid_fitresult, ~] = createLineInv(idy,idx,size(FullMask));
            plot(mid_fitresult, 'y');
            v_legend = [v_legend, 'center axis'];
        end
    end
   % legend(v_legend);
   % print(fig,[path,'central_axis'],'-dpng')
end