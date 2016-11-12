function [obj_chess, transtions] = getColorBoundary(obj_red, obj_green, obj_blue, obj_yellow, checker_vector, dim)
    raw_order = []; raw_label = [];
    obj_chess = [obj_red,obj_green,obj_blue,obj_yellow];
    bw = false(dim(1),dim(2));
    for j=1:size(obj_chess,2);
        for i=1:size(obj_chess(j).chess,2)
            raw_order = [raw_order,[j;i;obj_chess(j).chess(i).raw_centroid']];
            if (strcmp(obj_chess(j).name,'red'))
                if(strcmp(obj_chess(j).chess(i).background,'Black'))
                    raw_label = [raw_label,{'magenta'}];
                else
                    raw_label = [raw_label,{obj_chess(j).name}];
                end
            elseif (strcmp(obj_chess(j).name,'blue'))
                if(strcmp(obj_chess(j).chess(i).background,'Black'))
                    raw_label = [raw_label,{'cyan'}];
                else
                    raw_label = [raw_label,{obj_chess(j).name}];
                end
            else
                raw_label = [raw_label,{obj_chess(j).name}];
            end  
        end
        bw = bw | obj_chess(j).raw_bw;
    end
    split = sum(bw);
    split = find(split~=0);
    if size(split,2) == 0
        obj.isEmpty = true;
        return;
    end
    [x, y] = compute_projections(bw, split);
    [raw_order(3,:), idx] = sort(raw_order(3,:));
    raw_order(1,idx) = raw_order(1,idx);
    raw_order(2,idx) = raw_order(2,idx);
    raw_order(4,idx) = raw_order(4,idx);
    raw_label = raw_label(idx);
    raw_label = [raw_label;raw_label];
    for i=1:size(raw_label,2)
        center_x = raw_order(3,i);
        if(center_x<x(1,2))
            raw_label(2,i) = {'Left'};
        end
        if(center_x>x(3,1))
            raw_label(2,i) = {'Right'};
        end
        if(center_x>x(2,1) && center_x<x(2,2))
            raw_label(2,i) = {'Center'};
        end
    end
    check = checker_vector(2,:,:);
    last = 0;
    idx_curr = raw_order(1:2,1)';
    % -1 = Left
    % 1 = Right
    side = -1;
    transtions = [idx_curr,side];
    side = -1*side;
    for i=1:size(raw_label,2)-1
        curr_color = name2code(raw_label(1,i));
        succ_color = name2code(raw_label(1,i+1));
        ii = find(check(:,:,1) == succ_color(1) & check(:,:,2) == succ_color(2) & check(:,:,3) == succ_color(3));
        jj = find(check(:,:,1) == curr_color(1) & check(:,:,2) == curr_color(2) & check(:,:,3) == curr_color(3));
%         raw_label(1,i)
%         raw_label(1,i+1)
        if(~strcmp(raw_label(2,i),raw_label(2,i+1)) && last ~=0)
            order = -1*last;
        elseif(mod(jj-1,6)==5 && mod(ii-1,6)==0)
            order = 1;
        elseif(mod(jj-1,6)==0 && mod(ii-1,6)==5)
            order = -1;
        elseif(mod(jj-1,6)<mod(ii-1,6))
            order = 1;
        else
            order = -1;
        end
        if(last ~= 0 && last~=order)   
            idx_curr = raw_order(1:2,i)';
            idx_succ = raw_order(1:2,i+1)';
            transtions = [transtions;[idx_curr,side]]; %obj_curr id_chess_curr, obj_succ id_chess_succ
            side = -1*side;
            transtions = [transtions;[idx_succ,side]]; %obj_curr id_chess_curr, obj_succ id_chess_succ
            side = -1*side;
        end
        last = order;
    end
    idx_succ = raw_order(1:2,end)';
    transtions = [transtions;[idx_succ,side]];
end