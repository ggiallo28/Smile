function [obj_chess, transitions, Container] = getColorBoundary(obj_red, obj_green, obj_blue, obj_yellow, checker_vector, Container)
    dim = Container.img_dim;
    
    raw_order = []; raw_label = [];
    obj_chess = [obj_red,obj_green,obj_blue,obj_yellow];
    bw = false(dim(1),dim(2));
    for j=1:size(obj_chess,2)
        for i=1:size(obj_chess(j).chess,2)
            if ( ~obj_chess(j).isEmpty )
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
        end
        bw = bw | obj_chess(j).raw_bw;
    end
    bw = bwareaopen(imdilate(bw,strel('square',20)),5*Container.Threshold);
    split = medfilt2(sum(bw),[1,10]);
    split = find(split~=0);
    if size(split,2) == 0
        obj.isEmpty = true;
        return;
    end
    [x, y] = compute_projections(bw, split);
    assert(size(x,1) == 3, 'Non riesco a dividere in Sinistra Centro e Destra');
    [raw_order(3,:), idx] = sort(raw_order(3,:));
    raw_order(1,:) = raw_order(1,idx);
    raw_order(2,:) = raw_order(2,idx);
    raw_order(4,:) = raw_order(4,idx);
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
    transitions = [idx_curr,side];
    transitions_label = {}; idt = 1;
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
            transitions = [transitions;[idx_curr,side]]; %obj_curr id_chess_curr, obj_succ id_chess_succ
            transitions_label(idt,1) = raw_label(2,i);
            side = -1*side;
            transitions = [transitions;[idx_succ,side]]; %obj_curr id_chess_curr, obj_succ id_chess_succ
            transitions_label(idt,2) = raw_label(2,i+1);
            side = -1*side;
            idt = idt +1;
        end
        last = order;
    end
    idx_succ = raw_order(1:2,end)';
    transitions = [transitions;[idx_succ,side]];
    vector_sizes = [];
    for i=1:size(transitions_label,1)
        if(strcmp(transitions_label(i,1),'Left') && strcmp(transitions_label(i,2),'Left'))
            if i == 1
                vector_sizes = [vector_sizes;Container.size_small];
            end
            vector_sizes = [vector_sizes;Container.size_small;Container.size_medium];
        end
        if(strcmp(transitions_label(i,1),'Left') && strcmp(transitions_label(i,2),'Center'))
            if i == 1
                vector_sizes = [vector_sizes;Container.size_medium];
            end
            vector_sizes = [vector_sizes;Container.size_medium;Container.size_big];
        end
        if(strcmp(transitions_label(i,1),'Center') && strcmp(transitions_label(i,2),'Right'))
            vector_sizes = [vector_sizes;Container.size_big;Container.size_medium];
            if i == size(transitions_label,1)
                vector_sizes = [vector_sizes;Container.size_medium];
            end
        end
        if(strcmp(transitions_label(i,1),'Right') && strcmp(transitions_label(i,2),'Right'))
            vector_sizes = [vector_sizes;Container.size_medium;Container.size_small];
            if i == size(transitions_label,1)
                vector_sizes = [vector_sizes;Container.size_small];
            end
        end
    end
    transitions = [transitions,vector_sizes];
    Container.raw_order = raw_order;
    Container.raw_label = raw_label;
end