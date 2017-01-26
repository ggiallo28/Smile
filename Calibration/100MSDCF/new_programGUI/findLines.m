function Lines = findLines( idxLines,  obj_chess, order)
    Lines = cell(1,1); k = 1;
    for i=1:size(idxLines,2)
        idx_obj_chess = order(1,idxLines(i));
        idx_chess = order(2,idxLines(i)); 
        Lines{k} = obj_chess(idx_obj_chess).chess(idx_chess).v_lines{1}; k = k+1;
        Lines{k} = obj_chess(idx_obj_chess).chess(idx_chess).v_lines{2}; k = k+1;
    end
end

