function [x, y] = compute_projections(BW, split) 
    x = []; x(1,1) = split(1); k = 1;
    for i=2:size(split,2)
        if(abs(split(i-1)-split(i))>1)
            x(k,2) = split(i-1);
            x(k+1,1) = split(i);
            k = k+1;
        end
    end
    x(k,2) = split(end); y = [];
    for i = 1:size(x,1)
        bw_cut = BW(:,x(i,1):x(i,2),:);
        bw_cut = sum(bw_cut');
        idx = find(bw_cut~=0);
        y = [y; min(idx),max(idx)];
    end
end