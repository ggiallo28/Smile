function [X, Y] = fix_one_miss(in_x, in_y, mask)
    idx = find(mask(:,1)==1);
    v = [];
    for i=1:size(in_x,2)
        in_x(idx,:) = circshift(in_x(idx,:)',1)';
        [~,idc] = find(in_x == 0,1);
        in_x_shifted = in_x; in_x_shifted(:,idc)=[];
        m = 0;
        for k=1:size(in_x_shifted,2)
            col = in_x_shifted(:,k);
            col = abs(diff(col));
            m = m+mean(col);
        end
        m = m/size(in_x_shifted,2);
        v = [v,m];
    end
    [~,shift_amount] = min(v);
    in_x(idx,:) = circshift(in_x(idx,:)',shift_amount)';
    [~,idc] = find(in_x == 0,1);
    in_x(:,idc) = [];
    in_y(idx,:) = circshift(in_y(idx,:)',shift_amount)';
    [~,idc] = find(in_y == 0,1);
    in_y(:,idc) = [];
    X = in_x; Y = in_y;
end