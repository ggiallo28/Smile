function [rj, cj] = findIntersections(horizontalCell, verticalCell, mask)
    vect_h = zeros(2,size(horizontalCell,2));
    for i=1:size(horizontalCell,2)
        if(~isempty(horizontalCell{i}))
            horizontal_image = false(size(mask));
            coeffs = coeffvalues(horizontalCell{i});
            x_t = 1:size(mask,2);
            y_t = round(polyval(coeffs,x_t));
            x_t(y_t<1 | y_t>size(mask,1)) = [];
            y_t(y_t<1 | y_t>size(mask,1)) = [];
            for j = 1:size(x_t,2)
                horizontal_image(y_t(j),x_t(j)) = 1;
            end
            horizontal_image = horizontal_image.*mask;
            w = sum(horizontal_image,2);
            idx = mean(find(w ~=0));
            vect_h(1,i) = idx;
            vect_h(2,i) = i;
        end
    end
    vect_h(:,vect_h(1,:)==0) = [];
    % Ordering line, in prima posizione quella più in alto
    for l=1:size(vect_h,2)-1
        for i=l:size(vect_h,2)
            if(vect_h(1,l)>vect_h(1,i))
                tmp = vect_h(:,l);
                vect_h(:,l) = vect_h(:,i);
                vect_h(:,i) = tmp;
            end
        end
    end
    vect_v = zeros(2,size(verticalCell,2));
    for i=1:size(verticalCell,2)
        if(~isempty(verticalCell{i}))
            vertical_image = false(size(mask));
            coeffs = coeffvalues(verticalCell{i});
            y_t = 1:0.01:size(mask,1);
            x_t = polyval(coeffs,y_t);
            for j = 1:size(x_t,2)
                vertical_image(round(y_t(j)),round(x_t(j))) = 1;
            end
            vertical_image = vertical_image.*mask;
            w = sum(vertical_image,1);
            idx = mean(find(w ~=0));
            vect_v(1,i) = idx;
            vect_v(2,i) = i;
        end
    end
    vect_v(:,vect_v(1,:)==0) = [];
    % Ordering line, in prima posizione quella più a sinistra
    for l=1:size(vect_v,2)-1
        for i=l:size(vect_v,2)
            if(vect_v(1,l)>vect_v(1,i))
                tmp = vect_v(:,l);
                vect_v(:,l) = vect_v(:,i);
                vect_v(:,i) = tmp;
            end
        end
    end
    rj = zeros(size(vect_h,2),size(vect_v,2));
    cj = zeros(size(vect_h,2),size(vect_v,2));
    for i=1:size(vect_h,2)
        for j=1:size(vect_v,2)
            id_horizontal = vect_h(2,i);
            id_vertical = vect_v(2,j);
            coeffs_vertical = coeffvalues(verticalCell{id_vertical});
            coeffs_horizontal = coeffvalues(horizontalCell{id_horizontal});
            y_t1 = 1:0.01:size(mask,1);
            x_t1 = polyval(coeffs_vertical,y_t1);
            y_t2 = polyval(coeffs_horizontal,x_t1);
            yy = abs(y_t1-y_t2);
            id = find(yy==min(yy));
            cj(i,j)  = x_t1(id);
            rj(i,j) = y_t1(id);
        end
    end
end