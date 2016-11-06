function [rj, cj] = orderPoints(P, horizontalCell, verticalCell, mask)
    vect_h = zeros(2,size(horizontalCell,2));
%    figure, imshow(mask), hold on;
    for i=1:size(horizontalCell,2)
        if(~isempty(horizontalCell{i}))
            horizontal_image = line2image(horizontalCell{i},size(mask));
            horizontal_image = horizontal_image.*mask;
            w = sum(horizontal_image,2);
            idx = mean(find(w ~=0));
            vect_h(1,i) = idx;
            vect_h(2,i) = i;
%            plot(horizontalCell{i})
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
            vertical_image = line2image(verticalCell{i},size(mask));
            vertical_image = vertical_image.*mask;
            w = sum(vertical_image,1);
            idx = mean(find(w ~=0));
            vect_v(1,i) = idx;
            vect_v(2,i) = i;
%            plot(verticalCell{i});
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
            f = verticalCell{id_vertical};
            f1 = horizontalCell{id_horizontal};
            proj = sum(mask)~=0;
            start = find(proj==1,1);
            until = start+sum(proj);
            xx = start:0.01:until;
            [~,ii] = min(f(xx(:)) - f1(xx(:)));
            try
              x_val = fzero(@(x)f(x) - f1(x),[0,xx(ii)]);
            catch
              x_val = fzero(@(x)f(x) - f1(x),[xx(ii),xx(end)]);
            end   
            cj(i,j)  = x_val;
            rj(i,j) = f(x_val);
%            scatter(cj(i,j),rj(i,j));
        end
    end
    P = [P,zeros(size(P,1),1)];
    for i=1:size(rj,1)
        for j=1:size(rj,2)
            D = 1000*ones(size(P,1),1);
            for k=1:size(P,1)
                if ~P(k,3)
                    D(k) = pdist([cj(i,j), rj(i,j); P(k,2), P(k,1)],'euclidean');
                end
            end
            [~, ii] = min(D);
            cj(i,j) = P(ii,2);
            rj(i,j) = P(ii,1);
            P(ii,3) = 1;
        end
    end
end