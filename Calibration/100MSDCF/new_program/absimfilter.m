function [ output_args,  img_filtered_x, img_filtered_y] = absimfilter(img, x_filter)
    img_filtered_x = imfilter(img,x_filter);
    neg = img_filtered_x; neg(neg>0) = 0; neg = abs(neg);
    pos = img_filtered_x; pos(pos<0) = 0; img_filtered_x = pos+neg;
    y_filter = x_filter';
    img_filtered_y = imfilter(img,y_filter);
    neg = img_filtered_y; neg(neg>0) = 0; neg = abs(neg);
    pos = img_filtered_y; pos(pos<0) = 0; img_filtered_y = pos+neg;
    output_args = sqrt(img_filtered_y.^2 + img_filtered_x.^2);
end

