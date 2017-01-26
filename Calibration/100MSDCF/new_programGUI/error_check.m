function obj_chess = error_check(obj_red, obj_green, obj_blue, obj_yellow)
    for s=1:2
        switch(s)
            case 1
                obj2split = obj_blue;
            case 2
                obj2split = obj_red;
        end
        for k=1:2   
            switch(str2double([num2str( s ),num2str( k )]))
                case 11 
                    obj2use = objBlobs('cyan');
                    color2find = 'Black';
                    color_value = [0, 255, 255];
                case 12
                    obj2use = objBlobs('blue');
                    color2find = 'White';
                    color_value = [0, 0, 255];
                case 21
                    obj2use = objBlobs('magenta');
                    color2find = 'Black';
                    color_value = [255, 0, 255];
                case 22
                    obj2use = objBlobs('red');
                    color2find = 'White';
                    color_value = [255, 0, 0];
            end
            idk = 1;
            for i=1:size(obj2split.chess,2)
                if(strcmp(obj2split.chess(i).background,color2find))
                    obj2use.isEmpty = false;
                    obj2use.chess(idk) = objChess();
                    obj2use.chess(idk) = obj2split.chess(i);
                    x_split = obj2split.bbox_x(i,:);
                    y_split = obj2split.bbox_y(i,:);
                    obj2use.bbox_x(idk,1:2) = x_split;
                    obj2use.bbox_y(idk,1:2) = y_split;
                    if(obj2use.true_color_mask == 0)
                        obj2use.true_color_mask = false(size(obj2split.color_mask));
                    end
                    obj2use.true_color_mask(y_split(1):y_split(2),x_split(1):x_split(2)) =...
                        obj2split.true_color_mask(y_split(1):y_split(2),x_split(1):x_split(2));
                    if(obj2use.color_mask == 0)
                        obj2use.color_mask = false(size(obj2split.color_mask));
                    end
                    obj2use.color_mask(y_split(1):y_split(2),x_split(1):x_split(2)) =...
                        obj2split.color_mask(y_split(1):y_split(2),x_split(1):x_split(2));
                    if(obj2use.inv_color_mask == 0)
                        obj2use.inv_color_mask = false(size(obj2split.color_mask));
                    end
                    obj2use.inv_color_mask(y_split(1):y_split(2),x_split(1):x_split(2)) =...
                         obj2split.inv_color_mask(y_split(1):y_split(2),x_split(1):x_split(2));
                    if(obj2use.black_white == 0)
                        obj2use.black_white = false(size(obj2split.color_mask));
                    end 
                    obj2use.black_white(y_split(1):y_split(2),x_split(1):x_split(2)) =...
                         obj2split.black_white(y_split(1):y_split(2),x_split(1):x_split(2));
                    if(obj2use.masked_rgb == 0)
                        obj2use.masked_rgb = uint8(zeros(size(obj2split.color_mask,1),size(obj2split.color_mask,2),3));
                    end
                    Rband = obj2use.masked_rgb(:,:,1); Gband = obj2use.masked_rgb(:,:,2); Bband = obj2use.masked_rgb(:,:,3);
                    Rband(obj2use.color_mask) = color_value(1); Gband(obj2use.color_mask) = color_value(2); Bband(obj2use.color_mask) = color_value(3);
                    obj2use.masked_rgb(:,:,1)=Rband; obj2use.masked_rgb(:,:,2)=Gband; obj2use.masked_rgb(:,:,3)=Bband; 
                    idk = idk +1;
                end
            end      
            switch(str2double([num2str( s ),num2str( k )]))
                case 11
                    obj_ciano = obj2use;
                case 12
                    obj_blue_splitted = obj2use;
                case 21
                    obj_magenta = obj2use;
                case 22
                    obj_red_splitted = obj2use;
            end
        end
    end
    obj_chess = [obj_red_splitted;obj_green;obj_blue_splitted;obj_yellow;obj_ciano;obj_magenta];
end