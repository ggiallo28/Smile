function fuse = show_result(Container, toSave, path)
%% INIT
    obj_chess = Container.obj_chess;
    I = Container.I;
%% LOGIN
    inv_BW = false(size(obj_chess(1).black_white)); col_BW = inv_BW; dual_BW = col_BW;
    colors_fuse = uint8(zeros(size(obj_chess(1).masked_rgb)));
    for l=1:size(obj_chess,1)
        inv_BW = inv_BW | obj_chess(l).black_white;
        col_BW = col_BW | obj_chess(l).color_mask;
        dual_BW = dual_BW| obj_chess(l).inv_color_mask;
        colors_fuse = colors_fuse + obj_chess(l).masked_rgb;
    end
    BW = uint8(255.*(inv_BW-col_BW));
    fuse = colors_fuse+repmat(BW,[1 1 3]);
    fig = figure; imshow([I;fuse]); hold on;
    for l=1:size(obj_chess,1)
        obj = obj_chess(l);
        if( ~obj.isEmpty )
            for i = 1:size(obj.chess,2)
                for j=1:size(obj.chess(i).center_x,2)
                    for k=1:size(obj.chess(i).center_x,1)
                        scatter(obj.chess(i).center_x(k,j),obj.chess(i).center_y(k,j)) % Riferimento Y verso il basso
                    end
                    plot(obj.chess(i).v_lines_centroid{j});
                    l1 = legend();
                    set(l1,'visible','off');
                end       
            end
        end
    end
    hold off
    if (toSave)
        print(fig,[path,'MySavedPlot'],'-dpng')
    end
    hold off;
end