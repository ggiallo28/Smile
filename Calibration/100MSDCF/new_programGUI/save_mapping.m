function points_math = save_mapping(Container)
    obj_chess = Container.obj_chess;
    points_math = cell(1,5);    
    idp = 1;
    
    for l=1:size(obj_chess,1)
        if ( ~obj_chess(l).isEmpty )
            for i=1:size(obj_chess(l).chess(1).intersections_x,1)
                for j=1:size(obj_chess(l).chess(1).intersections_x,2)
                    for k =1:size(obj_chess(l).chess,2)
                        if(strcmp(obj_chess(l).chess(k).position,'Left') && strcmp(obj_chess(l).chess(k).type,'Secondary'))
                            points_math(idp,1) = mat2cell([obj_chess(l).chess(k).intersections_x(i,j), obj_chess(l).chess(k).intersections_y(i,j)],1,2);
                        end
                        if(strcmp(obj_chess(l).chess(k).position,'Right') && strcmp(obj_chess(l).chess(k).type,'Primary'))
                            points_math(idp,4) = mat2cell([obj_chess(l).chess(k).intersections_x(i,j), obj_chess(l).chess(k).intersections_y(i,j)],1,2);
                        end
                        if(strcmp(obj_chess(l).chess(k).position,'Left') && strcmp(obj_chess(l).chess(k).type,'Primary'))
                            points_math(idp,2) = mat2cell([obj_chess(l).chess(k).intersections_x(i,j), obj_chess(l).chess(k).intersections_y(i,j)],1,2);
                        end
                        if(strcmp(obj_chess(l).chess(k).position,'Right') && strcmp(obj_chess(l).chess(k).type,'Secondary'))
                            points_math(idp,5) = mat2cell([obj_chess(l).chess(k).intersections_x(i,j), obj_chess(l).chess(k).intersections_y(i,j)],1,2);
                        end
                        if(strcmp(obj_chess(l).chess(k).position,'Center'))
                            points_math(idp,3) = mat2cell([obj_chess(l).chess(k).intersections_x(i,j), obj_chess(l).chess(k).intersections_y(i,j)],1,2);
                        end
                    end
                    idp = idp + 1;
                end
            end
        end
    end
end