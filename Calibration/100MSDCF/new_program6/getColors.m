function colors = getColors(i,j,obj_chess,w,c)
    % i = riga del pixel
    % j = colonna del pixel
    % obj_chess = oggetto contenente le maschere
    % w = ampiezza finestra
    
    lum = (0.2126 * c(1,1,1)) + (0.7152 * c(1,1,2)) + (0.0722 * c(1,1,3)); % 
    persistent red_h;
    persistent red;
    persistent green_h;
    persistent green;
    persistent blue_h;
    persistent white;
    persistent blue;
    persistent light_blue;
    persistent light_blue_h;
    persistent yellow_h;
    persistent yellow;
    persistent black;
    persistent violet;
    persistent violet_h;

    
    
%     if isempty(red_h) 
%         red_h = rgb2lab(reshape([1 0 0].*0.3,1,1,3)); 
%         red = rgb2lab(reshape([1 0 0],1,1,3)); 
%         green_h = rgb2lab(reshape([0 1 0].*0.3,1,1,3));  
%         green = rgb2lab(reshape([0 1 0],1,1,3)); 
%         blue_h =  rgb2lab(reshape([0 0 1].*0.3,1,1,3)); 
%         blue=  rgb2lab(reshape([0 0 1],1,1,3)); 
%         light_blue_h = rgb2lab(reshape([0 1 1].*0.3,1,1,3)); 
%         light_blue = rgb2lab(reshape([0 1 1],1,1,3)); 
%         yellow_h = rgb2lab(reshape([1 1 0].*0.3,1,1,3)); 
%         yellow = rgb2lab(reshape([1 1 0],1,1,3)); 
%         violet_h = rgb2lab(reshape([1 0 1].*0.3,1,1,3));  
%         violet = rgb2lab(reshape([1 0 1],1,1,3));  
%         white = rgb2lab(reshape([1 1 1],1,1,3)); 
%         black = rgb2lab(reshape([0 0 0],1,1,3)); 
%     end

    if isempty(red_h) 
        red_h = reshape([1 0 0],1,1,3); 
        red = reshape([1 0 0],1,1,3); 
        green_h = reshape([0 1 0],1,1,3);  
        green = reshape([0 1 0],1,1,3); 
        blue_h =  reshape([0 0 1],1,1,3); 
        blue=  reshape([0 0 1],1,1,3); 
        light_blue_h = reshape([0 1 1],1,1,3); 
        light_blue = reshape([0 1 1],1,1,3); 
        yellow_h = reshape([1 1 0],1,1,3); 
        yellow = reshape([1 1 0],1,1,3); 
        violet_h = reshape([1 0 1],1,1,3);  
        violet = reshape([1 0 1],1,1,3);  
        white = reshape([1 1 1],1,1,3); 
        black = reshape([0 0 0],1,1,3); 
    end
    
    colors = [];
    black_added = false;
    
    for l=1:size(obj_chess,1)
        for k=1:size(obj_chess(l).chess,2)
             mask = obj_chess(l).chess(k).ch_mask;
             color_mask = obj_chess(l).color_mask;
             inv_color_mask =  obj_chess(l).inv_color_mask;
             s_wid = round(j-w/2);
             s_hei = round(i-w/2);
             e_wid = round(j+w/2);
             e_hei = round(i+w/2);
             if(s_wid < 1)
                 s_wid = 1;
             end
             if(s_hei < 1)
                 s_hei = 1;
             end
             if(e_wid > size(mask,2))
                 e_wid = size(mask,2);
             end
             if(e_hei > size(mask,1))
                 e_hei = size(mask,1);
             end
             check = mask(s_hei:e_hei,s_wid:e_wid);
             check_color = color_mask(s_hei:e_hei,s_wid:e_wid);
             check_inv_color = inv_color_mask(s_hei:e_hei,s_wid:e_wid);
             if(mask(i,j) | sum(sum(check))>0)
                 if(strcmp(obj_chess(l).chess(k).background,'White'))
                     colors = [colors, white];
                 end
                 if(strcmp(obj_chess(l).chess(k).background,'Black'))
                     if ~black_added
                        colors = [colors, black];
                        black_added = true;
                    end
                 end
                 if(strcmp(obj_chess(l).name,'red'))
                     colors = [colors, red_h, red];
                 end
                 if(strcmp(obj_chess(l).name,'yellow'))
                     colors = [colors, yellow_h, yellow];
                 end
                 if(strcmp(obj_chess(l).name,'green'))
                     colors = [colors, green, green_h];
                 end
                 if(strcmp(obj_chess(l).name,'blue'))
                     colors = [colors, blue_h, blue];
                 end
             else
                 if ~black_added
                    colors = [colors black];
                    black_added = true;
                 end
             end   
        end       
    end
end