function [LTp, RTp, LBp, RBp, rct]  = getRotatedRectangle(color_square)
    % https://kyamagu.github.io/mexopencv/matlab/minAreaRect.html
    [idr,idc] = ind2sub(size(color_square),find(imfill(color_square,'holes')==1));
    rct = cv.minAreaRect([idr,idc]);
    %ROTATED RECTANGLE
    Cx = rct.center(2); Cy = rct.center(1); %the coordinates of your center point in world coordinates
    W = rct.size(2); % the width of your rectangle
    H = rct.size(1); % the height of your rectangle
    theta = -deg2rad(rct.angle); % the angle you wish to rotate
    %The offset of a corner in local coordinates (i.e. relative to the pivot point) (which corner will depend on the coordinate reference system used in your environment)
    Ox = W / 2;
    Oy = H / 2;
    % The rotated position of this corner in world coordinates 
%     imshow(color_square); hold on;
    P1 = [Cx + (Ox  * cos(theta)) - (Oy * sin(theta)), Cy + (Ox  * sin(theta)) + (Oy * cos(theta))];
    P2 = [Cx - (Ox  * cos(theta)) - (Oy * sin(theta)), Cy - (Ox  * sin(theta)) + (Oy * cos(theta))];
    P3 = [Cx - (Ox  * cos(theta)) + (Oy * sin(theta)), Cy - (Ox  * sin(theta)) - (Oy * cos(theta))];
    P4 = [Cx + (Ox  * cos(theta)) + (Oy * sin(theta)), Cy + (Ox  * sin(theta)) - (Oy * cos(theta))];
%     plot([P1(1) P2(1)],[P1(2) P2(2)]); l1 = legend(); set(l1,'visible','off');  
%     plot([P1(1) P4(1)],[P1(2) P4(2)]); l1 = legend(); set(l1,'visible','off');  
%     plot([P4(1) P3(1)],[P4(2) P3(2)]); l1 = legend(); set(l1,'visible','off');  
%     plot([P3(1) P2(1)],[P3(2) P2(2)]); l1 = legend(); set(l1,'visible','off');  
    [LTp, RTp, LBp, RBp] = sortPoints(P1,P2,P3,P4,Cy);  
end

