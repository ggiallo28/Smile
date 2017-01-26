function image = prepareHull(image)
     [LTp, RTp, LBp, RBp]  = getRotatedRectangle(image);
     [LTp, RTp, LBp, RBp]  =  roundPoints([LTp; RTp; LBp; RBp],size(image));
     image(LTp(2),LTp(1)) = 1; image(RTp(2),RTp(1)) = 1; image(LBp(2),LBp(1)) = 1; image(RBp(2),RBp(1)) = 1;
end

