function Container = generate_bwconvhull(Container)
    maskC = Container.maskC;
    isC = Container.isC;
    maskL1 = Container.maskL1;
    isL1 = Container.isL1;
    maskL2 = Container.maskL2;
    isL2 = Container.isL2;
    maskR1 = Container.maskR1;
    isR1 = Container.isR1;
    maskR2 = Container.maskR2;
    isR2 = Container.isR2;
    
    maskCI = false(size(maskC));
    maskL1I = false(size(maskC));
    maskL2I = false(size(maskC));
    maskR1I = false(size(maskC));
    maskR2I = false(size(maskC));
    for k=1:5
        switch(k)
            case 1
                if (~isC)
                    continue;
                end
                mask = maskC;
            case 2
                if (~isL1)
                    continue;
                end
                mask = maskL1;
            case 3
                if (~isL2)
                    continue;
                end 
                mask = maskL2;
            case 4
                if (~isR1)
                    continue;
                end 
                mask = maskR1;
            case 5
                if (~isR2)
                    continue;
                end 
                mask = maskR2;
        end
        maskI = bwconvhull(mask);
        switch(k)
            case 1
                maskCI = maskI;
            case 2
                maskL1I = maskI;
            case 3
                maskL2I = maskI;
            case 4
                maskR1I = maskI;
            case 5
                maskR2I = maskI;
        end   
    end
    Container.maskCI = maskCI;
    Container.maskL1I = maskL1I;
    Container.maskL2I = maskL2I;
    Container.maskR1I = maskR1I;
    Container.maskR2I = maskR2I;
end