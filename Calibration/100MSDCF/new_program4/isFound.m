function bool = isFound(a,b,dim)
    [idr,idc]=ind2sub(dim,b);
    tess = convhulln([idc,idr]);
    bool = false;
    for i=1:size(a,1)
        if(round(a(i))==b)
            bool = true;
            return;
        end
    end
end