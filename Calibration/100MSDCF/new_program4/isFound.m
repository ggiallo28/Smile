function bool = isFound(a,b,dim)
    bool = false;
    for i=1:size(a,1)
        if(round(a(i))==b)
            bool = true;
            return;
        end
    end
end