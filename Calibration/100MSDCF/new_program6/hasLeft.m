function bool_value = hasLeft(label,id)
    if(id==1)
        bool_value = false;
        return;
    end
    label_left = label(:,id-1);
    label_curr = label(:,id);
    if( strcmp(label_left(3),label_curr(3)) && strcmp(label_left(4),label_curr(4)) )
       bool_value = true; 
    else
       bool_value = false; 
    end
end
