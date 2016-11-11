function [ output_args, shift ] = allignleft( input_args )
    output_args = false(size(input_args));
    input_args_copy = input_args;
    id = find(input_args_copy(1,:) == 0);
    input_args_copy(:,id) = [];
    output_args(:,1:size(input_args_copy,2)) = input_args_copy;
    shift = 0;
    for i=1:size(input_args,2)
        if(input_args(1,i) == 0)
            shift = shift +1;
        else
            break;
        end
    end
end

