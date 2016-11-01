function output_args = allignleftdouble( input_args )
    output_args = zeros(size(input_args));
    input_args_copy = input_args;
    id = find(input_args_copy(1,:) == 0);
    input_args_copy(:,id) = [];
    output_args(:,1:size(input_args_copy,2)) = input_args_copy;
end

