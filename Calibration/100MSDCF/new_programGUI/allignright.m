function [ output_args, shift ] = allignright( input_args )
    output_args = false(size(input_args));
    input_args_copy = input_args;
    id = find(input_args_copy(1,:) == 0);
    input_args_copy(:,id) = [];
    start = size(input_args,2) - size(input_args_copy,2)+1;
    output_args(:,start:end) = input_args_copy;
    shift = start-1;
end

