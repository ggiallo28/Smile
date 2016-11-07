function code = name2code(name)
    if iscell(name)
        name = name{1};
    end
    name = lower(name);
    switch name
        case 'red'
            code = [255,0,0;255,0,255];
        case 'green'
            code = [0,255,0];
        case 'blue'
            code = [0,0,255;0,255,255];
        case 'violet'
            code = [255,0,255;255,0,0];
        case 'yellow'
            code = [255,255,0];
        case 'light_blue'
            code = [0,255,255;0,0,255];        
        case 'white'
            code = [255,255,255];
        case 'black'
            code = [0,0,0];
    end
end

