function num = axisdistance(color_fore,color_back,checher_vector)
    code_fore = name2code(color_fore);
    code_back = name2code(color_back);
    assert(size(code_fore,1)<=2 & size(code_back,1)<=2,'Indeterminazione nei colori')
    assert((size(code_fore,1) + size(code_back,1))<4,'Indeterminazione nei colori'); 
    % Recupero il colore di foreground usando quello di background
    if(size(code_fore,1) == 2)
        if(sum(code_back)==765) % 255+255+255 = allora il background è bianco
            idx = find(sum(code_fore,2) == min(sum(code_fore,2)));
            code_fore = code_fore(idx,:);
        end
        if(sum(code_back)==0)   % 0+0+0 = allora il background è nero
            idx = find(sum(code_fore,2) == max(sum(code_fore,2)));
            code_fore = code_fore(idx,:);
        end
    end
    pos_zero = find(sum(sum(checher_vector,3))==1152); % Somma per il valore centrale
    color_checher_vector = checher_vector(2,:,:);
    color_checher_vector = [color_checher_vector(:,:,1);color_checher_vector(:,:,2);color_checher_vector(:,:,3)];
    pos_color = find(color_checher_vector(1,:) == code_fore(1) & color_checher_vector(2,:) == code_fore(2) & color_checher_vector(3,:) == code_fore(3));
    num = pos_color-pos_zero;
end

