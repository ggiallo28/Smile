function num = axisdistance(color_fore,checher_vector)
    code_fore = name2code(color_fore);
    pos_zero = find(sum(sum(checher_vector,3))==1152); % Somma per il valore centrale
    color_checher_vector = checher_vector(2,:,:);
    color_checher_vector = [color_checher_vector(:,:,1);color_checher_vector(:,:,2);color_checher_vector(:,:,3)];
    pos_color = find(color_checher_vector(1,:) == code_fore(1) & color_checher_vector(2,:) == code_fore(2) & color_checher_vector(3,:) == code_fore(3));
    num = pos_color-pos_zero;
end

