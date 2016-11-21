function [I, I_BG, O, O_BG, BB] = normalize_image(orig, orig_bg, bb)
    if ~isempty(bb)
       [O, BB] = imcrop(orig, bb); 
    else
        figure, [O, BB] = imcrop(orig);
    end
    R = im2double(O(:,:,1));
    G = im2double(O(:,:,2));
    B = im2double(O(:,:,3));
    L = (0.2126 * R) + (0.7152 * G) + (0.0722 * B); 
    K = zeros(size(O));
    K(:,:,1) = R.*(2-L);
    K(:,:,2) = G.*(2-L);
    K(:,:,3) = B.*(2-L);
    I = im2uint8(K);
    O_BG = imcrop(orig_bg,BB);
    Rt = im2double(O_BG(:,:,1));
    Gt = im2double(O_BG(:,:,2));
    Bt = im2double(O_BG(:,:,3));
    L = (0.2126 * Rt) + (0.7152 * Gt) + (0.0722 * Bt); 
    K = zeros(size(O));
    K(:,:,1) = Rt.*(2-L);
    K(:,:,2) = Gt.*(2-L);
    K(:,:,3) = Bt.*(2-L);
    I_BG = im2uint8(K);
end