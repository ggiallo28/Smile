function [ output_args ] = UnsharpMasking( input_image, sigma)
input_image = double(input_image);
Blurred = imgaussfilt(input_image,sigma);
if(size(input_image,3) == 3)
    R = Blurred(:,:,1); G = Blurred(:,:,2); B = Blurred(:,:,3);
    var_pl1 = var(R(:));
    var_pl2 = var(G(:));
    var_pl3 = var(B(:));
    Rin = input_image(:,:,1); Gin = input_image(:,:,2); Bin = input_image(:,:,3);
    ms1pl_cov = cov(Rin(:),R(:));
    ms2pl_cov = cov(Gin(:),G(:));
    ms3pl_cov = cov(Bin(:),B(:));
    G1 = ms1pl_cov(1,2)/var_pl1;
    G2 = ms2pl_cov(1,2)/var_pl2;
    G3 = ms3pl_cov(1,2)/var_pl3;
    output_args = zeros(size(input_image));
    output_args(:,:,1) = input_image(:,:,1) + ( input_image(:,:,1) - Blurred(:,:,1)) .* G1;
    output_args(:,:,2) = input_image(:,:,2) + ( input_image(:,:,2) - Blurred(:,:,2)) .* G2;
    output_args(:,:,3) = input_image(:,:,3) + ( input_image(:,:,3) - Blurred(:,:,3)) .* G3;
else
    var_pl = var(Blurred(:));
    mspl_cov = cov(input_image(:),Blurred(:));
    G = mspl_cov(1,2)/var_pl;
    output_args = input_image + ( input_image - Blurred) .* G;
end
    output_args = uint8(output_args);
end