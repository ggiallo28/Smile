function bwD = mask_split(BW, inv_BW, x, y, Container)
toplot = true & ~Container.isGUI;
    BW = bwareaopen(BW, 50); % Parametro
    inv_BW = bwareaopen(inv_BW, 50); % Parametro
    BW_edge = edge(BW);
    assert(sum(sum(inv_BW&BW))==0);
    inv_BW_edge = edge(inv_BW);
    eedge = imdilate(BW_edge|inv_BW_edge,strel('square',3));
    eedge = imerode(eedge,strel('square',2));
    BW(eedge) = 0;
    inv_BW(eedge) = 0;
    bwD = BW | inv_BW;
    bwD = bwareaopen(bwD, 50); % Parametro

    for i = 1:size(x,1)
        chess = bwD(y(i,1):y(i,2),x(i,1):x(i,2));
        condition = true;
if toplot | Container.isGUI
        figure
end
        while condition
            CC = bwconncomp(chess);
            bb = regionprops(CC,'BoundingBox'); bboxCorr = zeros(1,size(bb,1));
            for iter=1:size(bb,1)
                BW_TMP = false(size(chess));
                BW_TMP(CC.PixelIdxList{iter}) = 1;
                s = regionprops(BW_TMP,'BoundingBox');
                img = imcrop(BW_TMP,s.BoundingBox);
                img = imopen(img,strel('square',5));
                toCompare = bwconvhull(img);
                bboxCorr(iter) = corr2(toCompare,imfill(img,'holes'));
                % Provare ad usare la Rectangularity
    %             figure, imshowpair(toCompare,img,'montage');
                bboxCorr(iter) = 1-sum(sum(abs(double(toCompare) - double(img))))/sum(sum(abs(double(img))));
            end
            if(find(bboxCorr<0.65)) %0.72
                idx = find(bboxCorr<0.65);
                tmp_chess = zeros(size(chess));
                for k=1:size(idx)
                    tmp_chess(CC.PixelIdxList{idx(k)}) = 1;
                    chess(CC.PixelIdxList{idx(k)}) = 0;
                    %tmp_chess = imerode(tmp_chess,strel('square',2));
                    tmp_chess = bwdist(~tmp_chess,'chessboard');
                    tmp_chess = tmp_chess>1;
if toplot | Container.isGUI
                    imshowpair(tmp_chess,chess,'falsecolor');
                    pause(0.2);
end
                    chess = chess | tmp_chess;
                    chess = bwareaopen(chess, 100);
                end      
            else
                condition = false;
            end        
        end
        bwD(y(i,1):y(i,2),x(i,1):x(i,2))=chess;
    end
% 
%     L = bwlabel(bwD);
%     RGB = label2rgb(L);
%     imshow(RGB);
close all;
end