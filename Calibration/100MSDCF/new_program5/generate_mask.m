function [maskC, isC, maskL1, isL1, maskL2, isL2, maskR1, isR1, maskR2, isR2] = generate_mask(obj_chess, label, order, positions, types, dimensions)
    maskC  = false(dimensions(1),dimensions(2));
    maskL1 = false(dimensions(1),dimensions(2));
    maskL2 = false(dimensions(1),dimensions(2));
    maskR1 = false(dimensions(1),dimensions(2));
    maskR2 = false(dimensions(1),dimensions(2));
    isC = false;
    isL1 = false;
    isL2 = false;
    isR1 = false;
    isR2 = false;

    for i=1:size(obj_chess,1)
        for j=1:size(obj_chess(i).chess,2)
            if(strcmp(obj_chess(i).chess(j).position,'Left') && strcmp(obj_chess(i).chess(j).type,'Secondary'))
                maskL2 = (maskL2 | obj_chess(i).chess(j).mask);
                isL2 = true;
            end
            if(strcmp(obj_chess(i).chess(j).position,'Right') && strcmp(obj_chess(i).chess(j).type,'Primary'))
                maskR1 = maskR1 | obj_chess(i).chess(j).mask;
                isR1 = true;
            end
            if(strcmp(obj_chess(i).chess(j).position,'Left') && strcmp(obj_chess(i).chess(j).type,'Primary'))
                maskL1 = maskL1 | obj_chess(i).chess(j).mask;
                isL1 = true;
            end
            if(strcmp(obj_chess(i).chess(j).position,'Right') && strcmp(obj_chess(i).chess(j).type,'Secondary'))
                maskR2 = maskR2 | obj_chess(i).chess(j).mask;
                isR2 = true;
            end
            if(strcmp(obj_chess(i).chess(j).position,'Center'))
                maskC = maskC | obj_chess(i).chess(j).mask;
                isC  = true;
            end
        end
    end
    %% Rifinitura Maschere
    if ( isC )
        % Conto quante linee ci devono essere in funzione dei colori
        idxLinesCenter = find(strcmp(label(3,:),positions(2)));
        numLinesCenter = size(idxLinesCenter,2)+1;
        % Uso il numero teorico di linee per unire quelle troppo vicine
        LinesCenter = findLines(idxLinesCenter,  obj_chess, order);
        LinesCenter = mergeNearestLines(LinesCenter, numLinesCenter, size(maskC));
        gap  = false(dimensions(1),dimensions(2));
        for i=1:size(LinesCenter,2)
            gap = gap | line2image(LinesCenter{i},size(maskC));
        end
        % Separo le checkerboard appartenenti ad una stessa tipologia di immagine
        gap = imdilate(gap,strel('disk',3));
        maskC = maskC & ~gap;
        maskC = imopen(maskC,strel('square',3));

    end
    if ( isL2 )

        idxLinesLeftSec = find(strcmp(label(3,:),positions(1)) & strcmp(label(4,:),types(2)));
        numLinesLeftSec = size(idxLinesLeftSec,2)+1;
        LinesLeftSec = findLines(idxLinesLeftSec,  obj_chess, order);
        LinesLeftSec = mergeNearestLines(LinesLeftSec, numLinesLeftSec, size(maskC));
        gap  = false(dimensions(1),dimensions(2));
        for i=1:size(LinesLeftSec,2)
            gap = gap | line2image(LinesLeftSec{i},size(maskC));
        end
        gap = imdilate(gap,strel('disk',3));
        maskL2 = maskL2 & ~gap;
        maskL2 = imopen(maskL2,strel('square',3));

    end
    if ( isL1 )

        idxLinesLeftPri = find(strcmp(label(3,:),positions(1)) & strcmp(label(4,:),types(1)));
        numLinesLeftPri = size(idxLinesLeftPri,2)+1;
        LinesLeftPri = findLines(idxLinesLeftPri,  obj_chess, order);
        LinesLeftPri = mergeNearestLines(LinesLeftPri, numLinesLeftPri, size(maskC));
        gap  = false(dimensions(1),dimensions(2));
        for i=1:size(LinesLeftPri,2)
            gap = gap | line2image(LinesLeftPri{i},size(maskC));
        end
        gap = imdilate(gap,strel('disk',3));
        maskL1 = maskL1 & ~gap;
        maskL1 = imopen(maskL1,strel('square',3));

    end
    if ( isR2 )

        idxLinesRightSec = find(strcmp(label(3,:),positions(3)) & strcmp(label(4,:),types(2)));
        numLinesRightSec = size(idxLinesRightSec,2)+1;
        LinesRightSec = findLines(idxLinesRightSec,  obj_chess, order);
        LinesRightSec = mergeNearestLines(LinesRightSec, numLinesRightSec, size(maskC));
        gap  = false(dimensions(1),dimensions(2));
        for i=1:size(LinesRightSec,2)
            gap = gap | line2image(LinesRightSec{i},size(maskC));
        end
        gap = imdilate(gap,strel('disk',3));
        maskR2 = maskR2 & ~gap;
        maskR2 = imopen(maskR2,strel('square',3));

    end
    if ( isR1 )

        idxLinesRightPri = find(strcmp(label(3,:),positions(3)) & strcmp(label(4,:),types(1)));
        numLinesRightPri = size(idxLinesRightPri,2)+1;
        LinesRightPri = findLines(idxLinesRightPri,  obj_chess, order);
        LinesRightPri = mergeNearestLines(LinesRightPri, numLinesRightPri, size(maskC));
        gap  = false(dimensions(1),dimensions(2));
        for i=1:size(LinesRightPri,2)
            gap = gap | line2image(LinesRightPri{i},size(maskC));
        end
        gap = imdilate(gap,strel('disk',3));
        maskR1 = maskR1 & ~gap;
        maskR1 = imopen(maskR1,strel('square',3));

    end
    %% Ripulisco maschere da oggetti piccoli
    if ( isL2 )
        statsL2 = cell2mat(struct2cell(regionprops(maskL2,'Area')));
        maskL2 = bwareaopen(maskL2,round(0.1*(mean(statsL2)-std(statsL2))));
        figure, imshow(maskL2);
    end
    if ( isL1 )
        statsL1 = cell2mat(struct2cell(regionprops(maskL1,'Area')));
        maskL1 = bwareaopen(maskL1,round(0.1*(mean(statsL1)-std(statsL1))));
        figure, imshow(maskL1);
    end
    if ( isC )
        statsC = cell2mat(struct2cell(regionprops(maskC,'Area')));
        maskC = bwareaopen(maskC,round(0.1*(mean(statsC)-std(statsC))));
        figure, imshow(maskC);
    end
    if ( isR1 )
        statsR1 = cell2mat(struct2cell(regionprops(maskR1,'Area')));
        maskR1 = bwareaopen(maskR1,round(0.1*(mean(statsR1)-std(statsR1))));
        figure, imshow(maskR1);
    end
    if ( isR2 )
        statsR2 = cell2mat(struct2cell(regionprops(maskR2,'Area')));
        maskR2 = bwareaopen(maskR2,round(0.1*(mean(statsR2)-std(statsR2))));
        figure, imshow(maskR2);
    end
end