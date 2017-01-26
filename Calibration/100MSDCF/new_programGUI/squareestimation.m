function [ mean_square, size_idx, mean_big, mean_medium, mean_small  ] = squareestimation( bw, a_or_b )
if a_or_b
    idx = sort(cell2mat(struct2cell(regionprops(bwconncomp(~bw,4),'Area'))),'descend')';
    idx = idx(2:end); % Rimuovo la prima, enorme, è sicuramente background
else
    idx = sort(cell2mat(struct2cell(regionprops(bwconncomp(bw,4),'Area'))),'descend')';
    idx = sort(cell2mat(struct2cell(regionprops(bwconncomp(bw,4),'Area'))),'descend' )';
    remove = round(size(idx,1)*0.30);
    idx = idx(1:end-remove);
end
    idx(idx<mean(idx)-0.8*std(idx) & idx<mean(idx)+0.8*std(idx)) = [];
    idy = 1:size(idx,1);
    idy = idy';
    fitresult = createLine(idy,idx);
    outliers = []; toCont = 1;
    opts = fitoptions( 'Method', 'LinearLeastSquares' , 'Exclude',outliers); ft = fittype( 'poly1' );
    while toCont >0
        fdata = feval(fitresult,idy);
        I = abs(fdata - idx) > 1.5*std(idx);
        toCont = size(find(I>0),1);
        outliers = excludedata(idy,idx,'indices',I);
        [fitresult, ~] = fit( idy, idx, ft, opts );
        plot(fitresult,'r-',idy,idx,'k.',outliers,'m*');
        idy(outliers) = [] ; idx(outliers) = [];
        pause(0.01);
    end
    step_dim_a = floor(size(idx)/3);
    mean_square = mean(idx);
    size_idx = size(idx,1);
    mean_big = mean(fitresult(idy(1:step_dim_a)));
    mean_medium = mean(fitresult(idy(step_dim_a+1:2*step_dim_a)));
    mean_small = mean(fitresult(idy(2*step_dim_a+1:3*step_dim_a)));
    RRGGBB = double(cat(3,bw,zeros(size(bw)),zeros(size(bw))));
    c = bwconncomp(~bw,4); figure
    for i=1:size(c.PixelIdxList,2)
        if(find(idx == size(c.PixelIdxList{i},1)))
            tmp = RRGGBB(:,:,2);
            tmp(c.PixelIdxList{i}) = 1;
            RRGGBB(:,:,2) = tmp;
            imshow(RRGGBB);
        end
    end
end

