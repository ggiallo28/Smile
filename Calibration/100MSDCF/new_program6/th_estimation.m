function Container = th_estimation(Container, bw)
    bw_close = imclose(bw,strel('square',5));
    idx = sort(cell2mat(struct2cell(regionprops(bwconncomp(~bw_close,4),'Area'))),'descend')';
    remove = round(size(idx,1)*0.30);
    idx = idx(2:end-remove); % Rimuovo la prima, enorme, è sicuramente background
    idy = 1:size(idx,1);
    idy = idy';
    fitresult = createLine(idy,idx);
    outliers = []; toCont = 1;
    opts = fitoptions( 'Method', 'LinearLeastSquares' , 'Exclude',outliers); ft = fittype( 'poly1' );
    while toCont >0
        fdata = feval(fitresult,idy);
        I = abs(fdata - idx) > 1.7*std(idx);
        toCont = size(find(I>0),1);
        outliers = excludedata(idy,idx,'indices',I);
        [fitresult, ~] = fit( idy, idx, ft, opts );
        plot(fitresult,'r-',idy,idx,'k.',outliers,'m*');
        idy(outliers) = [] ; idx(outliers) = [];
        pause(0.01);
    end
    step_dim = floor(size(idx)/3);
    mean_square1 = mean(idx);
    mean_big1 = mean(fitresult(idy(1:step_dim)));
    mean_medium1 = mean(fitresult(idy(step_dim+1:2*step_dim)));
    mean_small1 = mean(fitresult(idy(2*step_dim+1:3*step_dim)));
%     RRGGBB = double(cat(3,bw,zeros(size(bw)),zeros(size(bw))));
%     c = bwconncomp(~bw_close,4);
%     for i=1:size(c.PixelIdxList,2)
%         if(find(idx == size(c.PixelIdxList{i},1)))
%             tmp = RRGGBB(:,:,2);
%             tmp(c.PixelIdxList{i}) = 1;
%             RRGGBB(:,:,2) = tmp;
%             imshow(RRGGBB);
%         end
%     end
    bw_open = imerode(bw,strel('square',10));
    bw_open = bwareaopen(bw_open,200);
    idx = sort(cell2mat(struct2cell(regionprops(bwconncomp(bw_open,4),'Area'))),'descend' )';
    remove = round(size(idx,1)*0.30);
    idx = idx(1:end-remove);
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
    step_dim = floor(size(idx)/3);
    mean_square2 = mean(idx);
    mean_big2 = mean(fitresult(idy(1:step_dim)));
    mean_medium2 = mean(fitresult(idy(step_dim+1:2*step_dim)));
    mean_small2 = mean(fitresult(idy(2*step_dim+1:3*step_dim)));
%     RRGGBB = double(cat(3,bw_close,zeros(size(bw)),zeros(size(bw))));
%     c = bwconncomp(bw_open,4);
%     for i=1:size(c.PixelIdxList,2)
%         if(find(idx == size(c.PixelIdxList{i},1)))
%             tmp = RRGGBB(:,:,2);
%             tmp(c.PixelIdxList{i}) = 1;
%             RRGGBB(:,:,2) = tmp;
%             imshow(RRGGBB);
%         end
%     end
    mean_square = mean([mean_square2, mean_square1]);
    Container.size_square = mean_square;
    Container.size_big = mean([mean_big1, mean_big2]);
    Container.size_medium = mean([mean_medium2, mean_medium1]);
    Container.size_small = mean([mean_small1, mean_small2]);
    Container.Threshold = round(mean_square/Container.fraction);
end