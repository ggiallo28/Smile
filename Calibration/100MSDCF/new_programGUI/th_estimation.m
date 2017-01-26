function Container = th_estimation(Container, bw)
%     bw_split = bwareaopen(imdilate(bw,strel('square',20)),5*Container.Threshold);
%     split = medfilt2(sum(bw_split),[1,10]);
%     split = find(split~=0);
%     [x, y] = compute_projections(bw, split);
%     assert(size(x,1) == 3,'Errore th_estimation, non riesco a dividere');
%     
    bw_close = imclose(bw,strel('square',5));
    idx = sort(cell2mat(struct2cell(regionprops(bwconncomp(~bw_close,4),'Area'))),'descend')';
    idx = idx(2:end); % Rimuovo la prima, enorme, � sicuramente background
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
        if ~Container.isGUI
           plot(fitresult,'r-',idy,idx,'k.',outliers,'m*');
%         else
%            set(0,'DefaultFigureVisible','off');
%            plot(fitresult,'r-',idy,idx,'k.',outliers,'m*');
%            [Printed,~]=export_fig('-native');
%            close all;
%            set(0,'DefaultFigureVisible','on');
%            imshow(Printed,'Parent',Container.app.UIAxes6);
%            drawnow;
        end
        idy(outliers) = [] ; idx(outliers) = [];
%        pause(0.01);
    end
    step_dim_a = floor(size(idx)/3);
    mean_square1 = mean(idx);
    size_a = size(idx,1);
    mean_big1 = mean(fitresult(idy(1:step_dim_a)));
    mean_medium1 = mean(fitresult(idy(step_dim_a+1:2*step_dim_a)));
    mean_small1 = mean(fitresult(idy(2*step_dim_a+1:3*step_dim_a)));
%     if Container.isGUI
%         gui_axes = Container.app.UIAxes8;
%         RRGGBB = double(cat(3,bw,zeros(size(bw)),zeros(size(bw))));
%         c = bwconncomp(~bw_close,4);
%         figure
%         for i=1:size(c.PixelIdxList,2)
%             if(find(idx == size(c.PixelIdxList{i},1)))
%                 tmp = RRGGBB(:,:,2);
%                 tmp(c.PixelIdxList{i}) = 1;
%                 RRGGBB(:,:,2) = tmp;
%                 imshow(RRGGBB,'Parent',gui_axes);
%                 drawnow;
%             end
%         end
%     end
    % Provare a trattare separatamente i riflessi, si hanno problemi con 5
    % riflessi nella stima dei quadrati
    bw_open = imerode(bw,strel('square',10));
    bw_open = bwareaopen(bw_open,200);
    idx = sort(cell2mat(struct2cell(regionprops(bwconncomp(bw_open,4),'Area'))),'descend' )';
    remove = round(size(idx,1)*0.30);
    idx = idx(1:end-remove);
    idx(idx<mean(idx)-0.8*std(idx) & idx<mean(idx)+0.8*std(idx)) = [];
    idy = 1:size(idx,1);
    idy = idy';
    fitresult = createLine(idy,idx);
    outliers = []; toCont = 1;
    opts = fitoptions( 'Method', 'LinearLeastSquares' , 'Exclude',outliers); ft = fittype( 'poly1' );
    while toCont >0
        fdata = feval(fitresult,idy);
        I = abs(fdata - idx) > 1.3*std(idx);
        toCont = size(find(I>0),1);
        outliers = excludedata(idy,idx,'indices',I);
        [fitresult, ~] = fit( idy, idx, ft, opts );
        if ~Container.isGUI
            plot(fitresult,'r-',idy,idx,'k.',outliers,'m*');
        end
        id_tmp = idy;
        id_tmp(outliers) = [] ;
        if(size(id_tmp)<=2)
            break;
        end
        idy(outliers) = [] ; idx(outliers) = [];
    end
    step_dim_b = floor(size(idy)/3);
    mean_square2 = mean(idx);
    size_b = size(idx,1);
    mean_big2 = mean(fitresult(idy(1:step_dim_b)));
    mean_medium2 = mean(fitresult(idy(step_dim_b+1:2*step_dim_b)));
    mean_small2 = mean(fitresult(idy(2*step_dim_b+1:3*step_dim_b)));
%     if Container.isGUI
%         gui_axes = Container.app.UIAxes8;
%         RRGGBB = double(cat(3,bw_close,zeros(size(bw)),zeros(size(bw))));
%         c = bwconncomp(bw_open,4);
%         for i=1:size(c.PixelIdxList,2)
%             if(find(idx == size(c.PixelIdxList{i},1)))
%                 tmp = RRGGBB(:,:,2);
%                 tmp(c.PixelIdxList{i}) = 1;
%                 RRGGBB(:,:,2) = tmp;
%                 imshow(RRGGBB,'Parent',gui_axes);
%                 drawnow;
%             end
%         end
%     end
%% mean  
    ratio_a = size_a/(size_a+size_b);
    ratio_b = size_b/(size_a+size_b);
    step_a = step_dim_a/(step_dim_b+step_dim_a);
    step_b = step_dim_b/(step_dim_b+step_dim_a);
    if isnan(mean_square2)
        mean_square = mean_square1;
    elseif isnan(mean_square1)
        mean_square = mean_square2;
    else
        mean_square = 2*mean([ratio_b*mean_square2, ratio_a*mean_square1]);
    end
    Container.size_square = mean_square;
    Container.Threshold = round(mean_square/Container.fraction);
%% big   
    if isnan(mean_big1)
        Container.size_big = mean_big2;
    elseif isnan(mean_big2)
        Container.size_big = mean_big1;
    else
        Container.size_big = 2*mean([step_b*mean_big1, step_a*mean_big2]);
    end
%% mid
    if isnan(mean_medium1)
        Container.size_medium = mean_medium2;
    elseif isnan(mean_medium2)
        Container.size_medium = mean_medium1;
    else
       Container.size_medium = 2*mean([step_a*mean_medium2, step_b*mean_medium1]);
    end
%% small
    if isnan(mean_small1)
        Container.size_big = mean_small2;
    elseif isnan(mean_small2)
        Container.size_big = mean_small1;
    else
         Container.size_small = 2*mean([step_b*mean_small1, step_a*mean_small2]);
    end 
end