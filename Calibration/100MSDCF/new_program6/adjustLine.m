function [fitresult, isDone] = adjustLine(fitresult,idy,idx)
    % Iterative fitting
    outliers = []; toCont = 1;
    opts = fitoptions( 'Method', 'LinearLeastSquares' , 'Exclude',outliers); ft = fittype( 'poly1' );
    isDone = true;
    while toCont >0
        fdata = feval(fitresult,idy);
        I = abs(fdata - idx) > 1.5*std(idx);
        toCont = size(find(I>0),1);
        outliers = excludedata(idy,idx,'indices',I);
        [fitresult, ~] = fit( idy, idx, ft, opts );
%        imshow(color_square_edge_h); hold on; plot(fitresult_top,'r-',idy_top,idx_top,'k.',outliers,'m*');
        idy(outliers) = [] ; idx(outliers) = [];
        if(isempty(idx) || isempty(idy))
            isDone = false;
            return 
        end
    end
end

