function board = growCheckerboard(points, scores, Ix2, Iy2, Ixy, theta)

% Exit immediately if no corner points were found
if isempty(scores)
    if isempty(coder.target)
        board = struct('BoardIdx', zeros(3), 'BoardCoords', zeros(3,3,3), ...
            'Energy', Inf, 'isValid', 0);
    else
        board = vision.internal.calibration.checkerboard.Checkerboard;
    end
    return;
end



% only use corners with high scores as seeds to reduce computation
seedIdx = 1:size(points, 1);
[~, sortedIdx] = sort(scores(seedIdx), 'descend');
seedIdx = seedIdx(sortedIdx);
seedIdx = seedIdx(1:round(numel(seedIdx / 2)));

angleThreshold = 3 * pi / 16;

if isempty(coder.target)
    v1_matrix = [];
    v2_matrix = [];
    seedIdx_matrix = [];
    
    for i = seedIdx
        [v1, v2] = cornerOrientations(Ix2, Iy2, Ixy, round(points(i, :)));
        alpha1 = abs(atan2(v1(2), v1(1)));
        alpha2 = abs(atan2(v2(2), v2(1)));
        if abs(abs(alpha1 - pi) - theta) > angleThreshold && ...
                abs(abs(alpha2 - pi) - theta) > angleThreshold
            continue;
        else
            v1_matrix = [v1_matrix;v1]; %#ok<AGROW>
            v2_matrix = [v2_matrix;v2]; %#ok<AGROW>
            seedIdx_matrix = [seedIdx_matrix;i]; %#ok<AGROW>
        end
    end
    
    board = visionInitializeAndExpandCheckerboard(seedIdx_matrix,single(points),v1_matrix,v2_matrix);
else
    previousBoard = vision.internal.calibration.checkerboard.Checkerboard;
    currentBoard = vision.internal.calibration.checkerboard.Checkerboard;
    for i = 1:numel(seedIdx)
        [v1, v2] = cornerOrientations(Ix2, Iy2, Ixy, round(points(seedIdx(i), :)));
        alpha1 = abs(atan2(v1(2), v1(1)));
        alpha2 = abs(atan2(v2(2), v2(1)));
        if abs(abs(alpha1 - pi) - theta) > angleThreshold && ...
                abs(abs(alpha2 - pi) - theta) > angleThreshold
            continue;
        end
        
        currentBoard.initialize(seedIdx(i), points, v1, v2);
        expandBoardFully(currentBoard);
        if currentBoard.Energy < previousBoard.Energy            
            tmpBoard = previousBoard;
            previousBoard = currentBoard;
            currentBoard = tmpBoard;
        end
    end
    board = previousBoard;
end
end