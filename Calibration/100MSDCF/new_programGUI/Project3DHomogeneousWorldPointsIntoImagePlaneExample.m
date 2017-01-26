%% Compute Camera Matrix

% Copyright 2015 The MathWorks, Inc.

% Create a set of calibration images.
images = imageSet(fullfile(toolboxdir('vision'), 'visiondata', ...
  'calibration', 'slr'));

% Detect the checkerboard corners in the images.
[imagePoints, boardSize] = detectCheckerboardPoints(images.ImageLocation);

% Generate the world coordinates of the checkerboard corners in the
% pattern-centric coordinate system, with the upper-left corner at (0,0).
squareSize = 29; % in millimeters
worldPoints = generateCheckerboardPoints(boardSize, squareSize);

% Calibrate the camera.
cameraParams = estimateCameraParameters(imagePoints, worldPoints);

% Load image at new location.
imOrig = imread(fullfile(matlabroot, 'toolbox', 'vision', 'visiondata', ...
    'calibration', 'slr', 'image9.jpg'));
figure; imshow(imOrig);
title('Input Image');

% Undistort image.
im = undistortImage(imOrig, cameraParams);

% Find reference object in new image.
[imagePoints, boardSize] = detectCheckerboardPoints(im);

% Compute new extrinsics.
[rotationMatrix, translationVector] = extrinsics(...
  imagePoints, worldPoints, cameraParams);

% Calculate camera matrix
P = cameraMatrix(cameraParams, rotationMatrix, translationVector)

