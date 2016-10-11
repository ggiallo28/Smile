clear all, clc, close all;
% add all needed function paths
addpath .\coherenceFilter
addpath .\GLtree3DMex

load clown;
X = imread('_DSC0072.JPG');
X = imcrop(X);
I_gray = X;
%smooth the image by coherence filter:
filted_I = CoherenceFilter(I_gray,struct('T',5,'rho',2,'Scheme','I', 'sigma', 1));
%adjacent neighborhood  model:
L = graphSeg(filted_I, 0.5, 50, 2, 0);
%k-nearest neighborhood model:
Lnn = graphSeg(filted_I, 0.5/sqrt(3), 50, 10, 1);
%display:
subplot(3, 1, 1), imshow(I_gray, []), title('original image');
subplot(3, 1, 2), imshow(label2rgb(L)), title('adjacent neighborhood based segmentation');
subplot(3, 1, 3), imshow(label2rgb(Lnn)), title('k nearest neighborhood based segmentation');