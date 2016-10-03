clc; clear all; close all;
%% Costants
lengthMirrors = 65.5;       %cm
mirror2Pivot = 3;           %cm
offset2Mirror = 3.2;
pivot2pivot = 28.5;         %cm
lengthHead = 20;            %cm % vertical 2*radius
widthHead = 18;             %cm % horizontal 2*radius
%% Parameters
angleRightMirror = 73.7;    %deg
angleLeftMirror = 106.3;    %deg
distanceCamera = 200;       %cm
fovHCamera = 83;            %deg
headPosx = 0.5*pivot2pivot; %cm % x0,y0 ellipse centre coordinates
headPosy = 16;              %cm
%% Logic
lengthPivotMirrors = lengthMirrors+mirror2Pivot;
angleBetweenMirror = abs(angleRightMirror-angleLeftMirror);
lengthHypotenuse = sqrt(offset2Mirror^2 + mirror2Pivot^2);
angle2Pivot = rad2deg(asin(offset2Mirror/lengthHypotenuse));  %deg angolo in alto del triangolo sotto.
N = 360/angleBetweenMirror;
% Proietto l'ipotenusa che si forma sull'asse X e Y in modo da avere le coordinate del punto in basso a destra. Il sistema
% di riferimento è posto sul perno sinisto.
%   ---------------->                          -------------------|->
%   |\                                                           /|
%   | \                                                         / |   %%  Standardizzazione dei riferimenti: centro sul perno sinistro
%   |__\ <- Punto                                    Punto ->  /__|
%       |                                                     |
% Punti interni specchi
inPointLeft(1) =  lengthHypotenuse*cos(deg2rad(angleLeftMirror-angle2Pivot));                   %X
inPointLeft(2) =  lengthHypotenuse*sin(deg2rad(angleLeftMirror-angle2Pivot));                   %Y
inPointRight(1) = pivot2pivot + lengthHypotenuse*cos(deg2rad(angleRightMirror+angle2Pivot));    %X, considero la traslazione lungo X per normalizzare i riferimenti
inPointRight(2) = lengthHypotenuse*sin(deg2rad(angleRightMirror+angle2Pivot));                  %Y
distanceBetweenMirror(1) = pdist([inPointLeft;inPointRight],'euclidean'); %shortest
% Punti esterni specchi
exPointLeft(1) = inPointLeft(1) + lengthMirrors*cos(deg2rad(angleLeftMirror));                  %X
exPointLeft(2) = inPointLeft(2) + lengthMirrors*sin(deg2rad(angleLeftMirror));                  %Y
exPointRight(1) = inPointRight(1) + lengthMirrors*cos(deg2rad(angleRightMirror));               %X,considero della traslazione lungo X
exPointRight(2) = inPointRight(2) + lengthMirrors*sin(deg2rad(angleRightMirror));               %Y
distanceBetweenMirror(2) = pdist([exPointLeft;exPointRight],'euclidean');
mLeft = ((-inPointLeft(2)) - (-exPointLeft(2)))/(inPointLeft(1) - exPointLeft(1));
kLeft = (inPointLeft(1)*(-exPointLeft(2)) - exPointLeft(1)*(-inPointLeft(2)))/(inPointLeft(1) - exPointLeft(1));
mRight = ((-inPointRight(2)) - (-exPointRight(2)))/(inPointRight(1) - exPointRight(1));
kRight = (inPointRight(1)*(-exPointRight(2)) - exPointRight(1)*(-inPointRight(2)))/(inPointRight(1) - exPointRight(1));
[mirrorsCenter(1),mirrorsCenter(2)] = inters2rette(mLeft,kLeft,mRight,kRight);
cameraCenter(1) = 0.5*pivot2pivot; cameraCenter(2) = -distanceCamera; 
angleHead = rad2deg(atan2((-headPosy-(-mirrorsCenter(2))),(headPosx-mirrorsCenter(1))));

% Inverti le Y. plot([x1 x2], [y1 y2])
figure,hold on,
plot([0 pivot2pivot], [0 0],'b');
plot([0 inPointLeft(1)], [0 -inPointLeft(2)],'b');
plot([pivot2pivot inPointRight(1)], [0 -inPointRight(2)],'b');
plot([inPointLeft(1) exPointLeft(1)], [-inPointLeft(2) -exPointLeft(2)],'b');
plot([inPointRight(1) exPointRight(1)], [-inPointLeft(2) -exPointRight(2)],'b');
scatter(mirrorsCenter(1),mirrorsCenter(2));
scatter(cameraCenter(1),cameraCenter(2));
%% Plot Heads
el_x=headPosx+0.5*widthHead*cos(-pi:0.01:pi);
el_y=-headPosy+0.5*lengthHead*sin(-pi:0.01:pi);
plot(el_x,el_y,'b')
h2m(1) = abs(angleLeftMirror-angleHead); % Angle 2 Left
h2m(2) = abs(angleRightMirror-angleHead); % Angle 2 Right
k = 1; v = zeros(1,round(N)-1);
for i=1:round(N)-1
    v(i) = 2*h2m(k);
    if(k == 1)
        k=2;
    else
        k=1;
    end
end









