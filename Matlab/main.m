clc; clear all; close all;
%% Costants
lengthMirrors = 65.5;       %cm
mirror2Pivot = 3;           %cm
offset2Mirror = 3.2;
pivot2pivot = 28.5;         %cm
lengthHead = 20;            %cm % vertical 2*radius
widthHead = 18;             %cm % horizontal 2*radius
%% Parameters
angleRightMirror = 43.7;    %deg
angleLeftMirror = 106.3;    %deg
distanceCamera = 200;       %cm
fovHCamera = 83;            %deg
headPosx = 6+0.5*pivot2pivot; %cm % x0,y0 ellipse centre coordinates
headPosy = 36;              %cm
%% Logic
lengthPivotMirrors = lengthMirrors+mirror2Pivot;
angleBetweenMirror = abs(angleRightMirror-angleLeftMirror);
lengthHypotenuse = sqrt(offset2Mirror^2 + mirror2Pivot^2);
angle2Pivot = rad2deg(asin(offset2Mirror/lengthHypotenuse));  %deg angolo in alto del triangolo sotto.
N = 360/angleBetweenMirror;
% Proietto l'ipotenusa che si forma sull'asse X e Y in modo da avere le coordinate del punto in basso a destra. Il sistema
% di riferimento � posto sul perno sinisto.
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
angleHead = rad2deg(atan2((-headPosy-mirrorsCenter(2)),(headPosx-mirrorsCenter(1))));
if angleHead<0
    angleHead = -angleHead;
end
[el_x, el_y] = genHead(headPosx, -headPosy, widthHead, lengthHead);
v = genMirroring(angleLeftMirror, angleRightMirror, angleHead, N);
el_X = zeros(size(v,2),size(el_x,2));
el_Y = zeros(size(v,2),size(el_y,2));
for i=1:size(v,2);
    el_x_flip = el_x;
%     if v(3,i) == -1
%         el_x_flip = fliplr(el_x_flip); %% Verificare
%     elseif v(3,i) == 1
%         el_x_flip = fliprl(el_x_flip);
%     end
    [el_X(i,:),el_Y(i,:)] = rotate(el_x_flip,el_y,mirrorsCenter(1),mirrorsCenter(2),v(1,i));
end

%% Plot Inverti le Y. plot([x1 x2], [y1 y2])
figure,hold on,axis([-distanceCamera distanceCamera -distanceCamera distanceCamera])
plot([0 pivot2pivot], [0 0],'b');
plot([0 inPointLeft(1)], [0 -inPointLeft(2)],'b');
plot([pivot2pivot inPointRight(1)], [0 -inPointRight(2)],'b');
plot([inPointRight(1) exPointRight(1)], [-inPointRight(2) -exPointRight(2)],'b');
plot([inPointLeft(1) inPointRight(1)], [-inPointLeft(2) -inPointRight(2)],'b');
plot([inPointLeft(1) exPointLeft(1)], [-inPointLeft(2) -exPointLeft(2)],'b');

scatter(mirrorsCenter(1),mirrorsCenter(2));
scatter(cameraCenter(1),cameraCenter(2));
plot(el_x,el_y,'b')
plot([min(el_x) cameraCenter(1)], [median(el_y) cameraCenter(2)],'b');
plot([max(el_x) cameraCenter(1)], [median(el_y) cameraCenter(2)],'b');

for i=1:size(v,2);
    plot(el_X(i,:),el_Y(i,:));
    [p1, p2] = getTangentLine(el_X(i,:),el_Y(i,:), cameraCenter);
    plot([p1(1) cameraCenter(1)], [p1(2) cameraCenter(2)]);
    plot([p2(1) cameraCenter(1)], [p2(2) cameraCenter(2)]);
end

%% Disegnare elemento centrale
mul = 2;
for i=1:size(v,2);
    if v(2,i) == 1
        [xx,yy] = rotate(inPointLeft(1), -inPointLeft(2),mirrorsCenter(1),mirrorsCenter(2),mul*angleBetweenMirror);
        plot([inPointRight(1) xx],[-inPointRight(2) yy]);
    else
        [xx,yy] = rotate(inPointRight(1), -inPointRight(2),mirrorsCenter(1),mirrorsCenter(2),-mul*angleBetweenMirror);
        plot([ inPointLeft(1) xx],[-inPointLeft(2) yy]);
    end
    if mod(i,2) == 0 && i> 1
        mul = mul+2;
    end
end
%% TODO fare flip e disegnare figura al centro
%% TODO calcolare indice di copertura





