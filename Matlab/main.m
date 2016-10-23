clc; clear all; close all;
%% Costants
lengthMirrors = 65.5;       %cm
mirror2Pivot = 3;           %cm Ci sta una struttura ad angolo retto che collega lo specchio all'ase di rotazione, abbiamo lungo l'asse dello specchio uno sfasamento di 3 cm
offset2Mirror = 3.2;        %cm ed uno sfasamento di 3.2 cm in avanti.
pivot2pivot = 28.5;         %cm Distanza tra i due perni
lengthHead = 20;            %cm % vertical 2*radius
widthHead = 10;             %cm % horizontal 2*radius
%% Parameters
angleRightMirror = 74.64;    %deg
angleLeftMirror = 130.81;    %deg
distanceCamera = 200;       %cm
fovHCamera = 83;            %deg
headPosx = 0.5*pivot2pivot; %cm % x0,y0 ellipse centre coordinates
headPosy = 36;              %cm
%% Logic
lengthPivotMirrors = lengthMirrors+mirror2Pivot; % Lunghezza dello specchio più offset dovuto al fatto che lo specchio non finisce dove è posto l'asse di rotazione (pivot).
angleBetweenMirror = abs(angleRightMirror-angleLeftMirror); % Angolo tra gli specchi
lengthHypotenuse = sqrt(offset2Mirror^2 + mirror2Pivot^2); % Ipotenusa dello sfasamento
angle2Pivot = rad2deg(asin(offset2Mirror/lengthHypotenuse));  %deg angolo in alto del triangolo sotto. Angolo fisso dovuto alla struttura (entità dell'offset)
N = 360/angleBetweenMirror; % FORMULA 1: Dice quanti riflessi vengono prodotti, quelli visibili saranno un sottoinsieme. Da usare come verifica.
% Proietto l'ipotenusa che si forma sull'asse X e Y in modo da avere le coordinate del punto in basso a destra. Il sistema
% di riferimento è posto sul perno sinisto.
%   ---------------->                          -------------------|->
%   |\                                                           /|
%   | \                                                         / |   %%  Standardizzazione dei riferimenti: centro sul perno sinistro
%   |__\ <- Punto                                    Punto ->  /__|
%       |                                                     |
% Calcolo posizione punti interni specchi (estremi inferiori)
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
% Coefficienti angolari e costante delle rette passanti per i due specchi, mi serve per calcolare la posizione del centro di rotazione
% y = (x*y1 - x*y2 + x1*y2 - x2*y1)/(x1 - x2);
% y = x*(y1 - y2)/(x1 - x2) + (x1*y2 - x2*y1)/(x1 - x2);
% m = (y1 - y2)/(x1 - x2)
% k = (x1*y2 - x2*y1)/(x1 - x2)
mLeft = ((-inPointLeft(2)) - (-exPointLeft(2)))/(inPointLeft(1) - exPointLeft(1));
kLeft = (inPointLeft(1)*(-exPointLeft(2)) - exPointLeft(1)*(-inPointLeft(2)))/(inPointLeft(1) - exPointLeft(1));
mRight = ((-inPointRight(2)) - (-exPointRight(2)))/(inPointRight(1) - exPointRight(1));
kRight = (inPointRight(1)*(-exPointRight(2)) - exPointRight(1)*(-inPointRight(2)))/(inPointRight(1) - exPointRight(1));
[mirrorsCenter(1),mirrorsCenter(2)] = inters2rette(mLeft,kLeft,mRight,kRight);
% Sto ipotizzando la camera al centro (X) e posta ad una distanza (Y) definita precedentemente
cameraCenter(1) = 0.5*pivot2pivot; cameraCenter(2) = -distanceCamera; 
% Angolo della retta passante per il centro degli specchi e il centro della testa, serve per capire di quanto ruotare
% Deommentare per allineare la testa e la camera con il centro di rotazione degli specchi
% headPosx = mirrorsCenter(1);
% cameraCenter(1) = mirrorsCenter(1);
angleHead = rad2deg(atan2((-headPosy-mirrorsCenter(2)),(headPosx-mirrorsCenter(1))));
if angleHead<0
    angleHead = -angleHead;
end

% [el_x, el_y] = genHead(headPosx, -headPosy, widthHead, lengthHead);
% v = genMirroring(angleLeftMirror, angleRightMirror, angleHead, N);
% el_X = zeros(size(v,2),size(el_x,2));
% el_Y = zeros(size(v,2),size(el_y,2));
% distOR = sqrt((mirrorsCenter(1)-mean(el_x))^2 + (mirrorsCenter(2)-mean(el_y))^2);
% %% Flip per effetto mirroring
% for i=1:size(v,2);
%     el_x_flip = el_x;
%     el_y_flip = el_y;
%     if v(2,i) ~= 0
%         [el_x_flip,el_y_flip] = rotate(el_x_flip,el_y_flip,headPosx,-headPosy,180,'y');
%     end
%     [el_X(i,:),el_Y(i,:)] = rotate(el_x_flip,el_y_flip,mirrorsCenter(1),mirrorsCenter(2),v(1,i),'z');
%     new_center_x = mean(el_X(i,:));
%     new_center_y = mean(el_Y(i,:));
%     distROT = sqrt((mirrorsCenter(1)-mean(el_X(i,:)))^2 + (mirrorsCenter(2)-mean(el_Y(i,:)))^2);
%     assert(abs(distOR-distROT)<exp(-10));
% end

[el_x, el_y] = genHead(headPosx, -headPosy, widthHead, lengthHead);
el_X = zeros(floor(N)-1,size(el_x,2));
el_Y = zeros(floor(N)-1,size(el_y,2));
for i=1:size(el_x,2)
    anglePoint = rad2deg(atan2((el_y(i)-mirrorsCenter(2)),(el_x(i)-mirrorsCenter(1))));
    if anglePoint<0
        anglePoint = -anglePoint;
    end
    v = genMirroring(angleLeftMirror, angleRightMirror, anglePoint, N);
    for j=1:floor(N)-1
       [el_X(j,i),el_Y(j,i)] = rotate(el_x(i),el_y(i),mirrorsCenter(1),mirrorsCenter(2),v(1,j),'z'); 
    end 
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
% y = (x*y1 - x*y2 + x1*y2 - x2*y1)/(x1 - x2);
% y = x*(y1 - y2)/(x1 - x2) + (x1*y2 - x2*y1)/(x1 - x2);
y1 = -inPointLeft(2); y2 = -inPointRight(2);
x1 = inPointLeft(1); x2 = inPointRight(1);
mCoeff = (y1 - y2)/(x1 - x2);
kCoeff = (x1*y2 - x2*y1)/(x1 - x2);
xRange = x1:0.1:x2;
yRange = xRange*mCoeff + kCoeff;

XRange = zeros(floor(N)-1,size(xRange,2));
YRange = zeros(floor(N)-1,size(yRange,2));
for i=1:size(xRange,2)
    anglePoint = rad2deg(atan2((yRange(i)-mirrorsCenter(2)),(xRange(i)-mirrorsCenter(1))));
    if anglePoint<0
        anglePoint = -anglePoint;
    end
    v = genMirroring(angleLeftMirror, angleRightMirror, anglePoint, N);
    for j=1:floor(N)-1
       [XRange(j,i),YRange(j,i)] = rotate(xRange(i),yRange(i),mirrorsCenter(1),mirrorsCenter(2),v(1,j),'z'); 
    end 
end

for j=1:floor(N)-1
    plot(XRange(j,:),YRange(j,:))
end
%% TODO calcolare indice di copertura





