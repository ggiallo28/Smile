clear all; close all;
%x,y coordinates of circle
x=10;
y=10;
%radius r of circle
r=10;
%angle for each plot
ang = 0:0.01:2*pi;
%creating x,y array for the plot
xp=awgn(r*2*cos(ang),5); %add white noise with SNR 5:1
yp=awgn(r*sin(ang),5);

xp=xp+x;%displace circle by x,y amount
yp=yp+y;
%scatter plot
scatter(xp,yp);
zp = zeros(1,size(xp,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Custom Equation:
% a*x^2 + b*x*y + c*y^2 +d*x + e*y + f
% Se ometti c e d ottieni un'ellissi
% https://www.youtube.com/watch?v=VC_l5q4TTAA
%% Fit Ellipse
load('ellipse_data.mat');
l = 1;
i = find(xxxd == -1,l+1);
i = i(l:l+1);
idx = xxxd(i(1)+1:i(2)-1);
idy = yyyd(i(1)+1:i(2)-1);

% zz = zeros(size(idx));
% A = EllipseDirectFit([idx,idy]);
% x = linspace(0,1000000,1000);
% y = linspace(0,1000000,1000);
% EQ = A(1).*x.^2 + A(2).*x.*y + A(3).*y.^2 + A(4).*x + A(5).*y + A(6); 
% str = [num2str(A(1)),'x^2 + ',num2str(A(2)),'*x*y ',num2str(A(3)),'*y^2 + ',num2str(A(4)),'*x + ',num2str(A(5)),'*y ',num2str(A(6))];
% ezplot(str)

% Se l'asse maggiore è più piccolo della larghezza del cilindro allargalo
imshow(I)
[xx,idxx,IC] = unique(idx);
yy = idy(idxx);
%yy = smooth(xx,yy,0.2,'rlowess');
ellipse_t = fit_ellipse(idx,idy);
t=-pi:0.01:pi;
x=ellipse_t.X0_in+ellipse_t.a*cos(t+ellipse_t.phi);
y=ellipse_t.Y0_in+ellipse_t.b*sin(t);
hold on
plot(x,y)
hold on
plot(idx,idy);


el = imread('el.png');
el = el(:,:,1);
el = imclose(el,strel('disk',2));
imshow(el);
CCh = bwconncomp(el,8);
[idy,idx] = ind2sub(size(el),CCh.PixelIdxList{1});
ellipse_t = fit_ellipse(idx,idy);
t=-pi:0.01:pi;
x=ellipse_t.X0_in+ellipse_t.a*cos(t+10*ellipse_t.phi);
y=ellipse_t.Y0_in+ellipse_t.b*sin(t);
hold on
plot(x,y)




        