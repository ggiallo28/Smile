%%Computes the clusters of Pixels based upon their color.
%% The algorithm is based uopon binary tree quantization technique
%% described by Orchard and Bouman.

%% input arguments
%%Image = input image should be RGB.
%% K= Maximum possible clusters required.

%% output arguments
%% C= A cell array containing cluster structures. Each cluster has
%% following fields:
%% point = An integer array containing linear indices of points in
%% clusters.
%% M = Sum of RGB values of all the pixels in cluster
%% N=  Number of Pixels in the cluster
%% R=  The sum of squares of pixels in the cluster
%% vector= Principal eigen vector of cluster

%% References:
%% ORCHARD, M. T., AND BOUMAN, C. A. 1991. Color Quantiza-tion of Images. 
%%IEEE Transactions on Signal Processing 39, 12,2677–2690.



function C=binary_quant(Image,K)
%% Read the image and convert it into double.
I=imread(Image);
I=im2double(I);
si=1;
R_Comp=I(:,:,1);
G_Comp=I(:,:,2);
B_Comp=I(:,:,3);
%% Treat complete image as single cluster and compute different fields of
%% cluster.
R1=zeros(3,3);
R1(1,1)=sum(sum(R_Comp.*R_Comp));
R1(1,2)=sum(sum(R_Comp.*G_Comp));
R1(1,3)=sum(sum(R_Comp.*B_Comp));
R1(2,1)=sum(sum(G_Comp.*R_Comp));
R1(2,2)=sum(sum(G_Comp.*G_Comp));
R1(2,3)=sum(sum(G_Comp.*B_Comp));
R1(3,1)=sum(sum(B_Comp.*R_Comp));
R1(3,2)=sum(sum(B_Comp.*G_Comp));
R1(3,3)=sum(sum(B_Comp.*B_Comp));
R_Comp=R_Comp(:);
G_Comp=G_Comp(:);
B_Comp=B_Comp(:);
N1=size(I(:),1)/3;
M=sum(sum(I));
M1=zeros(1,3);
M1(1,1)=M(:,:,1);
M1(1,2)=M(:,:,2);
M1(1,3)=M(:,:,3);
%% Compute covariance of the Image.
R1_bar=R1-((M1*M1')/N1);
%% Compute eigen values and eigen vectors of the image.
[V,D]=eig(R1_bar);
point=1:1:size(I,1)*size(I,2);
%% Form the initial cluster.
cluster.R=R1;
cluster.M=M1;
cluster.N=N1;
cluster.point=point;
cluster.vector=V(:,1);
%% Place the cluster in Cell array.
C{si}=cluster;
%% Store eigen value in a seperate array.
val(si)=abs(D(1,1));
for i=2:1:K
    if(max(val(:))==-1)
        break;
    end
    %% Pick the cluster having highest eigen value and split it into two parts.
  index=find(max(val(:))==val);
  vector=C{index}.vector;
  point=C{index}.point;
  point=point(:);
  M=C{index}.M;
  N=C{index}.N;
  Q=M/N;
  R=C{index}.R;
  R1=zeros(3,3);
  M1=zeros(1,3);
         a= vector(1,1)*R_Comp(point(:)) + vector(2,1)*G_Comp(point(:)) + vector(3,1)*B_Comp(point(:));
         %% Pixels of the cluster are divided into two parts based on their
         %% closeness to a plane perpendicular to principal eigen vector
         %% and passing through mean.
         point1=point(find(a<=vector'*Q'));
         point2=point(find(a>vector'*Q'));
         %% Compute various statistical value for the newly formed
         %% clusters.
              R1(1,1)=sum(R_Comp(point1(:)).*R_Comp(point1(:)));
              R1(1,2)=sum(R_Comp(point1(:)).*G_Comp(point1(:)));
              R1(1,3)=sum(R_Comp(point1(:)).*B_Comp(point1(:)));
              R1(2,1)=sum(G_Comp(point1(:)).*R_Comp(point1(:)));
              R1(2,2)=sum(G_Comp(point1(:)).*G_Comp(point1(:)));
              R1(2,3)=sum(G_Comp(point1(:)).*B_Comp(point1(:)));
              R1(3,1)=sum(B_Comp(point1(:)).*R_Comp(point1(:)));
              R1(3,2)=sum(B_Comp(point1(:)).*G_Comp(point1(:)));
              R1(3,3)=sum(B_Comp(point1(:)).*B_Comp(point1(:)));
              M1(1,1)=sum(R_Comp(point1(:)));
              M1(1,2)=sum(G_Comp(point1(:)));
              M1(1,3)=sum(B_Comp(point1(:)));
  R2=R-R1;
  M2=M-M1;
  N1=size(point1,1);
  N2=size(point2,1);
  %% Place the clusters into cell array and delete the parent cluster.
  if(N1>0 && N2>0)
  cluster1.point=point1;
  cluster1.R=R1;
  cluster1.M=M1;
  cluster1.N=N1;
  R1_bar=R1-((M1*M1')/N1);
  [V1 D1]=eig(R1_bar);
  cluster1.vector=V1(:,1);
  cluster2.R=R2;
  cluster2.N=N2;
  cluster2.M=M2;
  cluster2.point=point2;
  R2_bar=R2-((M2*M2')/N2);
  [V2 D2]=eig(R2_bar);
  cluster2.vector=V2(:,1);
  val(si+1)=abs(D2(1,1));
  val(index)=abs(D1(1,1));
  C{index}=cluster1;
  C{si+1}=cluster2;
  si=si+1;
  else
     val(index)=-1;
     i=i-1;
  end
end
%% display original Image.
figure,imshow(I);
%% display newly formed clusters.
for i=1:1:si;
    temp=C{i}.point;
    I1=zeros(size(I,1),size(I,2),3);
       [y_p x_p]=ind2sub([size(I,1) size(I,2)],temp);
  for t=1:1:size(y_p,1)
           y=y_p(t);
           x=x_p(t);
        I1(y,x,1)=I(y,x,1);
        I1(y,x,2)=I(y,x,2);
        I1(y,x,3)=I(y,x,3);
    end
    figure,imshow(I1);
end
end 