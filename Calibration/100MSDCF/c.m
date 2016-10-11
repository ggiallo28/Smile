      % http://it.mathworks.com/matlabcentral/answers/260232-detecting-center-pixels-of-noisy-and-blurred-checkerboard-pattern      
      IM1 = imread('_DSC0077.JPG');
      IM1 = imcrop(IM1);
      %Removing the predicted rotation in the image
      IM1 = imresize(IM1, 8);
      IM1 = imrotate(IM1, 0.1559, 'bicubic');
      IM1 = imresize(IM1, 0.125);
      %Using spatial method to detect centers
      %Image enhancing for detection  
      IM_spat = double(adaptivethreshold(IM1,21,0.075,0));
      PSF = fspecial('gaussian',21,2);
      IM_spat = edgetaper(IM_spat,PSF);
      smooth_filter = fspecial('gaussian',21,1.5)
      IM_spat = imfilter(IM_spat,smooth_filter,'replicate');
      % Calculating peaks using FastPeakFind algorithm    
      peak_white=FastPeakFind(IM_spat);  
      IM_spat2 = imcomplement(IM_spat);
      peak_black=FastPeakFind(IM_spat2);
      PeakX=[peak_black(1:2:end);peak_white(1:2:end)];
      PeakY= [peak_black(2:2:end);peak_white(2:2:end)];
      % Delaunay Triangulation to further precise the centers         
      DT = delaunayTriangulation(PeakX, PeakY);
      DeX=DT.Points(:,1);
      DeY=DT.Points(:,2);
      figure, imshow(IM1); hold on
      plot(DeX,DeY,'r*');
      hold off;
      %Cropping the image a bit to remove corners extreme distortion
      IDX = find(DeX > 20 & DeX < 1250);
      DeX1 = DeX(IDX);
      DeY1 = DeY(IDX);
      IDY = find(DeY1 > 20 & DeY1 < 1000);
      DeX1 = DeX1(IDY);
      DeY1 = DeY1(IDY);
      %Squeezing the rows to remove repetition for meshgrid, because the centers 
      %detected are a bit up and down in the same square row
      uniX= unique(DeX1);
      countX=size(uniX);
      for i = 1:countX(1)
       I=find(DeX1==uniX(i));
       K = DeY1(I);
       countY=size(K);
        for J=1:countY(1)
          GenMat(i,J)=K(J);
        end   
      end
      j = 0;
      i = 1;
      GenMat1 = zeros(0,250);
      while(i <= (size(uniX,1)-1))
          buff = GenMat(i,:);
          I = find(buff ~= 0);
          buff1 = buff(I);
          iniI = i;
          while(uniX(i+1)-uniX(i)<=3)
              buff = GenMat(i+1,:);
              I = find(buff ~= 0);
              if size(buff1,2) > 200
                  break;
              end;
              buff1 = [buff1,buff(I)];
              i = i + 1;
              if i == (size(uniX,1))
                  break;
              end;
          end;    
          if size(buff1,2) < 200
              j = j + 1;
              nuniX(j) = uniX(iniI + round((i-iniI)/2));
              GenMat1(j,1:size(buff1,2)) = buff1;
          end;
          i = i + 1;
      end;
      SortGen1=sort(GenMat1,2);
      for i = 1:size(SortGen1,1)
          I1 = find(SortGen1(i,:) == 0);
          I2 = find(SortGen1(i,:) ~= 0);
          SortGen2 = unique(SortGen1(i,I2));
          I = zeros(1,250-size(SortGen2,2));
          SortGen1(i,:) = [SortGen2 I];
      end;
      Start = min(SortGen1(:,1));
      for i = 1:size(SortGen1,1)
          k = 1;
          myGen = 0;
          myGen(1,1) = SortGen1(i,1);
          for j = 2:size(SortGen1,2)
              k = k + 1;
              if ((SortGen1(i,j) - SortGen1(i,j-1)) > SQ_SIZE)
                  N = round((SortGen1(1,2) - SortGen1(1,2-1))/8) - 1;
                  for l = 1:N
                      myGen(1,k) = SortGen1(i,j-1) + N * SQ_SIZE;
                      k = k + 1;
                  end;
              end;
              myGen(i,k) = SortGen1(i,j);
          end;
          I = zeros(1,250-size(myGen,2));
          SortGen2(i,:) = [myGen I];
      end;
      SortGen1=sort(GenMat1,2);
      for i = 1:size(SortGen1,1)
          I1 = find(SortGen1(i,:) == 0);
          I2 = find(SortGen1(i,:) ~= 0);
          SortG2 = unique(SortGen1(i,I2));
          I = zeros(1,250-size(SortG2,2));
          SortGen1(i,:) = [SortG2 I];
      end;
      Start = min(SortGen1(:,1));
      End = max(max(SortGen1(:,:)));
      SortGen1(:,1) = Start;
      for i = 1:size(SortGen1,1)
          H = find(SortGen1(i,:) == max(SortGen1(i,:)));
          SortGen1(i,H) = End;
      end;
      res1 = [];
      for i = 1:size(SortGen1,1)
          res = diff((SortGen1(i,:)));
          res1 = [res1 res(find(res > 4 & res < 10))];
      end;
      SQ_SIZE = mean(res1);
      SQ_SIZE = 8;
      for i = 1:size(SortGen1,1)
          k = 1;
          myGen = 0;
          myGen(1,1) = SortGen1(i,1);
          for j = 2:size(SortGen1,2)
              k = k + 1;
              if ((SortGen1(i,j) == 0))
                  break;
              end;
              if ((SortGen1(i,j) - SortGen1(i,j-1)) > SQ_SIZE)
                  N = round((SortGen1(i,j) - SortGen1(i,j-1))/SQ_SIZE) - 1;
                  for l = 1:N
                      myGen(1,k) = SortGen1(i,j-1) + l * SQ_SIZE;
                      k = k + 1;
                  end;
              end;
              myGen(1,k) = SortGen1(i,j);
          end;
          I = zeros(1,250-size(myGen,2));
          SortGen2(i,:) = [myGen I];
      end;
      for i = 1:size(SortGen2,1)
          T = find(SortGen2(i,:) ~=0);
          B(i) = size(T,2);
      end;
      [A,C] = max(B);
      for i = 1:size(SortGen2,1)
          k = 1;
          for j = 1:size(SortGen2,2)
              if (abs(SortGen2(i,j) - SortGen2(C,j)) < SQ_SIZE/2)
                  I = SortGen2(i,j);
                  SortGen3(i,k) = I;
                  k = k + 1;
              else 
                  SortGen3(i,k) = SortGen2(C,j);
                  k = k + 1;
              end;
          end;
      end;
      XX = repmat(nuniX',[1,A]);
      YY = SortGen3(:,1:A);     
       figure;imshow(IM1);
        hold on;plot(XX(:),YY(:),'rx');