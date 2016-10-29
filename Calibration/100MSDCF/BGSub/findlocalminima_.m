%-------------------Introduction------------------------------%
%
% This function finds dominant peaks and valleys of an image 
% histogram. Minimum peak distance and averaging filter size can be 
% changed to suit your need. 
% Although, the function is originally written for mages, 
% there is no reason why this can't be used on one 
% dimensional data. 
%
% I is the image segment in uint8 format
% mpd is the minimum peak distance DEFAULT = 10;
% windowSize is the averaging filter size to 
% average the histogram  DEFAULT = 5
% refine is to refine the peaks based on area DEFAULT = 0
%
% USAGE: 
% [peaks minima_low minima_high] = findlocalminima(I,10,5,0);
% OR
% [peaks minima_low minima_high] = findlocalminima(I);
%
% This is the function that is used for peak detection in 
% my paper, 
% "De Silva, D.V.S.; Fernando, W.A.C.; Kodikaraarachchi, H.;
% Worrall, S.T.; Kondoz, A.M.; ,"Adaptive sharpening of depth
% maps for 3D-TV," Electronics Letters , vol.46, no.23, 
% pp.1546 -1548, November 11 2010
%
% Kindly cite this paper if you use this in your work. 
%
% Varuna De Silva, University of Surrey, Guildford, UK
% varunax@gmail.com
%
%-------------------------------------------------------------%
function [output_peak, output_minima_low, output_minima_high, hist_size, output_minima_mid] = findlocalminima(I,mpd,windowSize,refine,isPeriodic)

if nargin < 2
%minimum peak distance
mpd = 10;
%smoothing windowSize
windowSize = 5;
%Local minima refine
refine = 0;
end

%threshold is set at least 
thresh = size(I,1)*size(I,2)*0.01;

%Find the histogram
h = imhist(I);
imshow(I)
%Smooth the histogram to remove local variations
h1 = filter(ones(1,windowSize)/windowSize,1,h);

d = diff(h1);

sign_previous = [0;];
sign = zeros(255,1);

for i = 1:255
    if(d(i) ~= 0)
        sign(i,1) = d(i)/abs(d(i));
    else
        sign(i,1) = 0;
    end
end

sign_previous = [sign_previous ; sign];

%Find peaks with a minimum peakdistance of a given value

% Find peaks
p_pos = 0;
p_neg = 0;
pnum = 0;
for i=2:256
    if sign_previous(i)<sign_previous(i-1)
        if sign_previous(i) == 0
            p_pos = i-1;
        else if sign_previous(i)<0 && sign_previous(i-1)==0
                p_neg = i-1;
                pnum = pnum+1;
                locs(pnum) = round((p_pos+p_neg)/2);
            else
                pnum=pnum+1;
                locs(pnum)=i-1;
             end
        end
    end
end

peak = zeros(size(locs,2),1);
minima_low = zeros(size(locs,2),1);
minima_high = zeros(size(locs,2),1);
area = zeros(size(locs,2),1);

for i = 1:size(locs,2)
    peak(i) = locs(i);
end

%Find Local minima
num = 1;
lmin = [];

for i=2:256
    if sign_previous(i)>sign_previous(i-1)
        lmin(num)=i-1;
        num=num+1;
    end
end

%find local minima_low and minima_high
for i = 1:size(locs,2)
    ind_min = find(lmin-locs(i)<0);
    ind_max = find(lmin-locs(i)>0);
    
    size_min = size(ind_min,2);
    size_max = size(ind_max,2);
    
    minima_low(i) = 1;
    minima_high(i) = 255;

    if size_min>=1
        minima_low(i) = lmin(ind_min(size_min));
    end

    if size_max>=1
        minima_high(i) = lmin(ind_max(1));
    end
end

% refine local minima
if refine>0
    for i = 1:size(locs,2)-1
        if(minima_high(i)< minima_low(i+1))
        if( minima_high(i)~=minima_low(i+1) && min(h1(minima_high(i):minima_low(i+1))>0) )
            minima = floor((minima_high(i)+ minima_low(i+1))/2);
            minima_low(i+1) = minima;
            minima_high(i) = minima;
        end
        end
    end
end
   
num = 0;
outpeak = 0;

%refine peaks based on area threshold & mpd
for i = 1:size(locs,2)
    area(i) = sum(h1(minima_low(i):minima_high(i)));
    if outpeak >= 1
        if area(i)>thresh 
            if locs(i)-locs(outpeak)>mpd
                num = num+1;
                outpeak = i;
                output_peak(num) = peak(i);
                output_minima_low(num) = minima_low(i);
                output_minima_high(num) = minima_high(i);
                
            else
                if area(i)>area(outpeak) 
                    outpeak = i;
                    output_peak(num) = peak(i);
                    %output_minima_low(num) = minima_low(i);
                    output_minima_high(num) = minima_high(i); 
                else
                    output_minima_high(num) = minima_high(i);
                end
            end
        end
    else
        if area(i)>thresh
            num = num+1;
            outpeak = i;
            output_peak(num) = peak(i);
            output_minima_low(num) = minima_low(i);
            output_minima_high(num) = minima_high(i);
        end
    end
end
%% Check Periodicity
hist_size = size(h,1);
if  isPeriodic  
%     for i=1:size(output_peak,2)-1
%         for j=i+1:size(output_peak,2)
%             n1 = abs(output_peak(i)-output_peak(j));
%             n2 = hist_size-abs(output_peak(i)-output_peak(j));
%             dist = min(n1,n2);
%             if(dist<=mpd)
%                 if(dist == n1)
%                     output_peak(i) = mod(output_peak(i)+dist/2,hist_size);
%                     output_peak(j) = output_peak(i);
%                 else
%                     output_peak(j) = mod(output_peak(j)+dist/2,hist_size);
%                     output_peak(i) = output_peak(j);
%                 end
%             end
%         end
%     end
%     output_peak = unique(output_peak', 'rows')';
%     
%     for i=1:size(output_peak,2)
%         for j=1:size(output_minima_high,2)
%             n1 = abs(output_peak(i)-output_minima_high(j));
%             n2 = hist_size-abs(output_peak(i)-output_minima_high(j));
%             dist = min(n1,n2);
%             if(dist<=mpd/2)
%                 output_minima_high(j) = -1;
%             end
%         end
%      end
%     
%      for i=1:size(output_peak,2)
%         for j=1:size(output_minima_low,2)
%             n1 = abs(output_peak(i)-output_minima_low(j));
%             n2 = hist_size-abs(output_peak(i)-output_minima_low(j));
%             dist = min(n1,n2);
%             if(dist<=mpd/2)
%                  output_minima_low(j) = -1;
%             end
%         end
%      end  
%      
     output_minima_low(output_minima_low==-1) = [];
     output_minima_high(output_minima_high==-1) = [];
     output_minima_low(output_minima_low<=output_peak(1)) = [];
     output_minima_low(output_minima_low>=output_peak(end)) = [];
     output_minima_high(output_minima_high<=output_peak(1)) = [];
     output_minima_high(output_minima_high>=output_peak(end)) = [];
end
output_minima_mid = (max([output_minima_low;output_minima_high]) - min([output_minima_low;output_minima_high]))/2 + min([output_minima_low;output_minima_high]);
plot(h1)
xlabel('Pixel Value');
ylabel('Number of Pixels');
hold on
plot(output_peak,h1(output_peak),'*');
plot(output_minima_low,h1(output_minima_low),'or');
plot(output_minima_high,h1(output_minima_high),'ob');
plot(round(output_minima_mid),h1(round(output_minima_mid)),'og');
hold off
axis tight
legend('hist','peak','min_low','min_high','min_mid');

end


