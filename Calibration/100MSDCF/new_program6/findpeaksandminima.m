function [pks, output_mid, hist_size ] = findpeaksandminima(image,windowSize,mpd)
    h = imhist(image);
    
    h(40) = h(40)+10000;
    h(110) = h(110)+10000;
    
    h1 = filter(ones(1,windowSize)/windowSize,1,[h;h;h]);
    h1 = h1(size(h,1)+1:2*size(h,1)+1);
    [pks,locs] = findpeaks(h1,'MinPeakDistance',mpd);
    [~,locs_m] = findpeaks(-h1,'MinPeakDistance',mpd);
    locs(pks<1000) = [];
    hist_size = size(h,1);
    if(size(locs,1)==5)
        output_mid(1) = locs(1)+0.5*(locs(2)-locs(1));
        output_mid(2) = locs(2)+0.5*(locs(3)-locs(2));
        output_mid(3) = locs(3)+0.5*(locs(4)-locs(3));
        output_mid(4) = locs(4)+0.5*(locs(5)-locs(4));
    end
    if(size(locs,1)==4)
            output_mid(1) = locs(1)+0.5*(locs(2)-locs(1));
            output_mid(2) = locs(2)+0.5*(locs(3)-locs(2));
            output_mid(3) = locs(3)+0.5*(locs(4)-locs(3));
            output_mid(4) = mod(locs(4)+0.5*abs((256+locs(1))-locs(4)),hist_size);
            output_mid = sort(output_mid);
    end
    %for i=1:size(output_mid,

    plot(h1)
    xlabel('Pixel Value');
    ylabel('Number of Pixels');
    hold on
    plot(locs,h1(locs),'*');
    plot(locs_m,h1(locs_m),'or');
    assert(size(locs,1) == 4 | size(locs,1) == 5);
    plot(round(output_mid),h1(round(output_mid)),'og');
    hold off
    axis tight
    legend('hist','peak','min_low','min_mid');    
end

