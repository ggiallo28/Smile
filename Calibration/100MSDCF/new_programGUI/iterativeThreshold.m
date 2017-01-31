function P = iterativeThreshold(CC, GRIDv, sep_squares, Container, path)
toSave = false;
if toSave
    persistent blob_index;
    if isempty(blob_index)
        blob_index = 0;
    end
end
    P = [];
    %% Create Folder
    if toSave && ~exist([path,'iterative_img'], 'dir')
        curr_folder = pwd;
        cd(path)
        mkdir iterative_img;
        cd(curr_folder);
    end
    if toSave
        path_2 = [path,'iterative_img'];
    end
    for l=1:CC.NumObjects
        [cutR,cutC] = ind2sub(size(GRIDv),CC.PixelIdxList{l});
        tmp = sep_squares(min(cutR):max(cutR),min(cutC):max(cutC));
%        deb = zeros(size(GRIDv));
%	    deb(min(cutR):max(cutR),min(cutC):max(cutC)) = tmp;
%        imshowpair(deb,sep_squares);
        tt = 0.3*size(tmp,1)*size(tmp,2);
        condition = true;
        th = 0.9; store = cell(2,1);
        dot = findDots(tmp);
        %figure, imshow(tmp), hold on, scatter(dot(:,1),dot(:,2))
        P = [P;min(cutR)+dot(:,2),min(cutC)+dot(:,1)];    
        D = [];
        if toSave && ~exist([path_2,'/img',num2str(l)], 'dir')
            curr_folder = pwd;
            cd(path_2)
            mkdir(['img',num2str(l)]);
            cd(curr_folder);
        end
        if toSave
            path_3 = [path_2,'\img',num2str(l)];
        end
        if(isempty(dot))            
%            figure
            while(condition && th>0)
                bw_tmp = im2bw(tmp,th);
                th = th-0.005;
                CC_tmp = bwconncomp(bwareaopen(bw_tmp,10,4),4);
                if ~Container.isGUI
                    imshow(bw_tmp);
                end
%                     pause(0.1)
%                     th
                if(CC_tmp.NumObjects ==1 && size(CC_tmp.PixelIdxList{1},1)>tt)
                    condition = false;
                    [blob1y,blob1x] = ind2sub(size(tmp),store{1});
                    [blob2y,blob2x] = ind2sub(size(tmp),store{2});
                    D = zeros(5,size(blob1y,1)*(size(blob2y,1)-1)); idist = 1;
                    for k=1:size(blob1y,1)
                        for j=1:size(blob2y,1)
                            D(1,idist) = blob1x(k);
                            D(2,idist) = blob1y(k);
                            D(3,idist) = blob2x(j);
                            D(4,idist) = blob2y(j);
                            D(5,idist) = pdist([blob1x(k) blob1y(k); blob2x(j) blob2y(j)],'euclidean');
                            idist = idist+1;
                        end
                    end
                    if(isempty(D))
                         condition = true;
                    end
                end
                if(CC_tmp.NumObjects >= 2)
%   blob_index                     
% CC_tmp
% imshow(bw_tmp);
if toSave
                   imwrite(im2double(bw_tmp),[path_3,'/blob',num2str(blob_index),'.jpg']);
                   blob_index = blob_index +1;
end
                   store(1) = CC_tmp.PixelIdxList(1);
                   store(2) = CC_tmp.PixelIdxList(2);          
                end 
            end
            %% INVECE DIFARE STA MANFRINA USARE corner(bw_tmp,'Harris'
            if(isempty(D))
                YY = size(bw_tmp,1)*0.5;
                XX = size(bw_tmp,2)*0.5;
                P = [P;min(cutR)+YY,min(cutC)+XX];
            else
                [D(5,:), idd] = sort(D(5,:));
                D(1,:) = D(1,idd);
                D(2,:) = D(2,idd);
                D(3,:) = D(3,idd);
                D(4,:) = D(4,idd);
                if size(D,2)>=5
                    XX = round(mean([D(1,1:5),D(3,1:5)]));
                    YY = round(mean([D(2,1:5),D(4,1:5)]));
                else
                    XX = round(mean([D(1,:),D(3,:)]));
                    YY = round(mean([D(2,:),D(4,:)]));
                end
                P = [P;min(cutR)+YY,min(cutC)+XX];
            end
        end
    end
end

