close all; clear all; clc;
path = '../foto/Evaluation/';
checker_vector = reshape([[0,0,0;255,0,255];[0,0,0;0,255,255];[0,0,0;255,255,0];[255,255,255;255,0,0];[255,255,255;0,255,0];[255,255,255;0,0,255]],[2,6,3]);
checker_center = [0.5*size(checker_vector,2),0.5*size(checker_vector,2)+1];
folders = dir(path);
for ff=1:size(folders,1)
    if ~(strcmp(folders(ff).name,'..') || strcmp(folders(ff).name,'.'))
        path_2 = [path,folders(ff).name,'/'];
        folders_2 = dir(path_2);
        for ff2=1:size(folders_2,1)
            if ~(strcmp(folders_2(ff2).name,'..') || strcmp(folders_2(ff2).name,'.'))
                path_3 = [path_2,folders_2(ff2).name,'/'];
                disp([path_3,'....']);
                if exist([path_3,'confidence'], 'file') == 2
                    parpar.confidence = 3.2;
                else
                    parpar.confidence = 2.8;
                end
                if exist([path_3,'mpd'], 'file') == 2
                     parpar.mpd = 25;
                else
                     parpar.mpd = 30;
                end
                load([path_3,'bb.mat']);
                parpar.bb = bb;
                if exist([path_3,'toflip'], 'file')
                    parpar.checker_vector = fliplr(reshape([[0,0,0;255,0,255];[0,0,0;0,255,255];[0,0,0;255,255,0];[255,255,255;255,0,0];[255,255,255;0,255,0];[255,255,255;0,0,255]],[2,6,3]));
                else
                    parpar.checker_vector = reshape([[0,0,0;255,0,255];[0,0,0;0,255,255];[0,0,0;255,255,0];[255,255,255;255,0,0];[255,255,255;0,255,0];[255,255,255;0,0,255]],[2,6,3]);
                end
                % parpar.checker_vector = checker_vector;
                parpar.checker_center = checker_center;
                save([path_3,'parpar.mat'],'parpar');
            end
        end
    end
end