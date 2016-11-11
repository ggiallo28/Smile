function obj = fix_positions(x, obj, BW, toFix)
    for i = 1:size(x,1)
        X = obj.chess(i).center_x;
        Y = obj.chess(i).center_y;
        ind = find(obj.chess(i).center_x==0);
        [id_row,id_col]=ind2sub(size(X),ind);
        X(id_row,:) = 0;
        Y(id_row,:) = 0;
        axis = 1:5; axis(id_row) = [];
        for k=1:size(X,2)
            value = X(:,k); value(X(:,k)==0) = [];
            if(size(value,1) <=1)
                continue
            end
            X(:,k) = pchip(axis,value,1:5);
            value = Y(:,k); value(Y(:,k)==0) = [];
            if(size(value,1) <=1)
                continue
            end
            Y(:,k) = pchip(axis,value,1:5); 
        end
        if(size(find(X==0),1)>0 || size(find(Y==0),1)>0)
            continue
        end
        % trovo la posizione giusta degli zeri in funzione di quelle che sono
        % le tessere mancanti
        for k=1:size(id_row,1)
            row_Xchess = obj.chess(i).center_x(id_row(k),:);
            row_Ychess = obj.chess(i).center_y(id_row(k),:);
            row_Xchess(row_Xchess ==0) = [];
            row_Ychess(row_Ychess ==0) = [];
            row_X = X(id_row(k),:);
            index = 1:size(row_Xchess,2);
            for j=1:size(row_Xchess,2)
                curr_val_x = row_Xchess(j);
                curr_val_y = row_Ychess(j);
                row = abs(row_X-curr_val_x);
                idx = find(row== min(row));
                X(id_row(k),idx) = curr_val_x;
                Y(id_row(k),idx) = curr_val_y;
                index(index == idx) = [];
            end
            X(id_row(k),index) = 0;
            Y(id_row(k),index) = 0;
        end  
        obj.chess(i).center_x = X;
        obj.chess(i).center_y = Y;
        if toFix
            figure,imshow(BW); hold on;
            ind = find(obj.chess(i).center_x==0);
            while(size(ind,1)>0)
                %       mod(4,2)   % Even  number
                %       ans = 0
                %       mod(5,2)   % even number 
                %       ans =1
                [id_row,id_col]=ind2sub(size(obj.chess(i).center_x),ind(1)); %Prendo sempre il primo
                X = obj.chess(i).center_x(:,id_col);
                Y = obj.chess(i).center_y(id_row,:);
                xv = 1:size(X,1);
                yv = 1:size(Y,2);    
                X_pchip = X; X_pchip(X==0)=[];
                Y_pchip = Y; Y_pchip(Y==0)=[];
                if(size(X_pchip,1) == 1) % Se abbiamo solo un punto switcha sulle righe/colle
                    X = obj.chess(i).center_x(id_row,:);
                    xv = 1:size(X,1);   
                    X_pchip = X; X_pchip(X==0)=[];
                end
                if(size(Y_pchip,2) == 1)
                    Y = obj.chess(i).center_y(:,id_col);
                    yv = 1:size(Y,1);    
                    Y_pchip = Y; Y_pchip(Y==0)=[];
                end
                xv_pchip = xv; xv_pchip(X==0)=[];
                yv_pchip = yv; yv_pchip(Y==0)=[];
                if(size(yv_pchip,2) == 1 || size(xv_pchip,2) == 1) % Se quanto fatto prima non risolve cancella la colonna
                    continue;
                end
                val_x =  pchip(xv_pchip,X_pchip',xv);
                o = obj.chess(i).center_x;
                obj.chess(i).center_x(id_row,id_col) = val_x(id_row);
                o =  obj.chess(i).center_y;
                val_y = pchip(yv_pchip,Y_pchip,yv);
                obj.chess(i).center_y(id_row,id_col) = val_y(id_col); 
                ind = find(obj.chess(i).center_x==0);
                disp(['fix',num2str(i)]);
            end
            for k=1:size(obj.chess(i).center_x,1)
                for j=1:size(obj.chess(i).center_x,2)
                    scatter(obj.chess(i).center_x(k,j),obj.chess(i).center_y(k,j));
                end
            end 
        end
    end
end