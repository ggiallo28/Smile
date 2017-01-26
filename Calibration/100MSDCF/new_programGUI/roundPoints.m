function [P1, P2, P3, P4]  =  roundPoints(vect,dim)
    vect = round(vect);
    vect(vect<1) = 1;
    vect_x = vect(:,1);
    vect_x(vect_x>dim(2)) = dim(2);
    vect_y = vect(:,2);
    vect_y(vect_y>dim(1)) = dim(1);
    P1 = [vect_x(1), vect_y(1)];
    P2 = [vect_x(2), vect_y(2)];
    P3 = [vect_x(3), vect_y(3)];
    P4 = [vect_x(4), vect_y(4)];
end