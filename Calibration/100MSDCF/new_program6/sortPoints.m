function [LTp, RTp, LBp, RBp] = sortPoints(P1,P2,P3,P4,Cy)
    X = [P1(1) P2(1) P3(1) P4(1)];
    Y = [P1(2) P2(2) P3(2) P4(2)];
    iT = find(Y<Cy);
    iB = find(Y>Cy);
    
    [val, i] = min(X(iT));
    YY = Y(iT);
    LTp = [val, YY(i)];
    [val, i] = max(X(iT));
    YY = Y(iT);
    RTp = [val, YY(i)];
    [val, i] = min(X(iB));
    YY = Y(iB);
    LBp = [val, YY(i)];
    [val, i] = max(X(iB));
    YY = Y(iB);
    RBp = [val, YY(i)];
end