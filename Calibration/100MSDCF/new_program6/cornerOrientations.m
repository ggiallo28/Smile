function [v1, v2] = cornerOrientations(Ix2, Iy2, Ixy, p)
% The orientation vectors are the eigen vectors of the 
% structure tensor:
% [Ix^2  Ixy ]
% [Ixy   Iy^2]

a = Ix2(p(2), p(1));
b = Ixy(p(2), p(1));
c = Iy2(p(2), p(1));

% % Computing eigenvectors "by hand", because the eig() function behaves
% % differently in codegen.
% % Since the matrix is positive-semidefinite, its eigenvectors are
% % orthogonal. Compute the first eigenvector, then swap its elements and
% % negate the y-component to make the second one.
sm = a + c;
df = a - c;
adf = abs(df);
tb = b + b;
ab = abs(tb);

if adf > ab
    rt = adf * sqrt(1 + (ab/adf)^2);
elseif adf < ab
    rt = ab * sqrt(1 + (adf/ab)^2);
else
    rt = ab * sqrt(2);
end

if sm < 0
    sgn1 = -1;
else
    sgn1 = 1;
end

if df > 0
    cs = df + rt;
    sgn2 = 1;
else
    cs = df - rt;
    sgn2 = -1;
end

acs = abs(cs);
if acs > ab
    ct = -tb / cs;
    sn1 = 1 / sqrt(1 + ct * ct);
    cs1 = ct * sn1;
else
    if ab == single(0)
        cs1 = single(1);
        sn1 = single(0);
    else
        tn = -cs / tb;
        cs1 = 1 / sqrt(1 + tn * tn);
        sn1 = tn * cs1;
    end
end
if sgn1 == sgn2
    tn = cs1;
    cs1 = -sn1;
    sn1 = tn;
end

v1 = [-sn1, cs1];
v2 = [cs1, sn1];
end