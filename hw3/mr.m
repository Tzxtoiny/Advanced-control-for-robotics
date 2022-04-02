M = sym([1, 0,  0, 10;...
            0, 1,  0, 0;...
            0, 0, 1, -20;...
            0,0,0,1]);
Slist =  sym([0, 1,  0,  0, 0,    0;0, 0, 1, 0, 0, 0]')
syms t positive;
thetalist = sym([t ,t])  
% A = MatrixExp6(VecTose3(Slist(:, 1))* thetalist(1))
ex = FKinSpace(M,Slist,thetalist)
ex= simplify(ex)
pr = latex(ex)
pt = latex(ex(1:3,4))