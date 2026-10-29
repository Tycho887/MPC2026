% Testing qpsolver

Htest = 2*eye(2);
gtest = [-2;-4];

ltest = [0;0];
utest = [3;3];

Atest = [1 1];
bltest = -Inf;
butest = 2;

xinitTest = [0;0];

[x,info] = qpsolver(Htest,gtest,ltest,utest, Atest,bltest,butest,xinitTest)
