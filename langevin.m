function [X]  = langevin(Xp)
% 郎之万模型
% 输入Xp = [x,y,xv,yv]'
% 输出X = [x,y,xv,yv]'
% 注：输入输出均为列向量
%
% 三个变量的人工初始化
vA = 10;
dT = 0.032;
vB = 0.001;
% 求a,b,I2
a = exp(-vA * dT);
b = vB * (1 - a^2)^(1/2);
I2 = zeros(2,2);
I2(1,1) = 1;
I2(2,2) = 1;
% matA和matB初始化
matA = zeros(4,4);
matB = zeros(4,4);
% matA填值
matA(1,1) = 1;
matA(2,2) = 1;
matA([1,2],[3,4]) = kron(a*dT, I2);
matA([3,4],[3,4]) = kron(a, I2);
% matB填值
matB([1,2],[1,2]) = kron(b*dT, I2);
matB([3,4],[3,4]) = kron(b, I2);
% 过程噪声
u = randn(4,1);
% 测试
X = matA * Xp + matB * u;
return

end
