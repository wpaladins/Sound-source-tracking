function [p] = CostFunction(R,yk,yp)
% 计算粒子位置，与预测位置的适应度
p = exp(-1/(2*R) * (yk - yp) * (yk - yp)' ); % 1/(2*pi*R) * 
end