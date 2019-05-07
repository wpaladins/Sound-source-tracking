function [tdoaT] = tdoaT_generator(s,r1,r2,c)
% s 说话人位置，例如s=[1,2]
% r1,r2 麦克风1和麦克风2的位置，例如r1=[1.2,0.5];r2=[1.8,0.5]
% c 声速，默认为340m/s
% 调用时至少设置前三个参数

if nargin < 3
   help tdoaT_generator
   return
end

if nargin < 4
    c = 340;
end

d1 = sqrt((s(1) - r1(1) ) ^ 2 + (s(2) - r1(2) ) ^ 2);
d2 = sqrt((s(1) - r2(1) ) ^ 2 + (s(2) - r2(2) ) ^ 2);

tdoaT = (d1 - d2)/c;
end

