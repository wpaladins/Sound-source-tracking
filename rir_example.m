function [h1,h2]  = rir_example(sp)
% 房屋环境：5m x 5m x 3m
% 笛卡尔坐标系中：
% 说话人初始位置：sp(1, 2)
% 麦克风位置：(1.2, 0.5)和(1.8, 0.5)
% 其他参数：声音频率16000Hz、声音速度340m/s、墙面回响系数0.4、采样数量4096、高度1.5m
% 调用请指定说话人位置，例如rir_example([1, 2])
if nargin < 1
   help rir_example
   return
end
s = [sp,1.5];
c = 340;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
r1 = [1.2 0.5 1.5];         % Receiver position [x y z] (m)
r2 = [1.8 0.5 1.5];         % Receiver position [x y z] (m)
% s = [1 2 1.5];              % Source position [x y z] (m)
L = [5 4 6];                % Room dimensions [x y z] (m)
beta = 0.4;                 % Reverberation time (s)
n = 4096;                   % Number of samples
% 生成rir
h1 = rir_generator(c, fs, r1, s, L, beta, n);
h2 = rir_generator(c, fs, r2, s, L, beta, n);
% figure(1);
% plot(h1);
% figure(2);
% plot(h2);
% audiowrite('rir1.wav',h1,fs)
% audiowrite('rir2.wav',h2,fs)
end