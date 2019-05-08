% Environment: Matlab r2017a for mac
clc;
clear;

% 1. 参数
randn('seed',50);  %#ok<*RAND>

T = 50;
dT = 0.032;

rawWav = wavread('raw.wav'); %#ok<*DWVRD>，读取的wav文件，频率与Image模型统一，此处为16000HZ

% 测试方形区域内各处的P
A = zeros(26,26,3);
a = linspace(0,5,26);
for i=1:26
    A(:,i,1) = a(i);
    A(i,:,2) = a(i);
end

for i=1:26
    for j=1:26
        [h1,h2] = rir_example([A(i,j,1),A(i,j,2)] );
        conv1 = conv(rawWav,h1);
        conv2 = conv(rawWav,h2);
        tdoaT = tdoaT_generator([A(i,j,1),A(i,j,2)],[1.2,0.5],[1.8,0.5]); % 获取理论TDOA，麦克风位置：(1.2, 0.5)和(1.8, 0.5)
        A(i,j,3) = gcc_phat_P(conv1,conv2,tdoaT);
    end
end
figure(1);
mesh(A(:,:,1),A(:,:,2),A(:,:,3));