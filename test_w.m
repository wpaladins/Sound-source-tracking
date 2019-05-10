% Environment: Matlab r2017a for mac
clc;
clear;

% 1. 参数
randn('seed',50);  %#ok<*RAND>

T = 50;
dT = 0.032;
fs = 16000;
k = 100;

rawWav = wavread('raw.wav'); %#ok<*DWVRD>，读取的wav文件，频率与Image模型统一，此处为16000HZ

% 测试方形区域各处粒子的权重
A = zeros(k,k,3);
a = linspace(0,5,k);
for i=1:k
    A(:,i,1) = a(i);
    A(i,:,2) = a(i);
end

% 假设此状态 的位置为(1,3)
[h1,h2] = rir_example([1,3]);
conv1 = conv(rawWav,h1);
conv2 = conv(rawWav,h2);
[gccResult,Nd] = gcc_phat_w(conv1,conv2);

for i=1:k
    for j=1:k
        tdoaT = tdoaT_generator([A(i,j,1),A(i,j,2)],[1.2,0.5],[1.8,0.5]); % 获取理论TDOA，麦克风位置：(1.2, 0.5)和(1.8, 0.5)
        % disp(tdoaT);
        A(i,j,3) = particle_weight_generator(gccResult,Nd,fs,tdoaT);
    end
end
figure(1);
mesh(A(:,:,1),A(:,:,2),A(:,:,3));