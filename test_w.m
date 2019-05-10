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

% 麦克风位置
s1r1 = [2.2,0.5]; s1r2 = [2.8,0.5];
s2r1 = [2.2,4.5]; s2r2 = [2.8,4.5];
s3r1 = [0.5,2.2]; s3r2 = [0.5,2.8];
s4r1 = [4.5,2.2]; s4r2 = [4.5,2.8];

% 测试方形区域各处粒子的权重
A = zeros(k,k,3);
a = linspace(0,5,k);
for i=1:k
    A(:,i,1) = a(i);
    A(i,:,2) = a(i);
end

% 假设此状态 的位置为
temp = [3.5085,1.8896];
[h1,h2] = rir_example(temp,s1r1,s1r2);
conv1 = conv(rawWav,h1);
conv2 = conv(rawWav,h2);
[s1_gccResult,s1_Nd] = gcc_phat_w(conv1,conv2);
% 得到s2_gccResult,s2_Nd
[h1,h2] = rir_example(temp,s2r1,s2r2);
conv1 = conv(rawWav,h1);
conv2 = conv(rawWav,h2);
[s2_gccResult,s2_Nd] = gcc_phat_w(conv1,conv2);
% 得到s3_gccResult,s3_Nd
[h1,h2] = rir_example(temp,s3r1,s3r2);
conv1 = conv(rawWav,h1);
conv2 = conv(rawWav,h2);
[s3_gccResult,s3_Nd] = gcc_phat_w(conv1,conv2);
% 得到s4_gccResult,s4_Nd
[h1,h2] = rir_example(temp,s4r1,s4r2);
conv1 = conv(rawWav,h1);
conv2 = conv(rawWav,h2);
[s4_gccResult,s4_Nd] = gcc_phat_w(conv1,conv2);

for i=1:k
    for j=1:k
        s1_tdoaT = tdoaT_generator([A(i,j,1),A(i,j,2)],s1r1,s1r2);
        s2_tdoaT = tdoaT_generator([A(i,j,1),A(i,j,2)],s2r1,s2r2);
        s3_tdoaT = tdoaT_generator([A(i,j,1),A(i,j,2)],s3r1,s3r2);
        s4_tdoaT = tdoaT_generator([A(i,j,1),A(i,j,2)],s4r1,s4r2);
        
        A(i,j,3) = particle_weight_generator(s1_gccResult,s1_Nd,...
                                                s2_gccResult,s2_Nd,...
                                                s3_gccResult,s3_Nd,...
                                                s4_gccResult,s4_Nd,...
                                                fs,...
                                                s1_tdoaT,...
                                                s2_tdoaT,...
                                                s3_tdoaT,...
                                                s4_tdoaT);
    end
end
figure(1);
mesh(A(:,:,1),A(:,:,2),A(:,:,3));