% Environment: Matlab r2017a for mac
clc;
clear;

randn('seed',50);  %#ok<*RAND>

dT = 0.032;
SNR = 40; % 信噪比

% 麦克风位置
s1r1 = [2.2,0.5]; s1r2 = [2.8,0.5];
s2r1 = [2.2,4.5]; s2r2 = [2.8,4.5];
s3r1 = [0.5,2.2]; s3r2 = [0.5,2.8];
s4r1 = [4.5,2.2]; s4r2 = [4.5,2.8];

% 1. 分帧
rawWav = wavread('raw.wav'); %#ok<*DWVRD>，读取的wav文件，频率与Image模型统一，此处为16000HZ
fs = 16000;
fram_time = 0.032;
fram_step_time = 0.032;
win = 'hanning';
%分帧结果：xn_frams(:,1)为第一帧，以此类推
xn_frams = framing(rawWav,fs,fram_time,fram_step_time,win);

% 2. 产生所有的真实状态，以及观测值
% 2.1. 真实状态
T = 50; % 说话人位置改变次数，帧数
X = zeros(50,4);
R = 1.5;
t = linspace(pi,2*pi,T);
X(:,1:2) = [(2.5 + R*cos(t))',(3 + R*sin(t))'];

for i=1:T-1 % 速度只有T为1:T-1有
    X(i,3) = (X(i+1,1) - X(i,1) ) / dT; % x方向的速度
    X(i,4) = (X(i+1,2) - X(i,2) ) / dT; % y方向的速度
end % 速度完成

% 3. 粒子滤波
% 3.1. 各项参数
numSamples = 25;
Xpf=zeros(numSamples,T,4); % 行代表某一个粒子，列代表某一个时刻，值为这个粒子在当前时刻的 粒子滤波后（除初始化）的 状态
Xparticles=zeros(numSamples,T,4); % 行代表某一个粒子，列代表某一个时刻，值为这个粒子在当前时刻的 粒子滤波前的 状态
Zpre_pf=zeros(numSamples,T); % 行代表某一个粒子，列代表某一个时刻，值为这个粒子在当前时刻的 观测值
weight=zeros(numSamples,T); % 行代表某一个粒子，列代表某一个时刻，值为这个粒子在当前时刻的 权重
QQQ = 0.01; % 高斯滤波的权值的平方[---待确认---]

Tpf = zeros(1,T);

% 粒子初始化[---待确认---]这里使用了真实值
Xpf(:,1,:)=X(1,:)+sqrt(QQQ)*randn(numSamples,4); % 初始粒子状态，使用高斯滤波对真实状态处理产生

% 画图显示粒子分布情况
figure(1);
plot(Xpf(:,1,1),Xpf(:,1,2),'g.',X(1,1),X(1,2),'r.');
axis([0 5 0 5]);
saveas(1,'./jpg/1.jpg');

% 粒子滤波核心循环
for k=2:T
    % 得到s1_gccResult,s1_Nd
    [h1,h2] = rir_example(X(k,[1,2]),s1r1,s1r2);
    conv1 = conv(rawWav,h1);
    conv2 = conv(rawWav,h2);
    %高斯白噪声 awgn1 = awgn(conv1,SNR);
    %高斯白噪声 awgn2 = awgn(conv2,SNR);
    [s1_gccResult,s1_Nd] = gcc_phat_w(conv1,conv2);
    % 得到s2_gccResult,s2_Nd
    [h1,h2] = rir_example(X(k,[1,2]),s2r1,s2r2);
    conv1 = conv(rawWav,h1);
    conv2 = conv(rawWav,h2);
    %高斯白噪声 awgn1 = awgn(conv1,SNR);
    %高斯白噪声 awgn2 = awgn(conv2,SNR);
    [s2_gccResult,s2_Nd] = gcc_phat_w(conv1,conv2);
    % 得到s3_gccResult,s3_Nd
    [h1,h2] = rir_example(X(k,[1,2]),s3r1,s3r2);
    conv1 = conv(rawWav,h1);
    conv2 = conv(rawWav,h2);
    %高斯白噪声 awgn1 = awgn(conv1,SNR);
    %高斯白噪声 awgn2 = awgn(conv2,SNR);
    [s3_gccResult,s3_Nd] = gcc_phat_w(conv1,conv2);
    % 得到s4_gccResult,s4_Nd
    [h1,h2] = rir_example(X(k,[1,2]),s4r1,s4r2);
    conv1 = conv(rawWav,h1);
    conv2 = conv(rawWav,h2);
    %高斯白噪声 awgn1 = awgn(conv1,SNR);
    %高斯白噪声 awgn2 = awgn(conv2,SNR);
    [s4_gccResult,s4_Nd] = gcc_phat_w(conv1,conv2);
    
    % 通过对 上一时刻的粒子状态 使用 状态方程，得到 这一时刻的粒子状态
    for i=1:numSamples
        QQ=0.01; % 网[---待确定---]
        net=sqrt(QQ)*randn(4,1); % 网[---待确定---]
        temp = zeros(1,4); % temp为上一时刻 当前粒子的 状态
        temp(1) = Xpf(i,k-1,1);
        temp(2) = Xpf(i,k-1,2);
        temp(3) = Xpf(i,k-1,3);
        temp(4) = Xpf(i,k-1,4);
        Xparticles(i,k,:) = langevin(temp' )' + net';% 网[---待确定---]
    end
    
    % 画图查看粒子变化
    % 绿色：当前粒子位置；红色：当前真实位置；蓝色：上个真实位置；黄色：粒子运动前的位置
    figure(2)
    plot(Xparticles(:,k,1),Xparticles(:,k,2),'g.',X(k,1),X(k,2),'r.',X(k-1,1),X(k-1,2),'b.',Xpf(:,k-1,1),Xpf(:,k-1,2),'y.');
    axis([0 5 0 5]);
    title(num2str(k));
    jpg = strcat('./jpg/',num2str(k));
    jpg = strcat(jpg,'.jpg');
    saveas(2,jpg);
    
    tic;
    
    % 粒子权重处理
    for i=1:numSamples
        temp = zeros(1,2);
        temp(1) = Xparticles(i,k,1);
        temp(2) = Xparticles(i,k,2);
        s1_tdoaT = tdoaT_generator(temp,s1r1,s1r2);
        s2_tdoaT = tdoaT_generator(temp,s2r1,s2r2);
        s3_tdoaT = tdoaT_generator(temp,s3r1,s3r2);
        s4_tdoaT = tdoaT_generator(temp,s4r1,s4r2);
        
        % 更新粒子权重
        weight(i,k) = particle_weight_generator(s1_gccResult,s1_Nd,...
                                                s2_gccResult,s2_Nd,...
                                                s3_gccResult,s3_Nd,...
                                                s4_gccResult,s4_Nd,...
                                                fs,...
                                                s1_tdoaT,...
                                                s2_tdoaT,...
                                                s3_tdoaT,...
                                                s4_tdoaT);
    end
    
    weight(:,k)=weight(:,k)./sum(weight(:,k));

    % 粒子 重新采样
    outIndex = multinomialR(weight(:,k));
    
    % 产生粒子滤波后的 所有粒子
    Xpf(:,k,:)= Xparticles(outIndex,k,:);
    
    Tpf(k) = toc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 最终目标的计算
Xmean_x_pf=mean(Xpf(:,:,1));
Xmean_y_pf=mean(Xpf(:,:,2));

bins=20;
Xmap_x_pf=zeros(T,1);
Xmap_y_pf=zeros(T,1);
for k=1:T
    [p,pos]=hist(Xpf(:,k,1),bins); % [n,xout]=hist(Y,nbins); % nbins是间隔数，也就是说我们应该统计多少个间隔，n是每一个区间的个数，xout是区间的中心位置
    map=find(p==max(p));
    Xmap_x_pf(k,1)=pos(map(1));  
end
for k=1:T
    [p,pos]=hist(Xpf(:,k,2),bins);
    map=find(p==max(p));
    Xmap_y_pf(k,1)=pos(map(1));  
end
Xstd_x_pf = zeros(1,T);
Xstd_y_pf = zeros(1,T);

% 计算RMSE
Xdiff_pf = zeros(1,T);
sum = 0;
for i=1:T
    temp1 = X(i,1) - Xmean_x_pf(i);
    temp2 = X(i,2) - Xmean_y_pf(i);
    Xdiff_pf(i) = temp1^2 + temp2^2;
    sum = sum + Xdiff_pf(i);
end
RMSE_pf = sqrt(1/T * sum);

figure(12);clf;
k=1:1:T;
plot(k,X(:,1),'b',k,Xmean_x_pf,'r'); 
legend('True X value','Posterior mean estimate of X');
xlabel('Time frame','fontsize',15);
ylabel('Value','fontsize',15);
axis([0 50 0 5]);
title('The comparison between True X value and Posterior mean estimate of X','fontsize',15);
saveas(12,'./jpg/X估计值与真值.jpg'); % 保存

figure(22);clf;
k=1:1:T;
plot(k,X(:,2),'b',k,Xmean_y_pf,'r'); 
legend('True Y value','Posterior mean estimate of Y');
xlabel('Time frame','fontsize',15);
ylabel('Value','fontsize',15);
axis([0 50 0 5]);
title('The comparison between True Y value and Posterior mean estimate of Y','fontsize',15);
saveas(22,'./jpg/Y估计值与真值.jpg'); % 保存

figure(16);
k=1:1:T;
plot(k,Xdiff_pf,'-');
xlabel('Time frame','fontsize',15);
ylabel('Deviation','fontsize',15);
title('Deviation of each time frame');
axis([0,T,0,5] );
saveas(16,'./jpg/RMSE.jpg'); % 保存

% 保存数据
save('./mat/pf.mat','Tpf','Xmean_x_pf','Xmean_y_pf','Xdiff_pf','RMSE_pf');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 程序结束提醒
disp('Done');