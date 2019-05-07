% Environment: Matlab r2017a for mac
clc;
clear;

randn('seed',50);  %#ok<*RAND>

dT = 0.032;

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
R = 1.5;
t = linspace(pi,2*pi,T);
X = [(2.5 + R*cos(t))',(3 + R*sin(t))']; % 位置完成
for i=1:T-1 % 速度只有T为1:T-1有
    X(i,3) = (X(i+1,1) - X(i,1) ) / dT; % x方向的速度
    X(i,4) = (X(i+1,2) - X(i,2) ) / dT; % y方向的速度
end % 速度完成

% 2.2. 观测值
Z=zeros(T,1); % 观测值为四个候选TDOA中最大的一个
for k=1:T
    % 将X(k,:)作为说话人位置发给rir_example
    [h1,h2] = rir_example(X(k,[1,2]));
    % 将h1,h2分别与真实信号做卷积
    conv1 = conv(xn_frams(:,k),h1);
    conv2 = conv(xn_frams(:,k),h2);
    % 添加高斯白噪声
    % awgn1 = awgn(conv1,SNR);
    % awgn2 = awgn(conv2,SNR);
    tdoaT = tdoaT_generator(X(k,[1,2]),[1.2,0.5],[1.8,0.5]); % 获取理论TDOA，麦克风位置：(1.2, 0.5)和(1.8, 0.5)
    Z(k,1) = gcc_phat(conv1,conv2,tdoaT); % 高斯白噪声未添加
end

% 2.3. 真实状态信息展示
% figure;
% plot(X(:,1),X(:,2),'b.');
% axis([0 5 0 5]);
% 2.4. 观测值信息展示
% figure;
% plot(Z,'r');

% 3. 粒子滤波
% 3.1. 各项参数
numSamples = 100;
Xpf=zeros(numSamples,T,4); % 行代表某一个粒子，列代表某一个时刻，值为这个粒子在当前时刻的 粒子滤波后（除初始化）的 状态
Xparticles=zeros(numSamples,T,4); % 行代表某一个粒子，列代表某一个时刻，值为这个粒子在当前时刻的 粒子滤波前的 状态
Zpre_pf=zeros(numSamples,T); % 行代表某一个粒子，列代表某一个时刻，值为这个粒子在当前时刻的 观测值
weight=zeros(numSamples,T); % 行代表某一个粒子，列代表某一个时刻，值为这个粒子在当前时刻的 权重
QQQ = 0.1; % 高斯滤波的权值的平方[---待确认---]

% 粒子初始化[---待确认---]这里使用了真实值
Xpf(:,1,:)=X(1,:)+sqrt(QQQ)*randn(numSamples,4); % 初始粒子状态，使用高斯滤波对真实状态处理产生

% 画图显示粒子分布情况
figure(1);
plot(Xpf(:,1,1),Xpf(:,1,2),'g.',X(1,1),X(1,2),'r.');
axis([0 5 0 5]);
saveas(1,'./jpg/1.jpg');

for i=1:numSamples % 产生每个粒子对应的 观测值
    temp = zeros(1,2); % temp为粒子所在位置
    temp(1) = Xpf(i,1,1);
    temp(2) = Xpf(i,1,2);
    [h1,h2] = rir_example(temp);
    conv1 = conv(xn_frams(:,1),h1);
    conv2 = conv(xn_frams(:,1),h2);
    tdoaT = tdoaT_generator(temp,[1.2,0.5],[1.8,0.5]);
    Zpre_pf(i,1) = gcc_phat(conv1,conv2,tdoaT); % 无需噪声？[---待确认---]
end

% 粒子滤波核心循环
for k=2:T
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
    jpg = strcat('./jpg/',num2str(k));
    jpg = strcat(jpg,'.jpg');
    saveas(2,jpg);
    
    % 粒子权重处理
    for i=1:numSamples
        % 得到粒子在当前时刻的观测值
        temp = zeros(1,2);
        temp(1) = Xparticles(i,k,1);
        temp(2) = Xparticles(i,k,2);
        [h1,h2] = rir_example(temp);
        conv1 = conv(xn_frams(:,k),h1);
        conv2 = conv(xn_frams(:,k),h2);
        tdoaT = tdoaT_generator(temp,[1.2,0.5],[1.8,0.5]);
        Zpre_pf(i,k) = gcc_phat(conv1,conv2,tdoaT);
        % 求粒子权重
        weight(i,k) = exp(-.5*0.001^(-1)*(Z(k,1)- Zpre_pf(i,k))^2); % R[---待确定---]
    end
    weight(:,k)=weight(:,k)./sum(weight(:,k));
    
    % 粒子 重新采样
    outIndex = multinomialR(weight(:,k));
    
    % 产生粒子滤波后的 所有粒子
    Xpf(:,k,:)= Xparticles(outIndex,k,:);
    
    we = 1;
end

% 最终目标的计算

Xmean_pf=mean(Xpf);  
bins=20;
Xmap_pf=zeros(T,1);
for k=1:T
    [p,pos]=hist(Xpf(:,k,1),bins);
    map=find(p==max(p));
    Xmap_pf(k,1)=pos(map(1));  
end
for k=1:T
    Xstd_pf(1,k)=std(Xpf(:,k)-X(k,1)); 
end


% 程序结束提醒
disp('Done');