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
R = 1.5;
t = linspace(pi,2*pi,T);
X = [(2.5 + R*cos(t))',(3 + R*sin(t))']; % 位置完成
for i=1:T-1 % 速度只有T为1:T-1有
    X(i,3) = (X(i+1,1) - X(i,1) ) / dT; % x方向的速度
    X(i,4) = (X(i+1,2) - X(i,2) ) / dT; % y方向的速度
end % 速度完成

% 3. 粒子滤波
% 3.1. 各项参数
numSamples = 100;
Xpf=zeros(numSamples,T,4); % 行代表某一个粒子，列代表某一个时刻，值为这个粒子在当前时刻的 粒子滤波后（除初始化）的 状态
Xparticles=zeros(numSamples,T,4); % 行代表某一个粒子，列代表某一个时刻，值为这个粒子在当前时刻的 粒子滤波前的 状态
Zpre_pf=zeros(numSamples,T); % 行代表某一个粒子，列代表某一个时刻，值为这个粒子在当前时刻的 观测值
weight=zeros(numSamples,T); % 行代表某一个粒子，列代表某一个时刻，值为这个粒子在当前时刻的 权重
QQQ = 0.01; % 高斯滤波的权值的平方[---待确认---]

%%%%%%%%%%%%%%%%%%%%PSO%%%%%%%%%%%%%%%%%%%%
problem.CostFunction = @(R,x,y)CostFunction(R,x,y);
problem.nVar =5;      % Number of Unknown (Decision) Variables
problem.VarMin = 0;   % Lower Bound of Decision Variables最小位置[---待确定---]
problem.VarMax = 5;   % Upper Bound of Decision Variables最大位置[---待确定---]
% PSO程序的数值赋予
params.MaxIt = 20;          % Maximum Number of Iterations迭代次数[---待确定---]
params.nPop = numSamples;   % Population Size (Swarm Size)粒子数目
params.c1 = 2;              % Personal Acceleration Coefficient学习因子[---待确定---]
params.c2 = 2;              % Social Acceleration Coefficient学习因子[---待确定---]

CostFunction = problem.CostFunction;
nVar = problem.nVar;        % Number of Unknown (Decision) Variables
VarSize = [nVar 1];         % Matrix Size of Decision Variables
VarMin = problem.VarMin;	% Lower Bound of Decision Variables
VarMax = problem.VarMax;    % Upper Bound of Decision Variables
% PSO参数
MaxIt = params.MaxIt;   % Maximum Number of Iterations最大迭代次数
nPop = params.nPop;     % Population Size (Swarm Size)粒子的数目
c1 = params.c1;         % Personal Acceleration Coefficient
c2 = params.c2;         % Social Acceleration Coefficient
MaxVelocity = 1; % 最大速度[---待确定---]
MinVelocity = - MaxVelocity; % [---待确定---]

%初始化粒子的状态：
GlobalBestCost = 0;
GlobalBestPosition = [0,0];
particleBestCost = zeros(numSamples);
particleBestPosition = zeros(numSamples,2);

% PSO变量
iparticles = zeros(T,2);
particlePosition = zeros(numSamples,2);
particleVelocity = zeros(numSamples,2);
particleCost = zeros(numSamples);
%%%%%%%%%%%%%%%%%%%%PSO%%%%%%%%%%%%%%%%%%%%

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
        QQ=0.01; % 网
        net=sqrt(QQ)*randn(4,1); % 网
        temp = zeros(1,4); % temp为上一时刻 当前粒子的 状态
        temp(1) = Xpf(i,k-1,1);
        temp(2) = Xpf(i,k-1,2);
        temp(3) = Xpf(i,k-1,3);
        temp(4) = Xpf(i,k-1,4);
        Xparticles(i,k,:) = langevin(temp' )' + net';% 网
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
    
    %%%%%%%%%%%%%%%%%%%%接下来是PSO优化后的计算权重的步骤%%%%%%%%%%%
    % 对粒子滤波产生的状态估计进行处理
    for i=1:numSamples
        iparticles(k,1) = iparticles(k,1) + Xparticles(i,k,1) * weight(i,k); % x  
        iparticles(k,2) = iparticles(k,2) + Xparticles(i,k,2) * weight(i,k); % y
    end
    
    %初始化粒子的状态
    for i=1:nPop
        particlePosition(i,:) = Xparticles(i,k,1:2); % 将重采样后粒子的实际值带入到particle.position中去
        particleVelocity(i,:) = 2 * rand(1,2); % 粒子所有的初始速度为 均匀分布随机数*2
        particleCost(i) = CostFunction(nPop,particlePosition(i,:),iparticles(k,:));%对这些测量值计算权重
        % Update the Personal Best
        particleBestCost(i) = particleCost(i);
        particleBestPosition(i,:) = particlePosition(i,:);
        % Update Global Best
        if particleBestCost(i) > GlobalBestCost
            GlobalBestPosition = particleBestPosition(i,:);
            GlobalBestCost = particleBestCost(i);
        end
    end
    
    % 迭代
    for it=1:MaxIt
        for i=1:nPop
            %这里是由同一时刻的nPop个粒子的真实状态，套入函数中来计算nPop个粒子的测量值
            w1=0.5*exp(1-it/MaxIt);
            particleVelocity(i,1) = w1*particleVelocity(i,1) ... % x
                + c1*rand().*(particleBestPosition(i,1) - particlePosition(i,1)) ...
                + c2*rand().*(GlobalBestPosition(1) - particlePosition(i,1));
            particleVelocity(i,2) = w1*particleVelocity(i,2) ... % y
                + c1*rand().*(particleBestPosition(i,2) - particlePosition(i,2)) ...
                + c2*rand().*(GlobalBestPosition(2) - particlePosition(i,2));

            % Apply Velocity Limits确保粒子的速度在速度限制范围之内
            particleVelocity(i,1) = max(particleVelocity(i,1), MinVelocity);
            particleVelocity(i,1) = min(particleVelocity(i,1), MaxVelocity);
            particleVelocity(i,2) = max(particleVelocity(i,2), MinVelocity);
            particleVelocity(i,2) = min(particleVelocity(i,2), MaxVelocity);

            % Update Position
            particlePosition(i,1) = particlePosition(i,1) + particleVelocity(i,1);
            particlePosition(i,2) = particlePosition(i,2) + particleVelocity(i,2);

            % Apply Lower and Upper Bound Limits确保粒子的位置在范围之内
            particlePosition(i,1) = max(particlePosition(i,1), VarMin);
            particlePosition(i,1) = min(particlePosition(i,1), VarMax);
            particlePosition(i,2) = max(particlePosition(i,2), VarMin);
            particlePosition(i,2) = min(particlePosition(i,2), VarMax);

            % Evaluation%每个粒子经过运动后的新的权值
            particleCost(i) = CostFunction(nPop,particlePosition(i,:),iparticles(k,:));
            if  particleCost(i) > particleBestCost(i)
                particleBestPosition(i,:)= particlePosition(i,:);
                particleBestCost(i) = particleCost(i);
                
                % Update Global Best
                if particleBestCost(i) > GlobalBestCost
                    GlobalBestCost = particleBestCost(i);
                    GlobalBestPosition = particleBestPosition(i,:);
                end
            end
        end
    end
    
    % 带入Xpf(:,k,:)并重新计算权重
    for i=1:nPop
        temp = particlePosition(i,:);
        Xpf(i,k,1:2) = temp;
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
    
    % 重新置0
    GlobalBestCost = 0;
    GlobalBestPosition = [0,0];
    particleBestCost = zeros(numSamples);
    particleBestPosition = zeros(numSamples,2);
    
    weight(:,k)=weight(:,k)./sum(weight(:,k));

    % 粒子 重新采样
    outIndex = multinomialR(weight(:,k));
    
    % 产生粒子滤波后的 所有粒子
    Xpf(:,k,:)= Xpf(outIndex,k,:);
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
RMSE = sqrt(1/T * sum);

figure(12);clf;
k=1:1:T;
plot(k,X(:,1),'b',k,Xmean_x_pf,'r',k,Xmap_x_pf,'g'); 
legend('系统真实状态值','后验均值估计','最大后验概率估计');
xlabel('次数','fontsize',15);
ylabel('X状态估计','fontsize',15);
saveas(12,'./jpg/X估计值与真值.jpg'); % 保存

figure(22);clf;  
k=1:1:T;
plot(k,X(:,2),'b',k,Xmean_y_pf,'r',k,Xmap_y_pf,'g'); 
legend('系统真实状态值','后验均值估计','最大后验概率估计');
xlabel('次数','fontsize',15);
ylabel('Y状态估计','fontsize',15);
saveas(22,'./jpg/Y估计值与真值.jpg'); % 保存

figure(16);
k=1:1:T;
plot(k,Xdiff_pf,'-');
xlabel('次数');ylabel('状态估计误差');
titleRMSE = strcat('RMSE = ',num2str(RMSE) );
title(titleRMSE);
axis([0,T,0,5] );
saveas(16,'./jpg/RMSE.jpg'); % 保存

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 程序结束提醒
disp('Done');