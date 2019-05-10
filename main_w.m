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

% 2.3. 真实状态信息展示
% figure;
% plot(X(:,1),X(:,2),'b.');
% axis([0 5 0 5]);

% 3. 粒子滤波
% 3.1. 各项参数
numSamples = 100;
Xpf=zeros(numSamples,T,4); % 行代表某一个粒子，列代表某一个时刻，值为这个粒子在当前时刻的 粒子滤波后（除初始化）的 状态
Xparticles=zeros(numSamples,T,4); % 行代表某一个粒子，列代表某一个时刻，值为这个粒子在当前时刻的 粒子滤波前的 状态
Zpre_pf=zeros(numSamples,T); % 行代表某一个粒子，列代表某一个时刻，值为这个粒子在当前时刻的 观测值
weight=zeros(numSamples,T); % 行代表某一个粒子，列代表某一个时刻，值为这个粒子在当前时刻的 权重
QQQ = 0.01; % 高斯滤波的权值的平方[---待确认---]

% 粒子初始化[---待确认---]这里使用了真实值
Xpf(:,1,:)=X(1,:)+sqrt(QQQ)*randn(numSamples,4); % 初始粒子状态，使用高斯滤波对真实状态处理产生

% 画图显示粒子分布情况
figure(1);
plot(Xpf(:,1,1),Xpf(:,1,2),'g.',X(1,1),X(1,2),'r.');
axis([0 5 0 5]);
saveas(1,'./jpg/1.jpg');

% 粒子滤波核心循环
for k=2:T
    % 得到gccResult,Nd
    [h1,h2] = rir_example(X(k,[1,2]));
    conv1 = conv(rawWav,h1);
    conv2 = conv(rawWav,h2);
    [gccResult,Nd] = gcc_phat_w(conv1,conv2);
    disp('真实位置');
    disp(X(k,[1,2]));

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
        temp = zeros(1,2);
        temp(1) = Xparticles(i,k,1);
        temp(2) = Xparticles(i,k,2);
        disp('粒子位置');
        disp(temp);
        tdoaT = tdoaT_generator(temp,[1.2,0.5],[1.8,0.5]);
        % 更新粒子权重
        weight(i,k) = particle_weight_generator(gccResult,Nd,fs,tdoaT);
    end
    
    figure(3) % 权重图
    plot3(Xparticles(:,k,1),Xparticles(:,k,2),weight(:,k));
    axis([0 5 0 5]);
    jpg = strcat('./jpg/p',num2str(k));
    jpg = strcat(jpg,'.jpg');
    saveas(3,jpg);
    
    weight(:,k)=weight(:,k)./sum(weight(:,k));

    % 粒子 重新采样
    outIndex = multinomialR(weight(:,k));
    
    % 产生粒子滤波后的 所有粒子
    Xpf(:,k,:)= Xparticles(outIndex,k,:);
    
    we = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 最终目标的计算
Xmean_x_pf=mean(Xpf(:,:,1));
Xmean_y_pf=mean(Xpf(:,:,2));

% Xmean_pf=mean(Xpf);

% 
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
for k=1:T
    Xstd_x_pf(1,k)=std(Xpf(:,k,1)-X(k,1)); 
end
for k=1:T
    Xstd_y_pf(1,k)=std(Xpf(:,k,2)-X(k,1)); 
end

%--20190508 figure(11);clf;
%--20190508 subplot(221);
%--20190508 plot(v);
%--20190508 xlabel('时间');
%--20190508 ylabel('测量噪声','fontsize',15);
%--20190508 subplot(222);
%--20190508 plot(w);    
%--20190508 xlabel('时间');
%--20190508 ylabel('过程噪声','fontsize',15);
%--20190508 subplot(223);
%--20190508 plot(X);   
%--20190508 xlabel('时间','fontsize',15);
%--20190508 ylabel('状态X','fontsize',15);
%--20190508 subplot(224);
%--20190508 plot(Z);   
%--20190508 xlabel('时间','fontsize',15);
%--20190508 ylabel('观测Z','fontsize',15);


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


%--20190508 figure(13);
% subplot(121);
%--20190508 plot(Xmean_x_pf,X(:,1),'+');
%--20190508 xlabel('X后验均值估计','fontsize',15);
%--20190508 ylabel('X真值','fontsize',15)
%--20190508 hold on;
%--20190508 c=0:1:5;
%--20190508 plot(c,c,'r');
%--20190508 axis([0 5 0 5]);
%--20190508 hold off;

%--20190508 subplot(122);  
%--20190508 plot(Xmap_pf,X,'+')
%--20190508 ylabel('真值','fontsize',15)
%--20190508 xlabel('MAP估计','fontsize',15)
%--20190508 hold on;
%--20190508 c=-25:1:25;
%--20190508 plot(c,c,'r');  
%--20190508 axis([-25 25 -25 25]);
%--20190508 hold off;

%--20190508 figure(23);
% subplot(121);
%--20190508 plot(Xmean_y_pf,X(:,2),'+');
%--20190508 xlabel('Y后验均值估计','fontsize',15);
%--20190508 ylabel('Y真值','fontsize',15)
%--20190508 hold on;
%--20190508 c=0:1:5;
%--20190508 plot(c,c,'r');
%--20190508 axis([0 5 0 5]);
%--20190508 hold off;

%--20190508 subplot(122);  
%--20190508 plot(Xmap_pf,X,'+')
%--20190508 ylabel('真值','fontsize',15)
%--20190508 xlabel('MAP估计','fontsize',15)
%--20190508 hold on;
%--20190508 c=-25:1:25;
%--20190508 plot(c,c,'r');  
%--20190508 axis([-25 25 -25 25]);
%--20190508 hold off;
 
%--20190508 domain=zeros(numSamples,1);
%--20190508 range=zeros(numSamples,1);
%--20190508 bins=10;
%--20190508 support=[-20:1:20];


%--20190508 figure(14);hold on; 
%--20190508 xlabel('样本空间','fontsize',15);
%--20190508 ylabel('时间','fontsize',15);
%--20190508 zlabel('后验密度','fontsize',15);
%--20190508 vect=[0 1];
%--20190508 caxis(vect);
%--20190508 for k=1:T
  
%--20190508     [range,domain]=hist(Xpf(:,k),support);
   
%--20190508     waterfall(domain,k,range);
%--20190508 end
%--20190508 axis([-20 20 0 T 0 100]);
 
%--20190508 figure(15);
%--20190508 hold on; box on;
%--20190508 xlabel('样本空间','fontsize',15);
%--20190508 ylabel('后验密度','fontsize',15); 
%--20190508 k=30;   
%--20190508 [range,domain]=hist(Xpf(:,k),support);
%--20190508 plot(domain,range);
 
%--20190508 XXX=[X(k,1),X(k,1)];
%--20190508 YYY=[0,max(range)+10]
%--20190508 line(XXX,YYY,'Color','r');
%--20190508 axis([min(domain) max(domain) 0 max(range)+10]);
 
figure(16);
k=1:1:T;
plot(k,Xstd_x_pf,'-');
xlabel('次数');ylabel('X状态估计误差标准差');
axis([0,T,0,1]);
saveas(16,'./jpg/X状态估计误差标准差.jpg'); % 保存

figure(26);
k=1:1:T;
plot(k,Xstd_y_pf,'-');
xlabel('次数');ylabel('Y状态估计误差标准差');
axis([0,T,0,1]);
saveas(26,'./jpg/Y状态估计误差标准差.jpg'); % 保存

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 程序结束提醒
disp('Done');