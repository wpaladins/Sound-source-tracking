clc;
clear;

% 1. load
track = 2; % track = 1圆形轨迹；track = 2三角轨迹；track = <其他>直线轨迹
loadlocation = './mat/';
if track == 1
    loadlocation = strcat(loadlocation,'semicircle/');
elseif track == 2
    loadlocation = strcat(loadlocation,'triangle/');
else
    loadlocation = strcat(loadlocation,'straightLine/');
end
loadlocation_pf = strcat(loadlocation,'pf.mat');
loadlocation_psopf = strcat(loadlocation,'psopf.mat');
% 1.1. pf
load(loadlocation_pf);
% 1.2. psopf
load(loadlocation_psopf);
% 1.3. 麦克风位置
s1r1 = [2.2,0.5]; s1r2 = [2.8,0.5];
s2r1 = [2.2,4.5]; s2r2 = [2.8,4.5];
s3r1 = [0.5,2.2]; s3r2 = [0.5,2.8];
s4r1 = [4.5,2.2]; s4r2 = [4.5,2.8];
srx = [s1r1(1),s1r2(1),...
       s2r1(1),s2r2(1),...
       s3r1(1),s3r2(1),...
       s4r1(1),s4r2(1)];
sry = [s1r1(2),s1r2(2),...
       s2r1(2),s2r2(2),...
       s3r1(2),s3r2(2),...
       s4r1(2),s4r2(2)];
% 1.4. 运动轨迹
T = 50; % 说话人位置改变次数，帧数
X = zeros(50,4);
if track == 1
    R = 1.5;
    t = linspace(pi,2*pi,T);
    X(:,1:2) = [(2.5 + R*cos(t))',(3 + R*sin(t))'];
elseif track == 2
    T1 = round(T/2);
    a = [1,2]; b = [2.5,1]; c = [4,2];
    t1x = linspace(a(1),b(1),T1);
    t1y = linspace(a(2),b(2),T1);
    t2x = linspace(b(1),c(1),T - T1 + 1);
    t2y = linspace(b(2),c(2),T - T1 + 1);
    X(1:T1,1:2) = [t1x',t1y'];
    X(T1:T,1:2) = [t2x',t2y'];
else
    t = linspace(0.5,4.5,T);
    X(:,1:2) = [t',t'];
end

% 2. 画图
% 2.1. 时间
figure(1);
hold on;box on;
p1 = plot(1:T,Tpf,'-go','lineWidth',1); % SIR-GCF-PF
p2 = plot(1:T,Tpsopf,'-r^','lineWidth',1); % APSO-GCF-PF
legend([p1,p2],'SIR-GCF-PF time','APSO-GCF-PF time');
saveas(1,'./mat/jpg/时间.jpg'); % 保存
T_mean_pf = mean(Tpf);
T_mean_psopf = mean(Tpsopf);
disp('Mean Tpf');
disp(T_mean_pf);
disp('Mean Tpsopf')
disp(T_mean_psopf);


% 2.2. X真值与估计值
figure(2);clf;
hold on;box on;
p1 = plot(1:T,X(:,1)','-b.','lineWidth',1);
p2 = plot(1:T,Xmean_x_pf,'-go','lineWidth',1);
p3 = plot(1:T,Xmean_x_psopf,'-r^','lineWidth',1);
legend([p1,p2,p3],'True X value','SIR-GCF-PF posterior mean estimate of X','APSO-GCF-PF posterior mean estimate of X');
xlabel('Time frame','fontsize',15);
ylabel('Value','fontsize',15);
axis([0 50 0 5]);
title('The comparison between True X value and Posterior mean estimate of X','fontsize',15);
saveas(2,'./mat/jpg/X估计值与真值.jpg'); % 保存

% 2.3. Y真值与估计值
figure(3);clf;
hold on;box on;
p1 = plot(1:T,X(:,2)','-b.','lineWidth',1);
p2 = plot(1:T,Xmean_y_pf,'-go','lineWidth',1);
p3 = plot(1:T,Xmean_y_psopf,'-r^','lineWidth',1);
legend([p1,p2,p3],'True Y value','SIR-GCF-PF posterior mean estimate of Y','APSO-GCF-PF posterior mean estimate of Y');
xlabel('Time frame','fontsize',15);
ylabel('Value','fontsize',15);
axis([0 50 0 5]);
title('The comparison between True Y value and Posterior mean estimate of Y','fontsize',15);
saveas(3,'./mat/jpg/Y估计值与真值.jpg'); % 保存

% 2.4. 偏差图
figure(4);
hold on;box on;
p1 = plot(1:T,Xdiff_pf,'-go','lineWidth',1);
p2 = plot(1:T,Xdiff_psopf,'-r^','lineWidth',1);
legend([p1,p2],'SIR-GCF-PF deviation','APSO-GCF-PF deviation');
xlabel('Time frame','fontsize',15);
ylabel('Deviation','fontsize',15);
title('Deviation of each time frame');
axis([0,T,0,5] );
saveas(4,'./mat/jpg/偏差图.jpg'); % 保存

% 2.5. 运动轨迹和麦克风位置图
figure(5);
hold on;box on;
p1 = plot(X(:,1),X(:,2),'b.');
p2 = plot(srx,sry,'ro');
legend([p1,p2],'Motion track','Microphone');
axis([0 5 0 5]);
saveas(5,'./mat/jpg/运动轨迹和麦克风位置图.jpg'); % 保存

% 3. 显示RMSE
disp('SIR-GCF-PF RMSE');
disp(RMSE_pf);
disp('APSO-GCF-PF RMSE');
disp(RMSE_psopf);