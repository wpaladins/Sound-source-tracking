function [xn_frams] = framing(xn,fs,fram_time,fram_step_time,win)
%函数用来对一个信号进行分帧    
%                   [xn_frams] = framing(xn,fs,fram_time,fram_step_time,win)
%输入参数：
%                   xn:待分帧的序列
%                   fs:采样率
%                   fram_time:帧时长
%                   fram_step_time:帧移动的时长
%                   win:窗函数，可选择海明窗或汉明窗
%输出参数：
%                   xn_frams:一个帧数为列，帧长度为行的矩阵，是对xn的分帧结果

[row, col]=size(xn);
if row>col
    xn=xn';
else
end

ts=1/fs;

l=length(xn);%序列的总长度
xn_time=l*ts;%#ok<*NASGU> %序列的总时长
fram_length=ceil(fram_time/ts);%每一帧的长度
fram_step_length=ceil(fram_step_time/ts);%帧移的长度
if win=='hanning'
    win=hanning(fram_length);
elseif win=='hamming'
    win=hamming(fram_length);
end
%帧数的计算公式为：序列点数-帧长度+帧移）/帧移
numOfframs=(l-fram_length+fram_step_length)/fram_step_length;%计算得到帧数
%考虑到序列总长度和分帧结果不匹配，补0以满足分帧要求
numOfframs=ceil(numOfframs);
%反算为了满足分帧序列应该增加的长度
l_added=(numOfframs*fram_step_length-fram_step_length+fram_length);%得到序列应该有的长度
l=l_added-l;
xn=[xn,zeros(1,l)];%补0
xn_time=ceil(l*ts);
xn_frams=zeros(fram_length,numOfframs);%建立存放结果的矩阵

%开始分帧
for k=1:numOfframs
    dn=(k-1)*fram_step_length+(1:fram_length);
    xn_frams(:,k)=xn(dn).*win';
end