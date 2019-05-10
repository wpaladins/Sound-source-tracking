function [ weight ] = particle_weight_generator(gccResult,Nd,fs,tdoaT)
% 本函数 为 求粒子权重的函数

% 延迟点数
DS = floor(tdoaT * fs);

% GCC
RSKDS = 1000 * gccResult(Nd + DS);

disp('RSKDS');
disp(RSKDS);
xi0 = 0.00025;
gamma = 3;

weight = (max(RSKDS,xi0)) ^ gamma;

end

