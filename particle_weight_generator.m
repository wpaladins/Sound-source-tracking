function [weight] = particle_weight_generator(s1_gccResult,s1_Nd,...
                                                s2_gccResult,s2_Nd,...
                                                s3_gccResult,s3_Nd,...
                                                s4_gccResult,s4_Nd,...
                                                fs,...
                                                s1_tdoaT,...
                                                s2_tdoaT,...
                                                s3_tdoaT,...
                                                s4_tdoaT)
% 本函数 为 求粒子权重的函数

% 延迟点数
s1_DS = floor(s1_tdoaT * fs * 4);
s2_DS = floor(s2_tdoaT * fs * 4);
s3_DS = floor(s3_tdoaT * fs * 4);
s4_DS = floor(s4_tdoaT * fs * 4);

% GCC
s1_RSKDS = 1000 * s1_gccResult(s1_Nd + s1_DS);
s2_RSKDS = 1000 * s2_gccResult(s2_Nd + s2_DS);
s3_RSKDS = 1000 * s3_gccResult(s3_Nd + s3_DS);
s4_RSKDS = 1000 * s4_gccResult(s4_Nd + s4_DS);
RSKDS = (s1_RSKDS + s2_RSKDS + s3_RSKDS + s4_RSKDS)/4;

xi0 = 0;
gamma = 4;

weight = (max(RSKDS,xi0)) ^ gamma;

end

