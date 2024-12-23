clear; clc; close all;

% ===================== 参数设置 =====================
delta_t = 1e-9;        % 时间分辨率 1 ns
t_min = 0;
t_max = 15e-6 - delta_t; % 脉冲持续时间约为15 µs
t = t_min:delta_t:t_max; 
N_chirp = length(t);    % 脉冲长度对应的采样点数（15000）
k = 20e6/N_chirp;       % 频率增量因子，用于线性调频

% 实际数据长度信息（与数据文件匹配）
% 数据为200条记录，每条记录长度为500000点
num_records = 200;
N_data = 500000; 

% 数据对应实际物理参数（根据第一版处理逻辑）
% 下述参数根据第一版代码逻辑获得
delta_T = 500e-6; % 每条记录之间的时间间隔为 500 µs (由N_data/fs=500µs，可不严格使用)
delta_X = 3e8*1e-9/2/1.46; % 单点对应的距离分辨率
T_axis = delta_T*(1:num_records); 
X_axis = delta_X*(1:N_data);

% ===================== 加载数据 =====================
load('data_polX.mat', 'data_polX_segement');
load('data_polY.mat', 'data_polY_segement');
i_x = data_polX_segement; % 大小为 200×500000
i_y = data_polY_segement; % 大小为 200×500000

% ===================== 生成与第一版相同的线性调频脉冲与匹配滤波器 =====================
% 第一版中的fs定义：fs = 300e6 + k*t/delta_t
% 注：t/delta_t = t * 1e9，即将t从s转为ns，再乘以k
f_start = 300e6;         % 起始频率300 MHz
fs_inst = f_start + k*(t/delta_t); 
st = exp(1j*2*pi*fs_inst.*t);
ht = conj(fliplr(st));   % 匹配滤波器（时域反转共轭）

% ===================== 数据处理与相位提取 =====================
% 对每条记录的X、Y数据进行Hilbert变换后与ht卷积，实现脉冲压缩
i_x_hilbert = zeros(num_records, N_data);
i_y_hilbert = zeros(num_records, N_data);
for i = 1:num_records
    ix_h = hilbert(i_x(i,:));
    iy_h = hilbert(i_y(i,:));
    % 卷积匹配滤波（'same'保证输出与输入长度一致）
    i_x_hilbert(i,:) = conv(ix_h, ht, 'same');
    i_y_hilbert(i,:) = conv(iy_h, ht, 'same');
end

% 偏振合成
I = i_x_hilbert + i_y_hilbert;

% 将第一条记录作为参考，对后续记录进行相位归一化
% 与第一版逻辑一致：对倒序索引进行处理，使其他行以第一行为参考进行相对相位调整
for i = 1:num_records
    I(num_records - i + 1,:) = I(num_records - i + 1,:) .* conj(I(1,:));
end

% 距离相关相位差分处理（与第一版相同）
delta_x = round(3e8/2/1.46/(20e6)/delta_X);
for j = 1:(N_data - delta_x)
    I(:, N_data - j + 1) = I(:, N_data - j + 1) .* conj(I(:, N_data - j + 1 - delta_x));
end

% 提取相位
phase = angle(I);

% ===================== 结果可视化 =====================
figure(1)
pcolor(X_axis, T_axis, phase);
shading interp; 
colorbar; colormap(jet);
xlabel('X/m'); ylabel('T/s');
xlim([26200,26600]);
title('相位分布图');

figure(2)
i = 255515; % 选定某一距离点
plot(T_axis, phase(:,i));
xlabel('time/s'); ylabel('phase/rad');
title('特定位置处相位随时间的变化');
