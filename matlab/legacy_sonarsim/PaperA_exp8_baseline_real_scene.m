clear all; close all;
legacy_sonarsim_setup();
%% 基本常量
c0 = 1500;
fs = 2e6;
interp_factor = 8;
dt = 1/fs;
f0 = 0.45e+06;      % Transducer center frequency [Hz]

lambda                  = c0/f0;           % Wavelength [m]
probe.element_height    = 15e-3;            % Height of element [m]
probe.pitch             = 0.5*lambda;        % probe.pitch [m]
kerf                    = 0.05*lambda;        % gap between elements [m]
probe.element_width     = probe.pitch-kerf;% Width of element [m]
% lens_el                 = 20e-3;           % position of the elevation focus
probe.N                 = 128;             % Number of elements
pulse_len=5e-3; 
lightPos = [0, 0, 0];
rx_xyz = [linspace(-probe.N*probe.pitch/2,probe.N*probe.pitch/2,probe.N)',zeros(probe.N,2)];
velocity = [0];
snr_level = [10];
%% Define phantom
% Define some points in a phantom for the simulation
%% ═══════════════════════ 用户参数 ═══════════════════════════════════════
stl_file    = legacy_model_path('ship.stl');    % optional exact local mesh
stl_scale   = 1;                     % 🔶 比例因子 0.5=缩小一半
stl_pos     = [0 40 40];                % 📌 放置位置 (x ,z ,y[depth])，注意 y 被用于 depth
stl_rot_deg = [90 0 -45];               % [roll pitch yaw]°

%% ═══════════════════════ ① STL 处理 ════════════════════════════════════
target1TR = stlread(stl_file);
stl_V = target1TR.Points;
%%
stl_V = stl_scale * stl_V;             % 🔶 缩放
r = stl_rot_deg;
Rxm = [1 0 0;0 cosd(r(1)) -sind(r(1));0 sind(r(1)) cosd(r(1))];
Rym = [cosd(r(2)) 0 sind(r(2));0 1 0;-sind(r(2)) 0 cosd(r(2))];
Rzm = [cosd(r(3)) -sind(r(3)) 0;sind(r(3)) cosd(r(3)) 0;0 0 1];
Rot  = Rzm*Rym*Rxm;
% 📌 变换并重排为 (x, y[depth], z)
pos_target1 = (Rot * stl_V.').' + stl_pos;  % 原为 Nx3 的 [x z y]
pos_target1 = pos_target1(:, [1 3 2]);       % 改为 [x y z]，其中 y 为 depth
target1TR = triangulation(target1TR.ConnectivityList,pos_target1);
clear pos_target1
%% ═══════════════════════ 用户参数 ═══════════════════════════════════════
stl_file    = legacy_model_path('merged_model.stl');
stl_scale   = 0.03;                     % 🔶 比例因子 0.5=缩小一半
stl_pos     = [0 15 8];                % 📌 放置位置 (x ,z ,y[depth])，注意 y 被用于 depth
stl_rot_deg = [0 0 90];               % [roll pitch yaw]°

%% ═══════════════════════ ① STL 处理 ════════════════════════════════════
target2TR = stlread(stl_file);
stl_V = target2TR.Points;
%%
stl_V = stl_scale * stl_V;             % 🔶 缩放
r = stl_rot_deg;
Rxm = [1 0 0;0 cosd(r(1)) -sind(r(1));0 sind(r(1)) cosd(r(1))];
Rym = [cosd(r(2)) 0 sind(r(2));0 1 0;-sind(r(2)) 0 cosd(r(2))];
Rzm = [cosd(r(3)) -sind(r(3)) 0;sind(r(3)) cosd(r(3)) 0;0 0 1];
Rot  = Rzm*Rym*Rxm;
% 📌 变换并重排为 (x, y[depth], z)
pos_target2 = (Rot * stl_V.').' + stl_pos;  % 原为 Nx3 的 [x z y]
pos_target2 = pos_target2(:, [1 3 2]);       % 改为 [x y z]，其中 y 为 depth
target2TR = triangulation(target2TR.ConnectivityList,pos_target2);
clear pos_target2

%% ═══════════════════════ ② 海底散射体生成（以 y 为深度） ════════════════
seabed_len  = 100;  seabed_wid = 200;    % 📌 z×x 尺寸   m
seabed_y0   = 48;  seabed_amp = 0;     % 📌 平均深度 / 起伏（depth 在 y 方向）
N_seabed_scatter = 1e5;                % 海底散射点数

Nx = round( sqrt(N_seabed_scatter * seabed_wid / seabed_len) );
Nz = round( N_seabed_scatter / Nx );

x_lin = linspace(-seabed_wid/2, seabed_wid/2, Nx);   % 横向 x
z_lin = linspace(0, seabed_len, Nz);                % 纵向 z（"前方"）
[Xg,Zg] = meshgrid(x_lin, z_lin);

Yg = seabed_y0 + seabed_amp * ( ...
      0.5*sin(2*pi*Xg/seabed_wid) + 0.5*sin(2*pi*Zg/seabed_len) ) + ...
      0.1*randn(size(Xg));  % 📌 y 方向作为深度，加入扰动

% 📌 散射点为 [x, y[depth], z]
pos_seabed = [Xg(:) Yg(:) Zg(:)];

% --- 生成 ConnectivityList ---
[Nz, Nx] = size(Xg);
idx  = reshape(1:Nz*Nx, Nz, Nx);
v1 = idx(1:end-1,1:end-1); v2 = idx(1:end-1,2:end);
v3 = idx(2:end,1:end-1);   v4 = idx(2:end,2:end);
ConnectivityList = [ [v1(:) v2(:) v4(:)]; ...
                     [v1(:) v4(:) v3(:)] ];
% --- 构造三角网格对象 / 可视化 ---
seabedTR = triangulation(ConnectivityList, pos_seabed);
clear ConnectivityList pos_seabed
%%
model1 = struct(...
    'Points', target1TR.Points, ...
    'ConnectivityList', target1TR.ConnectivityList, ...
    'reflectivity', 0.7, ...
    'scattering_coeff', 1, ...
    'center_frequency', f0, ...
    'transmittance', 0.1, ...
    'isHollow', false, ...
    'internal_scattering', 0.5);

model2 = struct(...
    'Points', target2TR.Points, ...
    'ConnectivityList', target2TR.ConnectivityList, ...
    'reflectivity', 0.5, ...
    'scattering_coeff', 0.5, ...
    'center_frequency', f0, ...
    'transmittance', 0.1, ...
    'isHollow', false, ...
    'internal_scattering', 0.5);

model3 = struct(...
    'Points', seabedTR.Points, ...
    'ConnectivityList', seabedTR.ConnectivityList, ...
    'reflectivity', 0.1, ...
    'scattering_coeff', 0.9, ...
    'center_frequency', f0, ...
    'transmittance', 0.8, ...
    'isHollow', true, ...
    'internal_scattering', 0);
clear target1TR target2TR seabedTR
% 输入参数
models = {model1,model3};
%% -------------------------- 完整参数 --------------------------
opts = struct('tx_az_deg', 0, 'tx_el_deg', 45, 'rx_el_deg', 45,...
              'rx_mode','sector', 'rx_sector_span_deg', 0, 'rx_sector_center_deg', 0);
[tx_dir, rx_dir] = fls_preview_and_dirs(models, lightPos, rx_xyz, opts);
params = struct();
params.f0              = f0;
params.c               = c0;      % 声速
params.rx_dir          = rx_dir;      % 每通道接收指向 (Nr x 3) 或 [1x3] 通用
params.tx_dir          = tx_dir;
params.max_bounces     = 1;           % 1/2/3
params.samples_per_tri = 2;           % 每个三角形采样数（此实现用质心；>1可扩展随机重心）
params.eps             = 1e-5;        % 交并误差
params.seed            = uint64(20250821);

params.tx_cos_power = -1;
params.rx_cos_power = -1;
% BRDF / 语义
params.brdf_cos_power   = 1.0;      % cos^m 镜面指数
params.rough_sigma      = lambda*16;     % m, ~ λ/8 粗糙度（根据中心频率估一下λ）
params.amp_type         = 0;        % 0幅度 1强度
params.reflect_is_power = 1;        % 模型 reflectivity 为功率语义时设 1

% 多跳附加衰减
params.order2_extra = 0.7;
params.order3_extra = 0.4;

% 可选
params.max_arrivals = 1e6;  % 大场景时提高上限
out = ...
    stl_raytrace_arrivals_cuda3_5_1(models, lightPos, lightPos, params);
clear model1 model2 model3 models params

point_position = out.pos;
point_amplitudes = out.amp;
% point_position = [model1.Points;model3.Points];
% point_amplitudes = ones(size(point_position,1),1);
% [tx_dir, rx_dir] = fls_preview_and_dirs(models, lightPos, rx_xyz, opts);
% scatter3(point_position(:,1),point_position(:,2),point_position(:,3))
%% 激励脉冲
pulse.fractional_bandwidth = 2/9;        % probe bandwidth [1]
pulse.center_frequency = f0;
BW=pulse.fractional_bandwidth*f0; 
impulse_response = 1;
t_p=0:1/fs:pulse_len;
excitation = chirp(t_p,f0-BW/2,pulse_len,f0+BW/2).*hamming(numel(t_p))';
one_way_ir = conv(impulse_response,excitation);
MF = conj(flipud(hilbert(excitation(:))));                  % FIR 形式
two_way_ir = conv(one_way_ir,MF);
lag = length(two_way_ir)/2+1;   
%%
opts = struct('round','round','precision','single');
opts.array.use_resample = true;                   % 若 true 尝试调用 MATLAB resample

% opts.array.density_scale_tx = 1;   % 间距×2（密度约减半）
% opts.array.density_scale_rx = 1;   % 间距×3
% opts.array.max_sub_per_elem_tx = 128; % 每个TX阵元的虚拟子阵上限
% opts.array.max_sub_per_elem_rx = 256; % 每个RX阵元上限
% opts.array.rx.shape  = 'disk';
% opts.array.rx.size   = [probe.width/2 probe.height];
% opts.array.rx.steer = [45*ones(probe.N,1) zeros(probe.N,2)];
tic
[v, t] = sim_rx_from_scatterers_perTX_cuda(point_position,point_amplitudes,lightPos,rx_xyz,f0,c0,fs,excitation,opts);
toc
v = (v-min(v(:)))./(max(v(:))-min(v(:)));
v = (v-0.5)*2;
figure;
win = hamming(probe.N);
filt = SignalFilter('SampleRate',fs, ...
    'Passband', BW/2, ...
    'StopbandAtten', 60, ...
    'Type', 'lowpass', ...
    'SidebandMode', 'single', ...
    'OutputSameLength', true);
MF_mix = MF.*(exp(1j*2*pi*f0*[0:length(MF)-1]/fs)');
%%
for k = 1:length(snr_level)
    for j = 1:length(velocity)
        v_calc = resample(v,c0+2*velocity(j),c0);
        v_calc = awgn(v_calc,snr_level(k),'measured');
        %% 下变频
        v_mix = v_calc.*repmat(exp(1j*2*pi*f0*[0:length(v_calc(:,1))-1]/fs)',[1 probe.N]);
        %% FIR滤波器
        y_fir = filt.process(v_mix);
        MF_fir = filt.process(MF_mix);
        %% Decimation
        y_deci = resample(y_fir,1,16);
        MF_deci =  resample(MF_fir,1,16);
        % for i = 1:probe.N
        % y_deci(:,i) = decimate(y_fir(:,i),8);
        % end
        % MF_deci = decimate(MF_fir,8);
        fs_deci = fs/16;
        %% Match Filter
        sig_mf = zeros(size(y_deci));
        tic
        for m = 1:probe.N
            tmp = conv(y_deci(:,m), MF_deci, 'same');
            sig_mf(:,m) = tmp*win(m);
        end
        toc
        % TVG
        tvg = tvg_factory(t,'c',c0,'mode','custom','freq_kHz',f0/1e3,'G0_dB',0,'K',10);
        % 假设采样率 fs=500 kHz，记录 100 ms：
        tau = (0:length(sig_mf)-1)/fs_deci;  % τ轴
        G_amp = tvg.Gamp(tau)';  % 用于幅度/包络补偿
        for m = 1:probe.N
            sig_mf(:,m) = sig_mf(:,m).*G_amp;
        end

        data = sig_mf./max(abs(sig_mf(:)));
        rx = data;
        cfg = struct;
        cfg.rx_xyz = rx_xyz;
        cfg.c = c0;
        cfg.t0 = t-pulse_len/2;                      % 用你的实际 t_start
        cfg.fs = fs_deci;
        cfg.angles = (-65:0.1:65) * pi/180;   % 扫描角（rad）
        cfg.ranges = (4.5:0.01:6).';         % 扫描深度（m）
        cfg.BW = 200e3;
        cfg.f0 = 450e3;
        % cfg.w = hamming(probe.N);
        % 假设 rx 为 Nsamp×64 的原始复值数据
        cfg.fd.nfft   = 16384;          % 或 64 / 256，越大分数时延越精细，但更慢
        cfg.fd.sign   = 1;
        cfg.cf.enable = false;
        cfg.cf.mode = 'cf';
        cfg.cf.gamma = 0.1;
        [beam,ranges_out] = das_plane_wave_fd(cfg, rx);
                % toc
        % 可视化（强度dB）
        % 假设你已用前面的 DAS 得到 img (Nr×Na)
        ax = subplot(length(snr_level),length(velocity),(j-1)*3+k);
        % sector_plot(cfg.angles, cfg.ranges, img, 30,"rect");
        beam = compensate_range_shift(beam, 1/((cfg.fd.nfft/cfg.fs)/c0), c0, 0.035,'true');
        sector_plot(cfg.angles,ranges_out, abs(beam), 40,"sector");
        % title('Plane-wave DAS (integer delays)');
        hold on;
        % scatter3(point_position(:,1),point_position(:,3),point_position(:,2),8,'r','filled');
        % scatter3(asin(point_position(:,1)./sqrt(point_position(:,1).^2+point_position(:,3).^2))/pi*180,sqrt(point_position(:,1).^2+point_position(:,3).^2),point_position(:,2),12,'r','filled');
         % ylim([4.8 5.7]);
        %          view(3)
        % xlim([-1 1])
        % ylim([4.9 5.06])
    end
end




