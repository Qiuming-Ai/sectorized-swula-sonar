function BaselineDataMaker(velocity, snr_level, ping_num)
% BaselineDataMaker - Generate baseline data and save to HDF5 file
%
% Inputs:
%   velocity   - Velocity parameter (scalar)
%   snr_level  - SNR level in dB (scalar)
%   ping_num   - Number of pings (scalar)
%
% Output:
%   Creates an HDF5 file named 'BaselineData_YYYYMMDD_HHMMSS.h5'
%
% clear all; close all;
legacy_sonarsim_setup();
%% 基本常量
c0 = 1500;
fs = 2e6;
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
decimate_factor = 16;
%% Define phantom
% Define some points in a phantom for the simulation
point_position(1,:) = [0 0 498e-2];
point_position(end+1,:) = [0 0 500e-2];
point_position(end+1,:) = [0 0 504e-2];
point_position(end+1,:) = [0 0 512e-2];
point_position(end+1,:) = [0 0 528e-2];
point_position(end+1,:) = [0 0 550e-2];
point_position(end+1,:) = [512e-2*sind(2) 0 512e-2*cosd(2)];
point_position(end+1,:) = [512e-2*sind(-2) 0 512e-2*cosd(-2)];
point_position(end+1,:) = [512e-2*sind(6) 0 512e-2*cosd(6)];
point_position(end+1,:) = [512e-2*sind(-6) 0 512e-2*cosd(-6)];
point_position(end+1,:) = [512e-2*sind(12) 0 512e-2*cosd(12)];
point_position(end+1,:) = [512e-2*sind(-12) 0 512e-2*cosd(-12)];
point_position(end+1,:) = [512e-2*sind(20) 0 512e-2*cosd(20)];
point_position(end+1,:) = [512e-2*sind(-20) 0 512e-2*cosd(-20)];
point_position(end+1,:) = [512e-2*sind(36) 0 512e-2*cosd(36)];
point_position(end+1,:) = [512e-2*sind(-36) 0 512e-2*cosd(-36)];
point_position(end+1,:) = [512e-2*sind(42) 0 512e-2*cosd(42)];
point_position(end+1,:) = [512e-2*sind(-42) 0 512e-2*cosd(-42)];
point_position(end+1,:) = [512e-2*sind(55) 0 512e-2*cosd(55)];
point_position(end+1,:) = [512e-2*sind(-55) 0 512e-2*cosd(-55)];
% Set point amplitudes
point_amplitudes = ones(size(point_position,1),1);

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
[v, t] = sim_rx_from_scatterers_simple(point_position,point_amplitudes,lightPos,rx_xyz,f0,c0,fs,excitation,opts);
v = (v-min(v(:)))./(max(v(:))-min(v(:)));
v = (v-0.5)*2;
win = hamming(probe.N);
filt = SignalFilter('SampleRate',fs, ...
    'Passband', BW/2, ...
    'StopbandAtten', 60, ...
    'Type', 'lowpass', ...
    'SidebandMode', 'single', ...
    'OutputSameLength', true);
MF_mix = MF.*(exp(1j*2*pi*f0*[0:length(MF)-1]/fs)');
MF_fir = filt.process(MF_mix);
MF_deci =  resample(MF_fir,1,decimate_factor);
fs_deci = fs/decimate_factor;

% 生成时间戳
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

% 初始化 attributes
attributes.array_type = 'Baseline';
attributes.signal_type = 'Baseband';
attributes.signal_win = 'hamming';
attributes.bandwidth = BW;
attributes.sampling_frequency = fs;
attributes.center_frequency = f0;
attributes.decimate_factor = decimate_factor;
attributes.sector_num = 1;
attributes.match_filter_data = MF_deci;
attributes.receive_array_num = probe.N;
attributes.receive_array_position = rx_xyz;
attributes.sample_delay = t;
attributes.pulse_duration = pulse_len;
attributes.sound_velocity = c0;
attributes.ping_num = ping_num;
attributes.velocity = velocity;
attributes.snr_level = snr_level;
attributes.timestamp = timestamp;
attributes.scan_angle = (-60:1:60);
% 预分配 raw_data 存储
raw_data_cell = cell(ping_num, 1);

for i = 1:ping_num
    v_calc = resample(v,round((c0+2*velocity*2*(rand-0.5))*10),c0*10);
    v_calc = awgn(v_calc,snr_level,'measured');
    %% 下变频
    v_mix = v_calc.*repmat(exp(1j*2*pi*f0*[0:length(v_calc(:,1))-1]/fs)',[1 probe.N]);
    %% FIR滤波器
    y_fir = filt.process(v_mix);
    
    %% Decimation
    raw_data = resample(y_fir,1,decimate_factor);
    raw_data_cell{i} = raw_data;
end

%% 写入 HDF5 文件
filename = sprintf('BaselineData_%s.h5', timestamp);

% 创建 raw_data 组
group_name = '/raw_data';
assert(~isfile(filename), 'Refusing to overwrite existing generated HDF5: %s', filename);

% 写入每个 ping 的数据集
for i = 1:ping_num
    dataset_name = sprintf('%s/ping_%d', group_name, i);
    ping_data = raw_data_cell{i};
    
    % 检查是否为复数数据
    if isreal(ping_data)
        % 实数数据：直接写入
        h5create(filename, dataset_name, size(ping_data), 'Datatype', 'single');
        h5write(filename, dataset_name, single(ping_data));
        h5writeatt(filename, dataset_name, 'complex', 0);
    else
        % 复数数据：分别存储实部和虚部
        h5create(filename, [dataset_name '/real'], size(ping_data), 'Datatype', 'single');
        h5write(filename, [dataset_name '/real'], single(real(ping_data)));
        h5create(filename, [dataset_name '/imag'], size(ping_data), 'Datatype', 'single');
        h5write(filename, [dataset_name '/imag'], single(imag(ping_data)));
        h5writeatt(filename, dataset_name, 'complex', 1);
    end
end

% 写入 attributes 到 raw_data/.attributes 组
% 创建 .attributes 子组
attrs_group = [group_name '/.attributes'];
% 使用低级 HDF5 API 创建组
fid = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
try
    gid = H5G.create(fid, attrs_group, 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
    H5G.close(gid);
catch ME
    % 如果组已存在，则打开它
    if contains(ME.message, 'already exists') || contains(ME.identifier, 'exists')
        try
            gid = H5G.open(fid, attrs_group);
            H5G.close(gid);
        catch ME2
            H5F.close(fid);
            rethrow(ME2);
        end
    else
        H5F.close(fid);
        rethrow(ME);
    end
end
H5F.close(fid);

% 处理数值型 attributes（存储为 .attributes 组的属性）
if isfield(attributes, 'sampling_frequency')
    h5writeatt(filename, attrs_group, 'sampling_frequency', attributes.sampling_frequency);
end
if isfield(attributes, 'center_frequency')
    h5writeatt(filename, attrs_group, 'center_frequency', attributes.center_frequency);
end
if isfield(attributes, 'decimate_factor')
    h5writeatt(filename, attrs_group, 'decimate_factor', attributes.decimate_factor);
end
if isfield(attributes, 'sector_num')
    h5writeatt(filename, attrs_group, 'sector_num', attributes.sector_num);
end
if isfield(attributes, 'receive_array_num')
    h5writeatt(filename, attrs_group, 'receive_array_num', attributes.receive_array_num);
end
if isfield(attributes, 'sound_velocity')
    h5writeatt(filename, attrs_group, 'sound_velocity', attributes.sound_velocity);
end
if isfield(attributes, 'ping_num')
    h5writeatt(filename, attrs_group, 'ping_num', attributes.ping_num);
end
if isfield(attributes, 'velocity')
    h5writeatt(filename, attrs_group, 'velocity', attributes.velocity);
end
if isfield(attributes, 'snr_level')
    h5writeatt(filename, attrs_group, 'snr_level', attributes.snr_level);
end
if isfield(attributes, 'bandwidth')
    h5writeatt(filename, attrs_group, 'bandwidth', attributes.bandwidth);
end
if isfield(attributes, 'pulse_duration')
    h5writeatt(filename, attrs_group, 'pulse_duration', attributes.pulse_duration);
end
if isfield(attributes, 'sample_delay')
    h5writeatt(filename, attrs_group, 'sample_delay', attributes.sample_delay);
end
% 处理字符串型 attributes
if isfield(attributes, 'array_type')
    h5writeatt(filename, attrs_group, 'array_type', attributes.array_type);
end
if isfield(attributes, 'signal_type')
    h5writeatt(filename, attrs_group, 'signal_type', attributes.signal_type);
end
if isfield(attributes, 'signal_win')
    h5writeatt(filename, attrs_group, 'signal_win', attributes.signal_win);
end
if isfield(attributes, 'timestamp')
    h5writeatt(filename, attrs_group, 'timestamp', attributes.timestamp);
end

% 处理数组型 attributes（作为数据集存储在 .attributes 组下）
if isfield(attributes, 'receive_array_position')
    h5create(filename, [attrs_group '/receive_array_position'], size(attributes.receive_array_position), 'Datatype', 'double');
    h5write(filename, [attrs_group '/receive_array_position'], attributes.receive_array_position);
end

if isfield(attributes, 'scan_angle')
    h5create(filename, [attrs_group '/scan_angle'], size(attributes.scan_angle), 'Datatype', 'double');
    h5write(filename, [attrs_group '/scan_angle'], attributes.scan_angle);
end

if isfield(attributes, 'match_filter_data')
    match_filter_data = attributes.match_filter_data;
    if isreal(match_filter_data)
        % 实数数据：直接写入
        h5create(filename, [attrs_group '/match_filter_data'], size(match_filter_data), 'Datatype', 'single');
        h5write(filename, [attrs_group '/match_filter_data'], single(match_filter_data));
        h5writeatt(filename, [attrs_group '/match_filter_data'], 'complex', 0);
    else
        % 复数数据：分别存储实部和虚部
        h5create(filename, [attrs_group '/match_filter_data/real'], size(match_filter_data), 'Datatype', 'single');
        h5write(filename, [attrs_group '/match_filter_data/real'], single(real(match_filter_data)));
        h5create(filename, [attrs_group '/match_filter_data/imag'], size(match_filter_data), 'Datatype', 'single');
        h5write(filename, [attrs_group '/match_filter_data/imag'], single(imag(match_filter_data)));
        h5writeatt(filename, [attrs_group '/match_filter_data'], 'complex', 1);
    end
end

fprintf('HDF5 file created: %s\n', filename);
end



