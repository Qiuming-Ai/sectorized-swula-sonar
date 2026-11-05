function CDMDataMaker(velocity, snr_level, ping_num)
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
BW=100e3;
fmax = f0+BW/2;
fmin = f0-BW/2;
lambda                  = c0/fmin;           % Wavelength [m]
interp_factor = 8;
Nrx = 32;             % Number of elements
Ntx = 12;
sector_num = 6;
lightPos = [0, 0, 0];
rx_xyz = [linspace(-(Nrx-1)*lambda,(Nrx-1)*lambda,Nrx)',zeros(Nrx,2)];
tx_xyz = [linspace(-(Ntx-1)*lambda/4,(Ntx-1)*lambda/4,Ntx)',zeros(Ntx,2)];

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
% 子带参数
angles_deg = -[-47.5,-26,-8.5,8.5,26,47.5]';              % 3 个方位角
fc   = 450e3;              % kHz 级别
BW   =100e3;

[exc_nt,exc_ori, delays_samp] = gen_multicode_cazac_excitation( ...
    tx_xyz, angles_deg, fc, BW, interp_factor*fs, c0,'amp',[1.2 1.2 1 1 1.2 1.2]);
 MF = conj(flipud(hilbert(exc_ori)));                  % FIR 形式
pulse_len=length(MF(:,1))/(fs*interp_factor); 
% 假设你已有：
%%
opts = struct('round','round','precision','single');
opts.array.use_resample = true;                   % 若 true 尝试调用 MATLAB resample

[v, t] = sim_rx_from_scatterers_perTX(point_position,point_amplitudes,tx_xyz,rx_xyz,f0,c0,interp_factor*fs,exc_nt,opts);
v = (v-min(v(:)))./(max(v(:))-min(v(:)));
v = (v-0.5)*2;
MF = resample(MF,1,interp_factor);
v =  resample(v,1,interp_factor);
win = hamming(Nrx);
filt = SignalFilter('SampleRate',fs, ...
    'Passband', BW/2, ...
    'StopbandAtten', 60, ...
    'Type', 'lowpass', ...
    'SidebandMode', 'single', ...
    'OutputSameLength', true);

MF_mix = MF.*(exp(1j*2*pi*fc*[0:length(MF(:,1))-1]/(fs))');
for i = 1:sector_num
MF_fir(:,i) = filt.process(MF_mix(:,i));
end
MF_deci =  resample(MF_fir,1,decimate_factor);
% angles_div{1} = -60:0.1:-34-0.1;
% angles_div{2} = -34:0.1:-17-0.1;
% angles_div{3} = -17:0.1:0-0.1;
% angles_div{4} = 0:0.1:17-0.1;
% angles_div{5} = 17:0.1:34-0.1;
% angles_div{6} = 34:0.1:60;


% 生成时间戳
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

% 初始化 attributes
attributes.array_type = 'CDM';
attributes.signal_type = 'Baseband';
attributes.signal_win = 'hamming';
attributes.bandwidth = BW;
attributes.sampling_frequency = fs;
attributes.center_frequency = f0;
attributes.decimate_factor = decimate_factor;
attributes.sector_num = sector_num;
attributes.sector_div = [-60,-34,-17,0,17,34,60.1]; 
attributes.match_filter_data = MF_deci;
attributes.receive_array_num = Nrx;
attributes.receive_array_position = rx_xyz;
attributes.receive_array_win = win;
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
    v_calc_1 = awgn(v_calc,snr_level,'measured');
    v_calc_2 = awgn(v_calc,snr_level,'measured');
    v_calc_3 = awgn(v_calc,snr_level,'measured');
    v_calc_4 = awgn(v_calc,snr_level,'measured');
    v_calc = v_calc_1+v_calc_2+v_calc_3+v_calc_4;
    %% 下变频
    v_mix = v_calc.*repmat(exp(1j*2*pi*f0*[0:length(v_calc(:,1))-1]/fs)',[1 Nrx]);
    %% FIR滤波器
    y_fir = filt.process(v_mix);
    
    %% Decimation
    raw_data = resample(y_fir,1,decimate_factor);
    raw_data_cell{i} = raw_data;
end

%% 写入 HDF5 文件
filename = sprintf('CDMData_%s.h5', timestamp);

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
    h5create(filename, [attrs_group '/receive_array_position'], size(attributes.receive_array_position), 'Datatype', 'single');
    h5write(filename, [attrs_group '/receive_array_position'], attributes.receive_array_position);
end

if isfield(attributes, 'sector_div')
    h5create(filename, [attrs_group '/sector_div'], size(attributes.sector_div), 'Datatype', 'single');
    h5write(filename, [attrs_group '/sector_div'], attributes.sector_div);
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



