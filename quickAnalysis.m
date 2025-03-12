path_home = 'X:\2024\work\MossyFibres'; % path to the main code folder
addpath(path_home);
addpath([path_home '/Init']);
Initialize

folder_name = file;
specialNames = {
    '200130_13_21_13 FunctAcq';
    '200130_13_36_14 FunctAcq';
    '200130_13_49_09 FunctAcq';
    '200130_14_02_12 FunctAcq';
    '200130_14_15_24 FunctAcq';
    '200130_14_29_30 FunctAcq';
};

if ismember(folder_name, specialNames)
    special_process = 0;
else
    special_process = 1;
end

% Header info
data_path = [data_home '/' folder_name];
path_header = data_path;

try
    read2__HeaderInfo_patchExp;
catch
    disp('Trying another method to read HeaderInfo');
    read__HeaderInfo_patchExp;
end

% A session has 13 trials instead of 20 in the header
if strcmp(folder_name, '191209_13_44_12')
    Ntr = 13;
end

speed_path = [data_path '/Speed_Data'];
path_speed = [data_path '/Speed_Data'];
read__SpeedData;
Tshift = 0;
laser_times;

file_name = [data_path '/MI_wheel_eyecam', folder_name, '.mat'];
if special_process == 1
    file_name = [data_path, '/MI_wheel', folder_name, ' FunctAcq.mat'];
end

try
    MI_wheel = read_MIdata(file_name, Ntr, T, Tb_tr, Tf_tr, Tshift);
catch
    disp('Try read_MI_wheel in another way')
    file_name = [data_path, '/MI_wheel', folder_name, ' FunctAcq.mat'];
    MI_wheel = read_MIdata(file_name, Ntr, T, Tb_tr, Tf_tr, Tshift);
end

file_name = [data_path '/MI_whiskers', folder_name, '.mat'];
if special_process == 1 
    file_name = [data_path, '/MI_', folder_name, ' FunctAcq.mat'];
end

try
    MI_whisker = read_MIdata(file_name, Ntr, T, Tb_tr, Tf_tr, Tshift);
catch
    disp('Try read_MI_whisker in another way')
    file_name = [data_path, '/MI_whisker', folder_name, ' FunctAcq.mat'];
    MI_whisker = read_MIdata(file_name, Ntr, T, Tb_tr, Tf_tr, Tshift);
end

% - functional
path_functional = data_path;
prefix = 'patch';
bogus_frames = [];

% -- patch numbers to analyze
patch_numbers = [2,4];
if strcmp(folder_name, '171212_16_19_37')
    patch_numbers = [1,2, 4:12, 14:15];
elseif strcmp(folder_name, '180131_14_58_58')
    patch_numbers = [1:17,19:20];
elseif strcmp(folder_name, '180131_15_21_48')
    patch_numbers = [1:17,19:20];
elseif strcmp(folder_name, '180131_15_37_10')
    patch_numbers = [1:17,19:20];
elseif strcmp(folder_name, '191018_13_04_19')
    patch_numbers = [1:6];
elseif strcmp(folder_name, '191018_13_39_41')
    patch_numbers = [1:13];
elseif strcmp(folder_name, '191018_13_56_55')
    patch_numbers = [1:13];
elseif strcmp(folder_name, '191018_14_11_33')
    patch_numbers = [1:13];
elseif strcmp(folder_name, '191018_14_30_00')
    patch_numbers = [1:13];
elseif strcmp(folder_name, '191209_13_44_12')
    patch_numbers = [1:8];
elseif strcmp(folder_name, '191209_14_46_58')
    patch_numbers = [1:8];
elseif strcmp(folder_name, '191209_14_18_13')
    patch_numbers = [1:8];
elseif strcmp(folder_name, '191209_14_32_39')
    patch_numbers = [1:8];
elseif strcmp(folder_name, '191209_15_01_22')
    patch_numbers = [1:8];
elseif strcmp(folder_name, '191209_14_04_14')
    patch_numbers = [1:8];    
end

if special_process == 1 && strcmp(folder_name, '171212_16_19_37')
    bogus_frames = [2807:2814, 3399:3401];
    prefix = 'clean_patch';
    path_functional = fullfile(path_functional, 'Functional_data');
end

Nrois = [];
ROI_dff_all = [];
ROI_fmb_all = [];
f_baselines = [];
backgrounds = [];
Xall = [];
Yall = [];
Zall = [];
for patch_no = patch_numbers
    patch_no
    try
        read__FunctionalData_patchExp;
    catch
        disp('Trying another method to read FunctionalData');
        path_functional = fullfile(data_path, 'Functional_data');
        read__FunctionalData_patchExp;
    end
    Nrois = [Nrois N_miniroi];
    ROI_dff_all = cat(1, ROI_dff_all, ROI_dff);
    ROI_fmb_all = cat(1, ROI_fmb_all, ROI_fmb);
    f_baselines = [f_baselines; f_baseline];
    backgrounds = [backgrounds; background];
    
    Xall = cat(1, Xall, Xc + X0(patch_no));
    Yall = cat(1, Yall, Yc + Y0(patch_no));
    Zall = cat(1, Zall, ones(length(Xc), 1) * Z0(patch_no));
end

xyz = [Xall(:), Yall(:), Zall(:)];
dff_r = reshape(permute(ROI_dff_all, [3,2,1]), [], size(ROI_dff_all,1))';
dff_b = reshape(permute(ROI_fmb_all, [3,2,1]), [], size(ROI_fmb_all,1))';

% adjustment NaNs
if strcmp(folder_name, '200130_13_36_14 FunctAcq')
    dff_r = dff_r(:, 7415:end);
end

speed(isnan(speed)) = 0;
spd_r = fastsmooth(reshape(speed', 1, []), 10, 3, 1);

if strcmp(folder_name, '200130_13_36_14 FunctAcq')
    spd_r = spd_r(7415:end);
end

MI_whisker(isnan(MI_whisker)) = 0;
MI_whisker_r = fastsmooth(reshape(MI_whisker', 1, []), 10, 3, 1);

MI_wheel(isnan(MI_wheel)) = 0;
MI_wheel_r = fastsmooth(reshape(MI_wheel', 1, []), 10, 3, 1);

if strcmp(folder_name, '200130_13_36_14 FunctAcq')
    MI_whisker_r = MI_whisker_r(7415:end);
    MI_wheel_r = MI_wheel_r(7415:end); % remove no active recording
end

nnids = [];
tr = (0:size(dff_r, 2) - 1) * dt / 1e3;
if strcmp(folder_name(1:15), '170712_21_31_26')
    nnids = (tr > 300);
elseif strcmp(folder_name(1:15), '170711_17_47_08')
    nnids = (tr > 280);
elseif strcmp(folder_name(1:15), '171212_16_19_37')
    nnids = (tr > 390);
end

T_state = 3; % seconds
whl = MI_wheel;

whl_s = fastsmooth(reshape(whl', 1, []), 100, 3, 1);  
whl_r = fastsmooth(reshape(whl', 1, []), 10, 3, 1);
whl_rs = fastsmooth(reshape(whl', 1, []), T_state * 1e3 / dt, 3, 1);

if strcmp(folder_name, '200130_13_36_14 FunctAcq')
    whl_s = whl_s(7415:end);
    whl_r = whl_r(7415:end);
    whl_rs = whl_rs(7415:end);
end

whl_b = prctile(whl_r, 10);
whl_th = whl_b + 0.1;
run_r = zeros(size(whl_s));
run_r(whl_s <= whl_th) = 0;
run_r(whl_s >= whl_th) = 1;
run_r(nnids) = nan;

L_th = prctile(whl_rs(~isnan(run_r)), 90);
L_state = zeros(size(whl_rs));
L_state(whl_rs > L_th) = 1;
L_state(nnids) = nan;

T_state = 3; % seconds
wsl = MI_whisker;
wsl_s = fastsmooth(reshape(wsl', 1, []), 100, 3, 1);  
wsl_r = fastsmooth(reshape(wsl', 1, []), 15, 3, 1);
wsl_rs = fastsmooth(reshape(wsl', 1, []), T_state * 1e3 / dt, 3, 1);

if strcmp(folder_name, '200130_13_36_14 FunctAcq')
    wsl_s = wsl_s(7415:end);
    wsl_r = wsl_r(7415:end);
    wsl_rs = wsl_rs(7415:end);
end

wsl_b = prctile(wsl_r, 10);
wsl_th = wsl_b + 0.1;
run_r = zeros(size(wsl_s));
run_r(wsl_s <= wsl_th) = 0;
run_r(wsl_s >= wsl_th) = 1;
run_r(nnids) = nan;

W_th = prctile(wsl_rs(~isnan(run_r)), 90);
W_state = zeros(size(wsl_rs));
W_state(wsl_rs > W_th) = 1;
W_state(nnids) = nan;

if string(folder_name) == "191209_14_46_58"
    L_th = prctile(whl_rs(~isnan(run_r)), 70);
    W_th = prctile(wsl_rs(~isnan(run_r)), 70);
elseif string(folder_name) == "191209_13_44_12"
    L_th = prctile(whl_rs(~isnan(run_r)), 85);
    W_th = prctile(wsl_rs(~isnan(run_r)), 55);
elseif string(folder_name) == "191209_14_04_14"
    L_th = prctile(whl_rs(~isnan(run_r)), 99);
end

windowSize = 300;

QW_Met = (whl_rs < L_th) & (wsl_rs < W_th);
AS_Met = (whl_rs >= L_th) & (wsl_rs >= W_th);

Q_state = zeros(size(QW_Met));
A_state = zeros(size(AS_Met));

Q_segments = bwconncomp(QW_Met);
for i = 1:Q_segments.NumObjects
    segment_indices = Q_segments.PixelIdxList{i};
    if length(segment_indices) >= windowSize
        Q_state(segment_indices) = 1;
    end
end

A_segments = bwconncomp(AS_Met);
for i = 1:A_segments.NumObjects
    segment_indices = A_segments.PixelIdxList{i};
    if length(segment_indices) >= windowSize
        A_state(segment_indices) = 1;
    end
end

A_state(1:150) = false;
L_state = A_state;
my_t = (1:length(spd_r)) * dt;

%% load axon index
mouse_id = '';
axon_idx_mice = '';

if exist('file', "var")
    for i = 1:size(files, 1)
        if strcmp(files{i, 1}, file)
            mouse_id = files{i, 2};
            break;
        end
    end

    if ~isempty(mouse_id)
        axon_idx_mice = ['axon_idx_' mouse_id];
        filepath = fullfile([path_home '\preprocessing\Axon grouping\Axon_index\', axon_idx_mice, '.mat']);
        if exist(filepath, 'file')
            load(filepath);
            fprintf('loading axon index file: %s\n', filepath);
        else
            fprintf('axon index file does not exist: %s\n', filepath);
        end
    else
        fprintf('cannot find the mouse \n');
    end
else
    fprintf('please input file \n');
end

%% activity merged
dff0 = ROI_dff_all; % #MF x trial x time
Nmf0 = size(dff0, 1);
dff0_r = reshape(permute(dff0, [1,3,2]), Nmf0, []);
dff0_rz = (dff0_r - nanmean(dff0_r')') ./ nanstd(dff0_r')';

if strcmp(folder_name, '200130_13_36_14 FunctAcq')
    dff0_rz = dff0_rz(:, 7415:end);
end

group_ids = merged_groups;
Nmf = length(group_ids);
dff_r = zeros(Nmf, size(dff0_r, 2));
MFA_mlt = 0;
MFA_sig = 0;
incl_ids = [];

for i = 1:Nmf
    ids = group_ids{i};
    incl_ids = [incl_ids ids];
end
incl_ids = sort(unique(incl_ids));

for i = 1:Nmf
    j = group_ids{i};
    if numel(j) == 1
        MFA_sig = MFA_sig + 1;
    else
        MFA_mlt = MFA_mlt + 1;
    end
    dff_r(i, :) = nanmean(dff0_r(j, :), 1);
end

t_r = (0:size(dff_r, 2) - 1) * dt / 1e3;
fprintf('single MFA# = %d\n', MFA_sig);
fprintf('multiple MFA# = %d\n', MFA_mlt);

% adjustment for bad recording
if strcmp(folder_name, '200130_13_36_14 FunctAcq')
    dff_r = dff_r(:, 7415:end);
end

% z-scored
dff_rz = (dff_r - nanmean(dff_r')') ./ nanstd(dff_r')';
nans_t = any(isnan(dff_rz), 1);
valid_t = find(~nans_t);

Ngroup = length(group_ids);
Ntrial = size(ROI_dff_all, 2);
Ntime = size(ROI_dff_all, 3);
dff_g = zeros(Ngroup, Ntrial, Ntime);

for i = 1:Ngroup
    ids = group_ids{i};
    dff_g(i, :, :) = nanmean(ROI_dff_all(ids, :, :), 1);
end

dff_gz = zeros(size(dff_g));
for i = 1:Ngroup
    group_data = squeeze(dff_g(i, :, :));
    group_mean = nanmean(group_data(:));
    group_std = nanstd(group_data(:));
    dff_gz(i, :, :) = (group_data - group_mean) / group_std;
end

if strcmp(folder_name, '200130_13_36_14 FunctAcq')
    dff_gz = dff_gz(:, :, 7415:end);
end

G = cellfun(@(x) numel(x) > 1, group_ids);
subgroups = group_ids(G);
