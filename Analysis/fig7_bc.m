%% fig.7bc MFB/PM/NM distance

clear all;close all;clc;

NN_dist = []; NN_dist_pm = []; NN_dist_nsm = []; NN_dist_nm = [];
NN_dist_pm_shuffle = []; NN_dist_nsm_shuffle = [];
NN_dist_nm_shuffle = [];

z_pm_all = []; z_nm_all = []; z_nsm_all = [];
x_pm_all = []; x_nm_all = []; x_nsm_all = [];
y_pm_all = []; y_nm_all = []; y_nsm_all = [];
zc_wL_all = [];

sig_level = 2;

folder_names = {
      '171212_16_19_37';
      'superBernie';
      'superBill';
      'superJeremy';
    };

folderPath = 'X:\MFB\MFB_AH_2023\Correlation_data';
confile = dir(fullfile(folderPath, 'concat_*.mat'));


for file_i = 1:length(confile)

    dt = 10;
    file = folder_names{file_i};

    if strcmp(file, '171212_16_19_37') == 1
        quickAnalysis;
        dff_r = reshape(permute(ROI_dff_all, [3, 2, 1]), [], size(ROI_dff_all, 1))';
    else
        filePath = fullfile(folderPath, confile(file_i).name);
        load(filePath);
        dff_r = reshape(permute(ROI_dff_All, [3, 2, 1]), [], size(ROI_dff_All, 1))';
    end

    dff_rz = (dff_r - nanmean(dff_r')') ./ nanstd(dff_r')';

    [cc_wL, zc_wL] = bootstrap_cc(dff_rz', L_state', 100, 200);

    zc_wL_all = [zc_wL_all, zc_wL];

    % Info about xyz
    file_lists = {'171212_16_19_37', '191209_13_44_12', '200130_13_21_13 FunctAcq', '191018_14_11_33'};
    file = file_lists{file_i};
    quickAnalysis;

    dd_all = pdist2(xyz(:,:), xyz(:,:));
    dd_all(eye(size(dd_all)) == 1) = nan;

    % Find PM/NM
    pm = zc_wL > sig_level;
    [~, pm_ids] = find(pm == 1);

    nsm = abs(zc_wL) < sig_level;
    [~, nsm_ids] = find(nsm == 1);

    nm = zc_wL < -sig_level;
    [~, nm_ids] = find(nm == 1);

    % NN distance
    N = length(xyz);
    N_PM = length(pm_ids);
    N_NM = length(nm_ids);
    N_NSM = length(nsm_ids);

    for i = 1:length(pm_ids)
        distances = dd_all(pm_ids(i), pm_ids);
        distances(i) = Inf;
        NN_dist_pm = [NN_dist_pm, min(distances)];
    end

    for i = 1:length(nm_ids)
        distances = dd_all(nm_ids(i), nm_ids);
        distances(i) = Inf;
        NN_dist_nm = [NN_dist_nm, min(distances)];
    end

    for i = 1:length(nsm_ids)
        distances = dd_all(nsm_ids(i), nsm_ids);
        distances(i) = Inf;
        NN_dist_nsm = [NN_dist_nsm, min(distances)];
    end

    % Shuffle
    ix_NSM_shuffled = randsample(1:N, N_NSM);
    ix_PM_shuffled = randsample(setdiff(1:N, ix_NSM_shuffled), N_PM);
    ix_NM_shuffled = randsample(setdiff(1:N, ix_PM_shuffled), N_NM);

    for i = 1:length(ix_PM_shuffled)
        distances = dd_all(ix_PM_shuffled(i), ix_PM_shuffled);
        distances(i) = Inf;
        NN_dist_pm_shuffle = [NN_dist_pm_shuffle, min(distances)];
    end

    for i = 1:length(ix_NM_shuffled)
        distances = dd_all(ix_NM_shuffled(i), ix_NM_shuffled);
        distances(i) = Inf;
        NN_dist_nm_shuffle = [NN_dist_nm_shuffle, min(distances)];
    end

    for i = 1:length(ix_NSM_shuffled)
        distances = dd_all(ix_NSM_shuffled(i), ix_NSM_shuffled);
        distances(i) = Inf;
        NN_dist_nsm_shuffle = [NN_dist_nsm_shuffle, min(distances)];
    end

    xyz2 = xyz(:,:);
    z_pm = xyz2(pm_ids, 3);
    z_nm = xyz2(nm_ids, 3);
    z_nsm = xyz2(nsm_ids, 3);
    z_pm_all = [z_pm_all; z_pm];
    z_nm_all = [z_nm_all; z_nm];
    z_nsm_all = [z_nsm_all; z_nsm];

    x_pm = xyz2(pm_ids, 1);
    x_nm = xyz2(nm_ids, 1);
    x_nsm = xyz2(nsm_ids, 1);
    x_pm_all = [x_pm_all; x_pm];
    x_nm_all = [x_nm_all; x_nm];
    x_nsm_all = [x_nsm_all; x_nsm];

    y_pm = xyz2(pm_ids, 2);
    y_nm = xyz2(nm_ids, 2);
    y_nsm = xyz2(nsm_ids, 2);
    y_pm_all = [y_pm_all; y_pm];
    y_nm_all = [y_nm_all; y_nm];
    y_nsm_all = [y_nsm_all; y_nsm];
end


pm_all = zc_wL_all > sig_level;
[~, pm_ids_all] = find(pm_all == 1);

nm_all = zc_wL_all < -sig_level;
[~, nm_ids_all] = find(nm_all == 1);

nsm_all = abs(zc_wL_all)<sig_level;
[~, nsm_ids_all] = find(nsm_all == 1);


% plot 7b
subplot(3,1,1)
histogram(x_pm_all, 'BinWidth', 10, 'FaceColor', [1, 1, 1], 'EdgeColor', 'r', 'LineWidth', 1,'Normalization','count');
hold on;
histogram(x_nm_all, 'BinWidth', 10, 'FaceColor', [1, 1, 1], 'EdgeColor', 'b', 'LineWidth', 1,'Normalization','count');
xlabel('X (μm)');
ylabel('Number');
set(gca, 'LineWidth', 1, 'FontSize', 15, ...
    'FontName'   , 'Helvetica', ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'   , ...
    'ZMinorTick'  , 'off'  )

[hx,px] = kstest2(x_pm_all, x_nm_all);
fprintf('KS test result for X: h = %d\n', hx);
fprintf('P-value: %f\n', px);
ylim([0 60])
yticks([0 30 60])

subplot(3,1,2)
histogram(y_pm_all, 'BinWidth', 10, 'FaceColor', [1, 1, 1], 'EdgeColor', 'r', 'LineWidth', 1,'Normalization','count');
hold on;
histogram(y_nm_all, 'BinWidth', 10, 'FaceColor', [1, 1, 1], 'EdgeColor', 'b', 'LineWidth', 1,'Normalization','count');
xlabel('Y (μm)');
ylabel('Number');
set(gca, 'LineWidth', 1, 'FontSize', 15, ...
    'FontName'   , 'Helvetica', ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'   , ...
    'ZMinorTick'  , 'off'  )

[hy,py] = kstest2(y_pm_all, y_nm_all);
fprintf('KS test result for Y: h = %d\n', hy);
fprintf('P-value: %f\n', py);


set(gcf, 'Position', [100, 100, 450, 550]);
ylim([0 60])
yticks([0 30 60])

subplot(3,1,3)

h1 = histogram(z_pm_all, 'BinWidth', 5, 'FaceColor', [1, 1, 1], 'EdgeColor', ...
    'r', 'LineWidth', 1,'Normalization','count');
hold on;
h2 = histogram(z_nm_all, 'BinWidth', 5, 'FaceColor', [1, 1, 1], 'EdgeColor', ...
    'b', 'LineWidth', 1,'Normalization','count');
xlabel('Z (μm)');
ylabel('Number');
ylim([0 140])
yticks([0 75 140])
set(gca, 'LineWidth', 1, 'FontSize', 15, ...
    'FontName'   , 'Helvetica', ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'   , ...
    'ZMinorTick'  , 'off'  )

legend([h1, h2], 'PM', 'NM','location','best');
legend('boxoff');

[hz,pz] = kstest2(z_pm_all, z_nm_all);
fprintf('KS test result for Z: h = %d\n', hz);
fprintf('P-value: %f\n', pz);

mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
folderPath = fullfile(mfbFolderPath, 'Figures', 'Figure7', currentDate);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end


% print
fileName = ['XYZ_PM_and_NM.png'];
fullFilePath = fullfile(folderPath, fileName);
print(fullFilePath, '-dpng', '-r300');


%% 7c
figure;
h1 = histogram(NN_dist_pm, 'BinWidth', 5, 'FaceColor', [1, 1, 1], 'EdgeColor', 'r', 'LineWidth', 1); % RGB color for red
xlabel('Distance (μm)');
ylabel('Number');
xticks([0 10 20 30 40 50 60 70]);
yticks([0 90 180]);
ylim([0 180]);
xlim([-5 75]);
hold on;
h2 = histogram(NN_dist_pm_shuffle, 'BinWidth', 5, 'FaceColor', [1, 1, 1], 'EdgeColor', [0.5, 0.5, 0.5], 'LineWidth', 1); % 灰色

legend([h1, h2], 'PM', 'Shuffle');
legend('boxoff');

fprintf("NN_dist_PM size=%d\n", size(NN_dist_pm, 2));
fig = gcf;
fig.Units = 'pixels';
fig.Position = [100 100 450 200];
set(gca, 'LineWidth', 1, 'FontSize', 15, ...
    'FontName', 'Helvetica', ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'TickLength', [.02 .02], ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'ZMinorTick', 'off');

[p, h] = signrank(NN_dist_pm_shuffle, NN_dist_pm);
fprintf('sign rank P-value: %f\n', p);

fileName = 'NN_dist_PM.png';
fullFilePath = fullfile(folderPath, fileName);
print(fullFilePath, '-dpng', '-r300');

% NM
figure;
h3 = histogram(NN_dist_nm, 'BinWidth', 5, 'FaceColor', [1, 1, 1], 'EdgeColor', 'b', 'LineWidth', 1.5);
xlabel('Distance (μm)');
ylabel('Number');
xticks([0 10 20 30 40 50 60 70]);
yticks([0 90 180]);
ylim([0 180]);
xlim([-5 75]);
hold on;
h4 = histogram(NN_dist_nm_shuffle, 'BinWidth', 5, 'FaceColor', [1, 1, 1], 'EdgeColor', [0.5, 0.5, 0.5], 'LineWidth', 1.5);

legend([h3, h4], 'NM', 'Shuffle');
legend('boxoff');

fprintf("NN_dist_NM size=%d\n", size(NN_dist_nm, 2));
fig = gcf;
fig.Units = 'pixels';
fig.Position = [100 100 450 200];
set(gca, 'LineWidth', 1, 'FontSize', 15, ...
    'FontName', 'Helvetica', ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'TickLength', [.02 .02], ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'ZMinorTick', 'off');

[p, h] = signrank(NN_dist_nm_shuffle, NN_dist_nm);
fprintf('sign rank P-value: %f\n', p);

fileName = 'NN_dist_NM.png';
fullFilePath = fullfile(folderPath, fileName);
print(fullFilePath, '-dpng', '-r300');

% NSM
figure;
h5 = histogram(NN_dist_nsm, 'BinWidth', 5, 'FaceColor', [1, 1, 1], 'EdgeColor', [0, 0, 0], 'LineWidth', 1.5);
xlabel('Distance (μm)');
ylabel('Number');
xticks([0 10 20 30 40 50 60 70]);
yticks([0 90 180]);
ylim([0 180]);
xlim([-5 75]);
hold on;
h6 = histogram(NN_dist_nsm_shuffle, 'BinWidth', 5, 'FaceColor', [1, 1, 1], 'EdgeColor', [0.5, 0.5, 0.5], 'LineWidth', 1.5);

legend([h5, h6], 'NSM', 'Shuffle');
legend('boxoff');

fprintf("NN_dist_NSM size=%d\n", size(NN_dist_nsm, 2));
fig = gcf;
fig.Units = 'pixels';
fig.Position = [100 100 450 200];
set(gca, 'LineWidth', 1, 'FontSize', 15, ...
    'FontName', 'Helvetica', ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'TickLength', [.02 .02], ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off', ...
    'ZMinorTick', 'off');

[p, h] = signrank(NN_dist_nsm_shuffle, NN_dist_nsm);
fprintf('sign rank P-value: %f\n', p);

fileName = 'NN_dist_NSM.png';
fullFilePath = fullfile(folderPath, fileName);
print(fullFilePath, '-dpng', '-r300');

% Save data
filePath = fullfile('X:\MFB\MFB_AH_2023\spatial_dist_all.mat');
save(filePath, 'NN_dist_pm_shuffle', 'NN_dist_pm', 'NN_dist_nm_shuffle', 'NN_dist_nm', 'NN_dist_nsm_shuffle', 'NN_dist_nsm');



