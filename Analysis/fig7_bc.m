
%% fig. 7a  
file='191209_13_44_12';
quickAnalysis;

Xall = xyz(:,1); Yall = xyz(:,2); Zall = xyz(:,3);

block_size=100;
dff1=dff0;
Nmf1 = length(xyz);
dff_r = reshape(permute(dff1, [1,3,2]), Nmf1, []);
dff_rz = (dff_r - nanmean(dff_r')') ./ nanstd(dff_r')';

[cc_wL, zc_wL] = bootstrap_cc(dff_rz', L_state', block_size, 40);

sig_lev = 2;


pids = find(zc_wL > sig_lev);
nids = find(zc_wL < -sig_lev);
nsids = find(abs(zc_wL)<sig_level);


figure('Position',[360, 278, 800, 800]*.5)
hold on; rotate3d on;
grid on;

xyz_mx = [250,250, max(Zall)];
xyz_mn = [0,0, min(Zall)];

% plot the box
pids = logical(zc_wL>sig_lev);
nids = logical(zc_wL<-sig_lev);
nsids = logical(abs(zc_wL)<sig_lev);

%scatter3(Xall(:), Yall(:), Zall(:), 40, [1,1,1], 'LineWidth', 1, 'MarkerFaceColor',[.9,.9,.9]);

h1=scatter3(Xall(pids), Yall(pids), Zall(pids), 20, 'r', 'LineWidth', 2);
h2=scatter3(Xall(nids), Yall(nids), Zall(nids), 20, 'b', 'LineWidth', 2);
h3=scatter3(Xall(nsids), Yall(nsids), Zall(nsids), 20, 'k', 'LineWidth', 2);

legend([h1,h2,h3], {'PM', 'NM', 'NS'}, 'Location','northwest')
legend boxoff


xlim([0,xyz_mx(1)]);
ylim([0,xyz_mx(2)]);
zlim([xyz_mn(3)-10,xyz_mx(3)+10]);
        
zticks([-100, -80, -60, -40, -20, 0]);

xlabel('X (μm)'); ylabel('Y (μm)'); zlabel('Z (μm)');
view(-28,28);

set(gca, 'LineWidth', 1, 'FontSize', 15, ...
    'FontName'   , 'Helvetica', ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'   , ...
    'ZMinorTick'  , 'off'  )

mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
folderPath = fullfile(mfbFolderPath, 'Figures', 'Figure7', currentDate);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
fileName = ['Spatial_example_locomotion'];
fullFilePathPDF = fullfile(folderPath, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


%% fig.7bc MFB/PM/NM distance

clear all;close all;clc;

NN_dist = []; NN_dist_pm = []; NN_dist_nsm = []; NN_dist_nm = [];
NN_dist_pm_shuffle = []; NN_dist_nsm_shuffle = [];
NN_dist_nm_shuffle = [];

x_pm_all = []; x_nm_all = []; x_nsm_all = [];
y_pm_all = []; y_nm_all = []; y_nsm_all = [];
z_pm_all = []; z_nm_all = []; z_nsm_all = [];

sig_level = 2;

folder_names = {
      '171212_16_19_37';
      'superBernie';
      'superBill';
      'superJeremy';
    };

folderPath = 'X:\MFB\MFB_AH_2023\Correlation_data\4mice';
confile = dir(fullfile(folderPath, 'concat_*.mat'));


for file_i = 1:4%length(confile)

    dt = 10;
    file = folder_names{file_i}

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

%     for i = 1:length(pm_ids)
%         distances = dd_all(pm_ids(i), pm_ids);
%         distances(i) = Inf;
%         NN_dist_pm = [NN_dist_pm, min(distances)];
%     end
% 
%     for i = 1:length(nm_ids)
%         distances = dd_all(nm_ids(i), nm_ids);
%         distances(i) = Inf;
%         NN_dist_nm = [NN_dist_nm, min(distances)];
%     end
% 
%     for i = 1:length(nsm_ids)
%         distances = dd_all(nsm_ids(i), nsm_ids);
%         distances(i) = Inf;
%         NN_dist_nsm = [NN_dist_nsm, min(distances)];
%     end

    % Shuffle
    rand_idx = randperm(N);
    ix_NSM_shuffled = rand_idx(1:N_NSM);
    ix_PM_shuffled = rand_idx(N_NSM+1:N_NSM+N_PM);
    ix_NM_shuffled = rand_idx(N_NSM+N_PM+1:N);
% 
%     for i = 1:length(ix_PM_shuffled)
%         distances = dd_all(ix_PM_shuffled(i), ix_PM_shuffled);
%         distances(i) = Inf;
%         NN_dist_pm_shuffle = [NN_dist_pm_shuffle, min(distances)];
%     end
% 
%     for i = 1:length(ix_NM_shuffled)
%         distances = dd_all(ix_NM_shuffled(i), ix_NM_shuffled);
%         distances(i) = Inf;
%         NN_dist_nm_shuffle = [NN_dist_nm_shuffle, min(distances)];
%     end
% 
%     for i = 1:length(ix_NSM_shuffled)
%         distances = dd_all(ix_NSM_shuffled(i), ix_NSM_shuffled);
%         distances(i) = Inf;
%         NN_dist_nsm_shuffle = [NN_dist_nsm_shuffle, min(distances)];
%     end

%     if length(idx_pm) > 3
%         NN_dist_pm = [NN_dist_pm, min(dd_all(idx_pm, idx_pm) + diag(Inf(size(idx_pm))))];
%         NN_dist_pm_shuffle = [NN_dist_pm_shuffle, min(dd_all(ix_PM_shuffled, ix_PM_shuffled) + diag(Inf(size(ix_PM_shuffled))))];
%     end
% 
%     if length(idx_nm) > 3
%         NN_dist_nm = [NN_dist_nm, min(dd_all(idx_nm, idx_nm) + diag(Inf(size(idx_nm))))];
%         NN_dist_nm_shuffle = [NN_dist_nm_shuffle, min(dd_all(ix_NM_shuffled, ix_NM_shuffled) + diag(Inf(size(ix_NM_shuffled))))];
%     end
% 
%     if length(idx_nsm) > 3
%         NN_dist_nsm = [NN_dist_nsm, min(dd_all(idx_nsm, idx_nsm) + diag(Inf(size(idx_nsm))))];
%         NN_dist_nsm_shuffle = [NN_dist_nsm_shuffle, min(dd_all(ix_NSM_shuffled, ix_NSM_shuffled) + diag(Inf(size(ix_NSM_shuffled))))];
%     end


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
fullFilePathPDF = fullfile(folderPath, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


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

filePath = fullfile('X:\MFB\MFB_AH_2023\spatial_dist_all.mat');
save(filePath, 'NN_dist_pm_shuffle', 'NN_dist_pm', 'NN_dist_nm_shuffle', 'NN_dist_nm', 'NN_dist_nsm_shuffle', 'NN_dist_nsm');

%% fig 7c.Ripley function analysis
folder_names = {
      '171212_16_19_37';
      'superBernie';
      'superBill';
      'superJeremy';
    };

folderPath = 'X:\MFB\MFB_AH_2023\Correlation_data\4mice';
confile = dir(fullfile(folderPath, 'concat_*.mat'));
all_xyz = {}; 
all_zc = {}; 
all_labels = {};  

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

    all_xyz{file_i} = xyz;
    all_zc{file_i} = zc_wL;
    
    sig_level = 2;
    pm_idx = find(zc_wL > sig_level);
    nm_idx = find(zc_wL < -sig_level);
    nsm_idx = find(abs(zc_wL) < sig_level);
    all_labels{file_i}.pm = pm_idx;
    all_labels{file_i}.nm = nm_idx;
    all_labels{file_i}.nsm = nsm_idx;
end

%% stats
for s = 1:length(all_xyz) 
    xyz_session = all_xyz{s};
    labels = all_labels{s};
   
    [r_vals, K_overall, L_overall] = ripleyKL3D(xyz_session,0);
    
    if length(labels.pm) >= 2
       [~, K_pm, L_pm] = ripleyKL3D(xyz_session(labels.pm, :),0);
    else
       K_pm = nan(size(r_vals)); L_pm = nan(size(r_vals));
    end

    if length(labels.nm) >= 2
       [~, K_nm, L_nm] = ripleyKL3D(xyz_session(labels.nm, :),0);
    else
       K_nm = nan(size(r_vals)); L_nm = nan(size(r_vals));
    end

    if length(labels.nsm) >= 2
       [~, K_nsm, L_nsm] = ripleyKL3D(xyz_session(labels.nsm, :),0);
    else
       K_nsm = nan(size(r_vals)); L_nsm = nan(size(r_vals));
    end
    

    overall_results{s}.r = r_vals;
    overall_results{s}.K = K_overall;
    overall_results{s}.L = L_overall;
    overall_results{s}.K_pm = K_pm;
    overall_results{s}.L_pm = L_pm;
    overall_results{s}.K_nm = K_nm;
    overall_results{s}.L_nm = L_nm;
    overall_results{s}.K_nsm = K_nsm;
    overall_results{s}.L_nsm = L_nsm;
end

save("X:\MFB\Processed\Figure 7\RipleyKL.mat",'overall_results')

%% Aggregated Analysis using overall_results (L Function)
numSess = length(overall_results);
r_vals = overall_results{1}.r; 
n_r = length(r_vals);

L_all = zeros(numSess, n_r);
L_all_pm = nan(numSess, n_r);
L_all_nm = nan(numSess, n_r);
L_all_nsm = nan(numSess, n_r);

for s = 1:numSess
    L_all(s,:) = overall_results{s}.L;       % Overall L
    L_all_pm(s,:) = overall_results{s}.L_pm;   % PM L
    L_all_nm(s,:) = overall_results{s}.L_nm;   % NM L
    L_all_nsm(s,:) = overall_results{s}.L_nsm; % NSM L
end

meanL_overall = nanmean(L_all, 1);
seL_overall = nanstd(L_all, 0, 1) / sqrt(numSess);
meanL_pm = nanmean(L_all_pm, 1);
seL_pm = nanstd(L_all_pm, 0, 1) / sqrt(numSess);
meanL_nm = nanmean(L_all_nm, 1);
seL_nm = nanstd(L_all_nm, 0, 1) / sqrt(numSess);
meanL_nsm = nanmean(L_all_nsm, 1);
seL_nsm = nanstd(L_all_nsm, 0, 1) / sqrt(numSess);

meanL_overall(isnan(meanL_overall)) = 0;
seL_overall(isnan(seL_overall)) = 0;

meanL_pm(isnan(meanL_pm)) = 0;
seL_pm(isnan(seL_pm)) = 0;

meanL_nm(isnan(meanL_nm)) = 0;
seL_nm(isnan(seL_nm)) = 0;

meanL_nsm(isnan(meanL_nsm)) = 0;
seL_nsm(isnan(seL_nsm)) = 0;

L_csr = r_vals;

figure('Position',[100 100 250 750]);
% Overall L Function
% subplot(1,3,1);
% fill([r_vals, fliplr(r_vals)], [meanL_overall+seL_overall, fliplr(meanL_overall-seL_overall)], [0.7 1 0.7],...
%     'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility','off');
% hold on;
% h_obs = plot(r_vals, meanL_overall, '-', 'Color', [0 0.5 0], 'LineWidth', 2);
% h_csr = plot(r_vals, L_csr, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
% xlabel('r (μm)'); ylabel('L(r)');
% title('Overall');
% legend([h_obs, h_csr], {'Observed','CSR'}, 'Location','best','box','off');
% set(gca, 'FontSize', 13);

% PM L Function
subplot(3,1,1);
fill([r_vals, fliplr(r_vals)], [meanL_pm+seL_pm, fliplr(meanL_pm-seL_pm)], [1 0.8 0.8],...
    'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility','off');
hold on;

h_obs = plot(r_vals, meanL_pm, '-', 'Color', 'r', 'LineWidth', 2);
h_csr = plot(r_vals, L_csr, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlabel('r (μm)'); ylabel('L(r)');
title('PM');
legend([h_obs, h_csr], {'Observed','CSR'}, 'Location','northwest','box','off');
set(gca, 'FontSize', 13,'XLim',[0 50],'XTick',[0 50],'YLim',[0 50],'YTick',[0 50]);

% NM L Function
subplot(3,1,2);
fill([r_vals, fliplr(r_vals)], [meanL_nm+seL_nm, fliplr(meanL_nm-seL_nm)], [0.8 0.8 1],...
    'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility','off');
hold on;

h_obs = plot(r_vals, meanL_nm, '-', 'Color', 'b', 'LineWidth', 2);
h_csr = plot(r_vals, L_csr, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlabel('r (μm)'); ylabel('L(r)');
title('NM');
legend([h_obs, h_csr], {'Observed','CSR'}, 'Location','northwest','box','off');
set(gca, 'FontSize', 13,'XLim',[0 50],'XTick',[0 50],'YLim',[0 50],'YTick',[0 50]);

% NSM L Function
subplot(3,1,3);
fill([r_vals, fliplr(r_vals)], [meanL_nsm+seL_nsm, fliplr(meanL_nsm-seL_nsm)], [0.7 0.7 0.7],...
    'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility','off');
hold on;

h_obs = plot(r_vals, meanL_nsm, '-', 'Color', 'k', 'LineWidth', 2);
h_csr = plot(r_vals, L_csr, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
xlabel('r (μm)'); ylabel('L(r)');
title('NSM');
legend([h_obs, h_csr], {'Observed','CSR'}, 'Location','northwest','box','off');
set(gca, 'FontSize', 13,'XLim',[0 50],'XTick',[0 50],'YLim',[0 50],'YTick',[0 50]);

mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
saveFolder = fullfile(mfbFolderPath, 'Figures', 'Figure7', currentDate);
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

fileName = 'RipleyL_3D';
fullFilePathPDF = fullfile(saveFolder, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

