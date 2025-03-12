clear all;close all;clc
%% fig. 7a  
file = '191209_13_44_12';
quickAnalysis;

Xall = xyz(:,1); 
Yall = xyz(:,2); 
Zall = xyz(:,3);

block_size = 100;
dff1 = dff0;
Nmf1 = length(xyz);
dff_r = reshape(permute(dff1, [1,3,2]), Nmf1, []);
dff_rz = (dff_r - nanmean(dff_r')') ./ nanstd(dff_r')';

[cc_wL, zc_wL] = bootstrap_cc(dff_rz', L_state', block_size, 40);

sig_lev = 2;

pids = find(zc_wL > sig_lev);
nids = find(zc_wL < -sig_lev);
nsids = find(abs(zc_wL) < sig_lev);

figure('Position',[360, 278, 800, 800]*.5)
hold on; rotate3d on;
grid on;

xyz_mx = [250,250, max(Zall)];
xyz_mn = [0,0, min(Zall)];

pids = logical(zc_wL>sig_lev);
nids = logical(zc_wL<-sig_lev);
nsids = logical(abs(zc_wL)<sig_lev);

h1 = scatter3(Xall(pids), Yall(pids), Zall(pids), 20, 'r', 'LineWidth', 2);
h2 = scatter3(Xall(nids), Yall(nids), Zall(nids), 20, 'b', 'LineWidth', 2);
h3 = scatter3(Xall(nsids), Yall(nsids), Zall(nsids), 20, 'k', 'LineWidth', 2);

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

currentDate = datestr(now, 'yyyy-mm-dd');
savepath2 = fullfile(savepath, 'Figures', 'Figure7', currentDate);
if ~exist(savepath2, 'dir')
    mkdir(savepath2);
end
fileName = ['Spatial_example_locomotion'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% fig.7bc MFB/PM/NM distance

NN_dist = []; NN_dist_pm = []; NN_dist_nsm = []; NN_dist_nm = [];
NN_dist_pm_shuffle = []; NN_dist_nsm_shuffle = [];
NN_dist_nm_shuffle = [];

x_pm_all = []; x_nm_all = []; x_nsm_all = [];
y_pm_all = []; y_nm_all = []; y_nsm_all = [];
z_pm_all = []; z_nm_all = []; z_nsm_all = [];

sig_level = 2;

folder_names = {
      'Animal1';
      'Animal2'
      'Animal3';
      'Animal4';
    };

concat_folder = 'X:\MFB\Concatenation data';
confiles = dir(fullfile(concat_folder, 'concat_*.mat'));

for file_i = 1:4

    dt = 10;
    file = folder_names{file_i}

    if strcmp(file, 'Animal1') == 1
        file = '171212_16_19_37';
        quickAnalysis;
        dff_r = reshape(permute(ROI_dff_all, [3, 2, 1]), [], size(ROI_dff_all, 1))';
    else
        filePath = fullfile(concat_folder, confiles(file_i).name);
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

    % Shuffle
    rand_idx = randperm(N);
    ix_NSM_shuffled = rand_idx(1:N_NSM);
    ix_PM_shuffled = rand_idx(N_NSM+1:N_NSM+N_PM);
    ix_NM_shuffled = rand_idx(N_NSM+N_PM+1:N);

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

% fig 7b
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

% print
fileName = ['XYZ_PM_and_NM.png'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% fig 7c.Ripley function analysis

all_xyz = {}; 
all_zc = {}; 
all_labels = {};  

for file_i = 1:length(confiles)
    dt = 10;
    file = folder_names{file_i};
    
    if strcmp(file, 'Animal1') == 1
        file = '171212_16_19_37';
        quickAnalysis;
        dff_r = reshape(permute(ROI_dff_all, [3, 2, 1]), [], size(ROI_dff_all, 1))';
    else
        filePath = fullfile(concat_folder, confiles(file_i).name);
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

% stats
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

save([savepath '\Processed\Figure 7\RipleyKL.mat'],'overall_results')

%% L Function
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


fileName = 'RipleyL_3D';
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% fig 7 d,e,f
plot_correlationsSpatial = 1;
plot_correlationsSpatial_PMNM = 1;

sig_lev = 0;
sig_level = 2;
block_size = 100;
dt = 10;

folder_names = {
    'Animal1';
    'Animal2';
    'Animal3';
    'Animal4';
};

cc_all_qw = [];  
dd_all_qw = [];
cc_all_as = [];
dd_all_as = [];
cc_all_pm = [];
dd_all_pm = [];
cc_all_nm = [];
dd_all_nm = [];
cc_all_ns = [];
dd_all_ns = [];

for file_i = 1:length(folder_names)
    
    file = folder_names{file_i};
    
    % get MFB activity from concatenated session
    if strcmp(file, 'Animal1') == 1
        file = '171212_16_19_37';
        quickAnalysis
        dff_r = reshape(permute(ROI_dff_all, [3,2,1]), [], size(ROI_dff_all,1))';
    else
        filePath = fullfile(concat_folder, confiles(file_i).name);
        load(filePath)
        dff_r = reshape(permute(ROI_dff_All, [3,2,1]), [], size(ROI_dff_All,1))';
    end
    
    cc1 = corr(dff_r', 'rows', 'complete');
    cc1(eye(length(cc1)) == 1) = nan;
    
    Nmf1 = length(cc1);
    
    dff_rz = (dff_r - nanmean(dff_r')') ./ nanstd(dff_r')';
    
    if ~all(L_state == 0)
        [cc_wL, zc_wL] = bootstrap_cc(dff_rz', L_state', block_size, 200);
    else
        continue
    end
    
    pm = zc_wL > sig_level;
    [~, pm_ids] = find(pm == 1);
    
    nm = zc_wL < -sig_level;
    [~, nm_ids] = find(nm == 1);
    
    nsm = abs(zc_wL) < sig_level;
    [~, ns_ids] = find(nsm == 1);
    
    Nmf = Nmf1;
    
    if Nmf > 1
        
        dd = nan * ones(Nmf, Nmf);
        for i = 1:Nmf
            dd(i,:) = sqrt(nansum((xyz - xyz(i,:)).^2, 2));
        end
        
        [cc_as, zc_as] = bootstrap_cc_pw(dff_rz(:, L_state == 1), block_size, 200); % AS
        
        if length(cc_as) > 1
            cc_as(eye(size(cc_as)) == 1) = nan;
            cc_as(abs(zc_as) < sig_lev) = nan;
            cc_all_as = [cc_all_as; cc_as(~isnan(cc_as))];
            dd_all_as = [dd_all_as; dd(~isnan(cc_as))];
        end
        
        [cc_qw, zc_qw] = bootstrap_cc_pw(dff_rz(:, L_state == 0), block_size, 200); % QW
        
        if length(cc_qw) > 1
            cc_qw(eye(size(cc_qw)) == 1) = nan;
            cc_qw(abs(zc_qw) < sig_lev) = nan;
            cc_all_qw = [cc_all_qw; cc_qw(~isnan(cc_qw))];
            dd_all_qw = [dd_all_qw; dd(~isnan(cc_qw))];
        end
        
        [cc_pm, zc_pm] = bootstrap_cc_pw(dff_rz(pm_ids, L_state == 0), block_size, 200); % QW set L_state==0, AS set L_state==1 
        if length(cc_pm) > 1
            cc_pm(eye(size(cc_pm)) == 1) = nan;
            cc_pm(abs(zc_pm) < sig_lev) = nan;
            cc_all_pm = [cc_all_pm; cc_pm(~isnan(cc_pm))];
            dd_all_pm = [dd_all_pm; dd(~isnan(cc_pm))];
        end
        
        [cc_nm, zc_nm] = bootstrap_cc_pw(dff_rz(nm_ids, L_state == 0), block_size, 200); % QW set L_state==0, AS set L_state==1
        if length(cc_nm) > 1
            cc_nm(eye(size(cc_nm)) == 1) = nan;
            cc_nm(abs(zc_nm) < sig_lev) = nan;
            cc_all_nm = [cc_all_nm; cc_nm(~isnan(cc_nm))];
            dd_all_nm = [dd_all_nm; dd(~isnan(cc_nm))];
        end
        
        [cc_ns, zc_ns] = bootstrap_cc_pw(dff_rz(ns_ids, L_state == 0), block_size, 200); % QW set L_state==0, AS set L_state==1
        if length(cc_ns) > 1
            cc_ns(eye(size(cc_ns)) == 1) = nan;
            cc_ns(abs(zc_ns) < sig_lev) = nan;
            cc_all_ns = [cc_all_ns; cc_ns(~isnan(cc_ns))];
            dd_all_ns = [dd_all_ns; dd(~isnan(cc_ns))];
        end
        
    end
end

if plot_correlationsSpatial
    
    for kk = 1:2
        
        if kk == 1
            x = dd_all_qw;
            y = cc_all_qw;
            which_state = 'Quiet Wakefulness';
            title_color = 'c';
        elseif kk == 2
            x = dd_all_as;
            y = cc_all_as;
            which_state = 'Active State';
            title_color = 'm';
        end

        figure('Position', [160 360 600 400])
        hold on

        title(which_state, 'Color', title_color)
        
        xx = 0:10:200; 

        clids = find((x > 10) .* (y < 1));

        x = x(clids);
        y = y(clids);

        bins = -1:.1:1;
        yhh = nan([length(bins)-1, length(xx)-1]);
        
        ym = zeros(1, length(xx)-1);
        ye = zeros(1, length(xx)-1);
        ym2 = zeros(1, length(xx)-1);
        ye2 = zeros(1, length(xx)-1);
        for j = 1:length(xx)-1
            bin_ids = find((x >= xx(j)) .* (x < xx(j+1)));
            y_all{file_i}{j} = y(bin_ids);

            ym(j) = nanmean(y(bin_ids));
            ye(j) = nanstd(y(bin_ids)); % ./ sqrt(length(y(bin_ids)));

            ym2(j) = nanmean(abs(y(bin_ids)));
            ye2(j) = nanstd(abs(y(bin_ids))); % ./ sqrt(length(y(bin_ids)));
            
            [yh, yx] = histcounts(y(bin_ids), bins);
            yhh(:, j) = yh;
        end
        xc = xx(1:end-1) + diff(xx(1:2))/2;
        yxc = yx(1:end-1) + diff(yx(1:2))/2;
        
        if kk == 1
            ym_sta = ym; ye_sta = ye;
            ym2_sta = ym2; ye2_sta = ye2;
        elseif kk == 2
            ym_run = ym; ye_run = ye;
            ym2_run = ym2; ye2_run = ye2;
        end

        [a, b] = meshgrid(yxc, xc);
        surfc(a, b, yhh'); % / sum(yhh(:))
        colormap('redblue')
         
        rotate3d on
        view(-15, 60)
          
        ylabel({'MFB pairwise', 'distance (μm)'}, 'Rotation', -60)
        xlabel({'MFB pairwise correlation'}, 'rotation', 0)
        zlabel('# MFB pairs')

        set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', ...
            'TickDir', 'out', 'TickLength', [.01 .01])

        currentDate = datestr(now, 'yyyy-mm-dd');
        savepath2 = fullfile(savepath, 'Figures', 'Figure7', currentDate);
        if ~exist(savepath2, 'dir')
            mkdir(savepath2);
        end
        fileName = ['CC_spatial_' which_state];
        fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
        exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

    end
    
    figure('Position', [160 260 350 600])
    
    subplot(211); 
    title('Quiet Wakefulness', 'color', 'c')
    hold on
    
    h1 = plot(xc, ym_sta, 'k-', 'LineWidth', 2);
    errorbar(xc, ym_sta, ye_sta, 'k-')
   
    h2 = plot(xc, ym2_sta, 'r-', 'LineWidth', 2);
    errorbar(xc, ym2_sta, ye2_sta, 'r-')

    legend([h1, h2], {'Average CC', 'Average |CC|'}, 'Location', 'best')
    legend boxoff
    
    xlim([0, 200]);
    ylim([-1, 1])
   
    pos = get(gca, 'Position');
    pos(1) = 0.2; pos(3) = 0.75;
    set(gca, 'Position', pos)
    
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', ...
        'TickDir', 'out', 'TickLength', [.01 .01])

    subplot(212); 
    title('Active State', 'color', 'm')
    hold on
    
    plot(xc, ym_run, 'k-', 'LineWidth', 2);
    errorbar(xc, ym_run, ye_run, 'k-')
    
    plot(xc, ym2_run, 'r-', 'LineWidth', 2);
    errorbar(xc, ym2_run, ye2_run, 'r-')
    
    xlim([0, 200]);
    ylim([-1, 1])
    
    xlabel({'MFB pairwise distance (μm)'})
    ylabel({'MFB pairwise correlation'})
    set(get(gca, 'YLabel'), 'Position', [-30 1.5 0])

    pos = get(gca, 'Position');
    pos(1) = 0.2; pos(3) = 0.75;
    set(gca, 'Position', pos)
    
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', ...
        'TickDir', 'out', 'TickLength', [.01 .01])
    
    fileName = ['CC_spatial_avg_AQ'];
    fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
    exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

end

if plot_correlationsSpatial_PMNM
    
    for kk = 1:3
        
        if kk == 1
            x = dd_all_pm;
            y = cc_all_pm;
            which_state = 'PM';
        elseif kk == 2
            x = dd_all_nm;
            y = cc_all_nm;
            which_state = 'NM';
        elseif kk == 3
            x = dd_all_ns;
            y = cc_all_ns;
            which_state = 'NSM';
        end

        figure('Position', [160 360 600 400])
        hold on
        title(which_state)
        
        xx = 0:10:200; 

        clids = find((x > 10) .* (y < 1));

        x = x(clids);
        y = y(clids);

        bins = -1:.1:1;
        yhh = nan([length(bins)-1, length(xx)-1]);
        
        ym = zeros(1, length(xx)-1);
        ye = zeros(1, length(xx)-1);
        ym2 = zeros(1, length(xx)-1);
        ye2 = zeros(1, length(xx)-1);
        for j = 1:length(xx)-1
            bin_ids = find((x >= xx(j)) .* (x < xx(j+1)));
            y_all{file_i}{j} = y(bin_ids);

            ym(j) = nanmean(y(bin_ids));
            ye(j) = nanstd(y(bin_ids)); % ./ sqrt(length(y(bin_ids)));

            ym2(j) = nanmean(abs(y(bin_ids)));
            ye2(j) = nanstd(abs(y(bin_ids))); % ./ sqrt(length(y(bin_ids)));
            
            [yh, yx] = histcounts(y(bin_ids), bins);
            yhh(:, j) = yh;
        end
        xc = xx(1:end-1) + diff(xx(1:2))/2;
        yxc = yx(1:end-1) + diff(yx(1:2))/2;
        
        if kk == 1
            ym_pm = ym; ye_pm = ye;
            ym2_pm = ym2; ye2_pm = ye2;
        elseif kk == 2
            ym_nm = ym; ye_nm = ye;
            ym2_nm = ym2; ye2_nm = ye2;
        elseif kk == 3
            ym_ns = ym; ye_ns = ye;
            ym2_ns = ym2; ye2_ns = ye2;
        end

        [a, b] = meshgrid(yxc, xc);
        surfc(a, b, yhh');
        colormap('redblue')
         
        rotate3d on
        view(-15, 60)

        ylabel({'MFB pairwise', 'distance (μm)'}, 'Rotation', -60)
        xlabel({'MFB pairwise correlation'}, 'rotation', 0)
        zlabel('# MFB pairs')

        set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', ...
            'TickDir', 'out', 'TickLength', [.01 .01])
        
        fileName = ['CC_spatial_' which_state];
        fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
        exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
                
    end

    figure('Position', [160 60 350 750])
    
    subplot(311); 
    title('PM')
    hold on
    
    h1 = plot(xc, ym_pm, 'k-', 'LineWidth', 2);
    errorbar(xc, ym_pm, ye_pm, 'k-')
   
    h2 = plot(xc, ym2_pm, 'r-', 'LineWidth', 2);
    errorbar(xc, ym2_pm, ye2_pm, 'r-')

    legend([h1, h2], {'Average CC', 'Average |CC|'}, 'Location', 'best')
    legend boxoff
    
    xlim([0, 200]);
    ylim([-1, 1])
    
    pos = get(gca, 'Position');
    pos(1) = 0.2; pos(3) = 0.75;
    set(gca, 'Position', pos)
    
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', ...
        'TickDir', 'out', 'TickLength', [.01 .01])

    subplot(312); 
    title('NM')
    hold on
    
    h3 = plot(xc, ym_nm, 'k-', 'LineWidth', 2);
    errorbar(xc, ym_nm, ye_nm, 'k-')
   
    h4 = plot(xc, ym2_nm, 'b-', 'LineWidth', 2);
    errorbar(xc, ym2_nm, ye2_nm, 'b-')
    
    legend([h3, h4], {'Average CC', 'Average |CC|'}, 'Location', 'best')
    legend boxoff

    xlim([0, 200]);
    ylim([-1, 1])
    ylabel({'MFB pairwise correlation'})

    pos = get(gca, 'Position');
    pos(1) = 0.2; pos(3) = 0.75;
    set(gca, 'Position', pos)
    
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', ...
        'TickDir', 'out', 'TickLength', [.01 .01])
    
    subplot(313); 
    title('NSM')
    hold on
    
    h5 = plot(xc, ym_ns, 'k-', 'LineWidth', 2);
    errorbar(xc, ym_ns, ye_ns, 'k-')
    
    h6 = plot(xc, ym2_ns, '-', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2);
    errorbar(xc, ym2_ns, ye2_ns, '-', 'color', [0.5, 0.5, 0.5])

    legend([h5, h6], {'Average CC', 'Average |CC|'}, 'Location', 'best')
    legend boxoff
    
    xlim([0, 200]);
    ylim([-1, 1])
    
    xlabel({'MFB pairwise distance (μm)'})
    set(get(gca, 'YLabel'), 'Position', [-30 0 0])

    pos = get(gca, 'Position');
    pos(1) = 0.2; pos(3) = 0.75;
    set(gca, 'Position', pos)
    
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', ...
        'TickDir', 'out', 'TickLength', [.01 .01])

    fileName = ['CC_spatial_avg'];
    fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
    exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

end