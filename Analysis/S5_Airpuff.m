%
%   This script performs analysis of peak latency in response 
%   to a puff stimulus across multiple experimental datasets. It includes 
%   several sections:
%
%     1. S5.a: Displays separate subplots for low, medium, and high latency groups,
%        as well as an average trace, allowing visual comparison across groups.
%     2. S5.b-e: For each experimental folder (specified in the folder list),
%        the script:
%          - Loads and preprocesses data using the "quickAnalysis" routine.
%          - Calculates the latency of peak and onset activity relative to a puff 
%            stimulus (using the "peak_latency" function).
%          - Generates heatmaps of normalized ΔF/F data, sorted by onset or peak 
%            time.
%          - Plots histograms and scatter plots to display the distribution of 
%            peak/onset latencies.
%     3. S5.f: Groups the experimental datasets into four categories (Nigel, Bernie, 
%        Bill, and Jeremy), computes the mean and standard deviation of the time 
%        to peak for each group, and generates a boxplot to compare these metrics.
%
% Note:
%   - Parameters such as "tpuff", "t_before", and "t_after" are set based on 
%     experimental conditions.
%
% Author: Yizhou Xie
% Date: 2024/7/6
% Updated: 2025/3/5

all_folder_names = {
    '171212_16_19_37'
    '191018_13_39_41';
    '191018_13_56_55';
    '200130_13_21_13 FunctAcq';
    '200130_14_29_30 FunctAcq';
    '191209_13_44_12';
    '191209_14_18_13';
    '191209_14_46_58'; 
};

%% S5.b-e
low_latency_t2p_all = [];

for file_i = 1:length(folder_names)
    dt = 10;
    file = folder_names{file_i};
    quickAnalysis;
    
    plot__peakLatency = 1;
    
    if plot__peakLatency
        ColorCodes
        tpuff = 5e3;
        
        after_puffs = [7, 7];
        Ntot = size(dff_gz, 1);
        
        % latency of the peak activity after puff
        for ijk = 1:2
            after_puff = after_puffs(ijk);
            
            Tb1 = 4.5e3;
            Tb2 = tpuff;
            Ta1 = tpuff;
            Ta2 = after_puff * 1e3;
            
            time_wind_before = Tb2 - Tb1;
            time_wind_after = Ta2 - Ta1;
            
            tis_before = find((T >= Tb1) .* (T < Tb2));
            tis_after = find((T >= Ta1) .* (T <= Ta2));
            
            tt = T(tis_after) / 1e3;
            
            % -
            spd_th = 0.1;
            spdm = nanmean(speed(:, tis_before), 2);
            
            ROI_dff_mfa = dff_gz(:, :, :);
            all_dff_m = squeeze(nanmean(ROI_dff_mfa(:, spdm < spd_th, :), 2));
            
            disp(['fraction filtered: ', num2str(sum(spdm(:) < spd_th) / numel(spdm) * 100), '%'])
            
            all_dff_m_tp = all_dff_m';
            zz = zeros(size(all_dff_m_tp));
            mu = nanmean(all_dff_m_tp, 1);
            sigma = nanstd(all_dff_m_tp, 0, 1);
            
            for i = 1:size(all_dff_m_tp, 2)
                zz(:, i) = (all_dff_m_tp(:, i) - mu(i)) / sigma(i);
            end
            zz = zz';
            
            [peak_time, onset_time, rel_inc, maxima, minima] = peak_latency(zz, tis_before, tis_after);
            
            if ijk == 1
                [a, b] = sort(onset_time);
                fig_suffix = 'onsetSorted_';
            elseif ijk == 2
                [a, b] = sort(peak_time);
                fig_suffix = 'peakSorted_';
            end
            
            figure('Position', [360 20 369 678])
            hold on
            imagesc(zz(b, :), [0 5])
            plot([5e3 / dt, 5e3 / dt], [0, Ntot], '--', 'LineWidth', 2, 'color', 'w')
            
            if ijk == 1
                plot((tpuff + a) / dt, 1:Ntot, '-', 'color', 'g', 'LineWidth', 2)
            elseif ijk == 2
                plot((tpuff + a) / dt, 1:Ntot, '-', 'color', 'y', 'LineWidth', 2)
            end
            
            xticks((3:10) * 1e3 / dt)
            xticklabels({'-2'; '-1'; '0'; '1'; '2'; '3'; ''; ''})
            
            xlabel('Time relative to puff (s)')
            ylabel('Sorted MFA #')
            
            xlim([3e3, 7e3] / dt)
            ylim([.5, Ntot + .5])
            
            set(gca, 'LineWidth', 1, 'FontSize', 20, ...
                'FontName', 'Helvetica', ...
                'Box', 'off', ...
                'TickDir', 'out', ...
                'TickLength', [.02 .02], ...
                'XMinorTick', 'off', ...
                'YMinorTick', 'off', ...
                'ydir', 'reverse')
            
            mfbFolderPath = 'X:\MFB';
            currentDate = datestr(now, 'yyyy-mm-dd');
            folderPath = fullfile(mfbFolderPath, 'Figures\Supp.Figures\S5\', currentDate);
            if ~exist(folderPath, 'dir')
                mkdir(folderPath);
            end
            
            fileName = [fig_suffix 'HeatMap_dff_' num2str(after_puff) '_' char(file)];
            fullFilePathPDF = fullfile(folderPath, [fileName, '.pdf']);
            exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
            
            % peak / onset time ranges
            if ijk == 1
                my_times = onset_time;
                fig_suffix = 'OnsetTimes_';
                xlbl = 'Onset time (ms)';
                bins = 0:10:(Ta2 - tpuff);
            elseif ijk == 2
                my_times = peak_time;
                fig_suffix = 'PeakTimes_';
                xlbl = 'Time to peak (ms)';
                bins = 0:20:(Ta2 - tpuff);
            end
            
            figure('Position', [360, 278, 600, 250], 'Visible', 'on'); 
            hold on
            axis off;
            axes('YAxisLocation', 'left', ...
                'Color', 'none', ...
                'YColor', 'none', 'XColor', 'k'); 
            hold on
            
            cl = 'k';
            plot(my_times, rand(size(my_times)), '^', ...
                'Markersize', 10, 'color', cl, 'MarkerFaceColor', 'none')
            xlabel(xlbl);
            
            set(gca, 'LineWidth', 1, 'FontSize', 20, ...
                'FontName', 'Helvetica', ...
                'Box', 'off', ...
                'TickDir', 'out', ...
                'TickLength', [.02 .02], ...
                'XMinorTick', 'off', ...
                'YMinorTick', 'off')
            
            [hst1, xx] = histcounts(my_times, bins, 'Normalization', 'count');
            hst1 = fastsmooth(hst1, 3, 3, 1);
            xx = xx(1:end-1) + diff(xx(1:2)) / 2;
            
            cl1 = [.5, .5, .5];
            
            figure('Position', [360, 278, 600, 250], 'Visible', 'on'); 
            hold on
            histogram(my_times, bins, 'Normalization', 'count', ...
                'LineWidth', .1, 'DisplayStyle', 'bar', 'FaceAlpha', .1875, 'FaceColor', cl1, 'EdgeColor', cl1)
            
            plot(xx, hst1, 'LineWidth', 2, 'Color', cl1)
            ylabel('#');
            
            xlabel(xlbl);
            
            set(gca, 'LineWidth', 1, 'FontSize', 20, ...
                'FontName', 'Helvetica', ...
                'Box', 'off', ...
                'TickDir', 'out', ...
                'TickLength', [.02 .02], ...
                'XMinorTick', 'off', ...
                'YMinorTick', 'off')
        end
    end
    
    % S5.b
    tpuff = 5e3;
    
    t_before = 0.5 * 1e3;
    t_after = 2 * 1e3;
    
    Tb1 = tpuff - t_before;
    Tb2 = tpuff;
    Ta1 = tpuff;
    Ta2 = tpuff + t_after;
    
    tis_before = find((T >= Tb1) .* (T < Tb2));
    tis_after = find((T >= Ta1) .* (T <= Ta2));
    
    spd_th = .1;
    
    spdm = nanmean(speed(:, tis_before), 2);
    
    x = squeeze(nanmean(dff_gz(:, spdm < spd_th, :), 2));
    
    x_z = (x - nanmean(x')') ./ nanstd(x')';
    
    [peak_time, onset_time, rel_inc, maxima, minima] = peak_latency(x_z, tis_before, tis_after);
    
    dt = 10;
    x1 = 0.5e3 / dt;
    x2 = (tpuff + 10e3) / dt;
    
    [ai, bi] = sort(rel_inc, 'ascend'); %descend for PM ascend for NM
    bi(isnan(ai)) = [];
    
    figure('Position', [100, 100, 800, 500])
    hold on
    
    low_latency = 600;
    mid_latency = 1200;
    colors = { [0.4, 0.75, 0.6], ...
               [0.95, 0.8, 0.4], ...
               [0.9, 0.5, 0.5] };
    
    for ij = 1:size(bi, 2)
        ii = bi(ij);
        xmx = max(x_z(ii, x1:x2));
        plot([tpuff, tpuff], [0, xmx], 'k--', 'LineWidth', 2)
        plot([tpuff, tpuff] + t_after, [0, xmx], 'k-', 'LineWidth', 2)
        
        if peak_time(ii) <= low_latency
            color = colors{1};
        elseif peak_time(ii) <= mid_latency
            color = colors{2};
        else
            color = colors{3};
        end
        
        plot(T, x_z(ii, :), 'Color', color, 'LineWidth', 1)
    end
    
    xlim([T(x1), T(x2)])
    set(gca, 'XColor', 'none', 'YColor', 'none', 'Box', 'off')
    xlabel('')
    ylabel('')
    set(gca, 'LineWidth', 1, 'FontSize', 20, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01])
    avg_x = nanmean(x_z(bi(:), :), 1);
    plot(T, avg_x, 'k-', 'LineWidth', 2)
    
    low_count = sum(peak_time <= low_latency);
    mid_count = sum(peak_time > low_latency & peak_time <= mid_latency);
    high_count = sum(peak_time > mid_latency);
    
    x_start = T(x1) + 1000;
    y_start = max(ylim) - 0.2 * range(ylim);
    x_length = 1000;
    y_length = 1;
    line([x_start, x_start + x_length], [y_start, y_start], 'Color', 'k', 'LineWidth', 3);
    line([x_start, x_start], [y_start, y_start + y_length], 'Color', 'k', 'LineWidth', 3);
    text(x_start + x_length / 2, y_start - 0.05 * range(ylim), '1s', 'HorizontalAlignment', 'center', 'FontSize', 17);
    text(x_start - 0.03 * range(xlim), y_start + y_length / 2, '\DeltaF/F = 1', 'HorizontalAlignment', 'right', 'FontSize', 17);
    
    mfbFolderPath = 'X:\MFB';
    currentDate = datestr(now, 'yyyy-mm-dd');
    folderPath = fullfile(mfbFolderPath, 'Figures', 'Supp.Figures/S5', currentDate);
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end
    
    fileName = ['peak_latency_' char(file)];
    fullFilePathPDF = fullfile(folderPath, [fileName, '.pdf']);
    exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
    
    % save low latency index and time info
    low_latency_idx = peak_time <= low_latency;
    low_latency_t2p = peak_time(low_latency_idx);
    low_latency_t2p_all = [low_latency_t2p_all low_latency_t2p];
end

%% S5.a
figure('Position', [50, -100, 500, 400])

low_latency = 600;
mid_latency = 1200;
colors = { [0.4, 0.75, 0.6], 
           [0.95, 0.8, 0.4],
           [0.9, 0.5, 0.5] };

for i = 1:4
    subplot(2, 2, i)
    hold on
    
    for ij = 1:size(bi, 2)
        ii = bi(ij);
        if i == 1 && peak_time(ii) <= low_latency
            plot(T, x_z(ii, :), 'Color', colors{1}, 'LineWidth', 1)
        elseif i == 2 && peak_time(ii) > low_latency && peak_time(ii) <= mid_latency
            plot(T, x_z(ii, :), 'Color', colors{2}, 'LineWidth', 1)
        elseif i == 3 && peak_time(ii) > mid_latency
            plot(T, x_z(ii, :), 'Color', colors{3}, 'LineWidth', 1)
        elseif i == 4
            avg_x = nanmean(x_z(bi(1:length(bi)), :), 1);
            plot(T, avg_x, 'k-', 'LineWidth', 2)
        end
    end
    
    plot([tpuff, tpuff], [-10 10], 'k--')
    plot([tpuff + t_after, tpuff + t_after], [-10 10], 'k-')
    
    xlim([T(x1), T(x2)])
    ylim([-10 10])
    
    set(gca, 'XColor', 'none', 'YColor', 'none', 'LineWidth', 1, ...
        'FontSize', 12, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01])
    
    if i == 1
        x_start = T(x1) + 1500;
        y_start = -9;
        x_length = 1000;
        y_length = 2;
        line([x_start, x_start + x_length], [y_start, y_start], 'Color', 'k', 'LineWidth', 2);
        line([x_start, x_start], [y_start, y_start + y_length], 'Color', 'k', 'LineWidth', 2);
        text(x_start + x_length / 2, y_start - 2, '1s', 'HorizontalAlignment', 'center', 'FontSize', 13);
        text(x_start - 0.08 * range(xlim), y_start + y_length / 2, '2 \DeltaF/F', 'HorizontalAlignment', 'right', 'FontSize', 13);
    end
    
    titles = {'Low (<600ms)', 'Medium (600-1200ms)', 'High (>1200ms)', 'Average'};
    title(titles{i})
end

fileName = ['peak_latency_separate'];
fullFilePathPDF = fullfile(folderPath, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% S5.f
time_to_peak = struct('Nigel', [], 'Bernie', [], 'Bill', [], 'Jeremy', []);

Nigel = {'171212_16_19_37'};
Bernie = {'191209_13_44_12', '191209_14_04_14', '191209_14_32_39', '191209_15_01_22', '191209_14_18_13', '191209_14_46_58'};
Bill = {'200130_13_21_13 FunctAcq', '200130_13_36_14 FunctAcq', '200130_13_49_09 FunctAcq', '200130_14_15_24 FunctAcq', '200130_14_29_30 FunctAcq'};
Jeremy = {'191018_13_39_41', '191018_13_56_55', '191018_14_11_33'};

for file_i = 1:length(folder_names)
    file = folder_names{file_i};
    quickAnalysis;
    
    tpuff = 5e3;
    t_before = 0.5 * 1e3;
    t_after = 2 * 1e3;
    
    Tb1 = tpuff - t_before;
    Tb2 = tpuff;
    Ta1 = tpuff;
    Ta2 = tpuff + t_after;
    
    tis_before = find((T >= Tb1) .* (T < Tb2));
    tis_after = find((T >= Ta1) .* (T <= Ta2));
    
    spd_th = .1;
    spdm = nanmean(speed(:, tis_before), 2);
    ROI_dff_mfa = dff_gz(:,:,:);
    ['fraction filtered: ', num2str(sum(spdm(:) < spd_th) / numel(spdm) * 100), '%']
    x = squeeze(nanmean(ROI_dff_mfa(:, spdm < spd_th, :), 2));
    [peak_time, ~, ~] = peak_latency(x, tis_before, tis_after);
    
    if ismember(file, Nigel)
        time_to_peak.Nigel = [time_to_peak.Nigel peak_time];
    elseif ismember(file, Bernie)
        time_to_peak.Bernie = [time_to_peak.Bernie peak_time];
    elseif ismember(file, Bill)
        time_to_peak.Bill = [time_to_peak.Bill peak_time];
    elseif ismember(file, Jeremy)
        time_to_peak.Jeremy = [time_to_peak.Jeremy peak_time];
    end
end

mean_peak_time_nigel = nanmean(time_to_peak.Nigel);
mean_peak_time_bernie = nanmean(time_to_peak.Bernie);
mean_peak_time_bill = nanmean(time_to_peak.Bill);
mean_peak_time_jeremy = nanmean(time_to_peak.Jeremy);

std_peak_time_nigel = nanstd(time_to_peak.Nigel);
std_peak_time_bernie = nanstd(time_to_peak.Bernie);
std_peak_time_bill = nanstd(time_to_peak.Bill);
std_peak_time_jeremy = nanstd(time_to_peak.Jeremy);

fprintf('Nigel mean time to peak: %.2f ms, std: %.2f ms\n', mean_peak_time_nigel, std_peak_time_nigel);
fprintf('Bernie mean time to peak: %.2f ms, std: %.2f ms\n', mean_peak_time_bernie, std_peak_time_bernie);
fprintf('Jeremy mean time to peak: %.2f ms, std: %.2f ms\n', mean_peak_time_jeremy, std_peak_time_jeremy);
fprintf('Bill mean time to peak: %.2f ms, std: %.2f ms\n', mean_peak_time_bill, std_peak_time_bill);

%% boxplot
figure;
group_labels = [repmat({'Animal 1'}, length(time_to_peak.Nigel), 1); ...
                repmat({'Animal 2'}, length(time_to_peak.Bernie), 1); ...
                repmat({'Animal 3'}, length(time_to_peak.Jeremy), 1); ...
                repmat({'Animal 4'}, length(time_to_peak.Bill), 1)];
h = boxplot([time_to_peak.Nigel'; time_to_peak.Bernie'; time_to_peak.Jeremy'; time_to_peak.Bill'], group_labels, 'colors', 'k');
ylabel('Time to Peak (ms)');
yticks([0 1000 2000]);

colors = [1, 0.4, 0.4;
          0.7, 0.1, 0.7;
          1, 0.7, 0;
          0.1, 0.5, 0.8];

hBox = findobj(gca, 'Tag', 'Box');
hBox = flipud(hBox);
for j = 1:length(hBox)
    patch(get(hBox(j), 'XData'), get(hBox(j), 'YData'), colors(j, :), ...
          'FaceAlpha', 0.5, 'EdgeColor', colors(j, :));
end

set(gca, 'LineWidth', 1, 'FontSize', 17, ...
    'FontName', 'Helvetica', ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'TickLength', [.02 .02], ...
    'ZMinorTick', 'off');

mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
folderPath = fullfile(mfbFolderPath, 'Figures', 'Supp.Figures/S5/', currentDate);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

fileName = ['Puff_latency_box'];
fullFilePathPDF = fullfile(folderPath, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
