%% load data
file ='200130_13_36_14 FunctAcq';
quickAnalysis
dff0 = ROI_dff_all;
Nmf0 = size(dff0,1);
dff0_r = reshape(permute(dff0,[1,3,2]), Nmf0, []);
dff0_rz = (dff0_r - nanmean(dff0_r')') ./ nanstd(dff0_r')';

%% S2a
cc0 = corr(dff0_r', 'rows', 'complete','type','Spearman'); 
cc0(eye(length(cc0))==1)=nan;
cc_th = prctile(cc0(~isnan(cc0)), 99); % 99th percentile

histogram(cc0(~isnan(cc0)),'FaceColor', [0 0 0],'BinWidth',0.02,'normalization','count')
xline(cc_th, '--r', 'LineWidth', 2, 'Label', '99th percentile');
xlabel('Spearman CC');
ylabel('Count');
yticks([0 6000 12000])
set(gca, 'LineWidth', 1, 'FontSize', 18, ...
    'FontName'   , 'Helvetica', ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'ZMinorTick'  , 'off'  ...
    )

currentDate = datestr(now, 'yyyy-mm-dd');
savepath2 = fullfile(savepath, 'Figures', 'Supp.Figures/S2/', currentDate);
if ~exist(savepath2, 'dir')
    mkdir(savepath2);
end
fileName = ['example correlation distribution'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

groups_cc = {};
pregroups1 = {};

for i = 1:Nmf0
    pregroups1{i} = i;

    for j = 1:Nmf0
        if i ~= j && cc0(i, j) > cc_th
            pregroups1{i} = [pregroups1{i}, j]; 
        end
    end
end

for i = 1:Nmf0
    if numel(pregroups1{i}) <= 2 
        groups_cc{i} = pregroups1{i};
        continue;
    end

    current_group = pregroups1{i};
    primary = current_group(1);
    cc_members = primary;

    for j = 2:length(current_group)
        for k = j+1:length(current_group)
            if cc0(current_group(j), current_group(k)) < cc_th
                if cc0(primary, current_group(j)) >= cc0(primary, current_group(k)) 
                    current_group(k) = [];
                else
                    current_group(j) = [];
                end
                break; 
            end
        end
        if numel(current_group) < 3
            break; 
        end
    end
    groups_cc{i} = current_group;
end


%% LDR
load(['X:\MFB\Smooth_data\', file, '.mat']);
load(['X:\MFB\LDR\', file, '.mat']);

ldr_th = 1.5; 
pregroups2 = cell(1, Nmf0);

for i = 1:Nmf0
    primary = i;
    pregroups2{i} = primary;
    for j = 2:length(groups_cc{i})
        if ldr(primary, groups_cc{i}(j)) <= ldr_th
            pregroups2{i} = [pregroups2{i}, groups_cc{i}(j)];
        end
    end
end

groups_ldr = cell(1, Nmf0);

for i = 1:Nmf0
    if numel(pregroups2{i}) <= 1
        groups_ldr{i} = pregroups2{i};
        continue;
    end

    current_group = pregroups2{i};
    primary = current_group(1); 
    ldr_members = primary;

    for j = 2:length(current_group)
        isgroup = true;
        for k = 2:length(current_group)
            if j ~= k && ldr(current_group(j), current_group(k)) > ldr_th
                isgroup = false;
                if ldr(primary, current_group(j)) < ldr(primary, current_group(k))
                    ldr_members = [ldr_members, current_group(j)];
                else
                    ldr_members = [ldr_members, current_group(k)];
                end
                break;
            end
        end
    
        if isgroup && ~ismember(current_group(j), ldr_members)
            ldr_members = [ldr_members, current_group(j)];
        end
    end
    groups_ldr{i} = ldr_members; 
end

plot(cc0,ldr,'o','color',[0.4 0.4 0.4])
yline(ldr_th, '--r', 'LineWidth', 2, 'Label', 'LDR = 1.5');
xlabel('Spearman CC');
ylabel('LDR');
ylim([0 10])
set(gca, 'LineWidth', 1, 'FontSize', 20, ...
    'FontName'   , 'Helvetica', ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'ZMinorTick'  , 'off'  ...
    )
fileName = ['LDR threshold'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


%% Supp 2c
m = [92,88]%[226 244];
fig = figure('visible', 'on', 'position', [100 100 830 150]);

x_limits = [1, size(dff0_r, 2)];
all_data = dff0_r(m,:);
y_limits = [min(all_data(:)), max(all_data(:))];

t = tiledlayout(length(m), 1, 'TileSpacing', 'none', 'Padding', 'compact');

for plot_i = 1:length(m)
    ax = nexttile;
    plot(dff0_r(m(plot_i), :), 'k');
    xlim(x_limits);
    ylim(y_limits);
    set(ax, 'LineWidth', 1, 'FontSize', 15, 'FontName', 'Helvetica',...
    'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02]);
    axis(ax, 'off');
end
linkaxes(findall(fig, 'type', 'axes'), 'xy');
ax = nexttile(length(m));
hold(ax, 'on');
x_start = x_limits(2) - 1100;
y_start = y_limits(1) - 1.4;
line(ax, [x_start, x_start+1000], [y_start, y_start],...
'Color', 'k', 'LineWidth', 2, 'Clipping', 'off');
line(ax, [x_start, x_start], [y_start, y_start+1],...
'Color', 'k', 'LineWidth', 2, 'Clipping', 'off');
text(ax, x_start+500, y_start, '10 s',...
'FontSize', 12, 'HorizontalAlignment', 'center',...
'VerticalAlignment', 'top');
text(ax, x_start-150, y_start+0.2, '1 \DeltaF/F',...
'FontSize', 12, 'HorizontalAlignment', 'right',...
'VerticalAlignment', 'bottom');
hold(ax, 'off');
set(findall(fig, 'type', 'axes'), 'XLim', x_limits, 'YLim', y_limits);
fileName = 'example_with_scale_bar';
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
