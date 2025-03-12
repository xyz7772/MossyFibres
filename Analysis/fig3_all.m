clear all;close all;clc
%% fig 3a heat map of activity
concat_folder = 'X:\MFB\Concatenation data'; % path to concatenated data folder
load([concat_folder '\concat_Animal2.mat'])

r2w = [linspace(1, 1, 128)', linspace(0, 1, 128)', linspace(0, 1, 128)'];
w2b = [linspace(1, 0, 128)', linspace(1, 0, 128)', linspace(1, 1, 128)'];
r2b = [r2w; w2b];

mycm = flipud(r2b);
suffix = 'changeSorted';
Nmf = size(dff_r,1);
Ntot = Nmf;
dt = 10;
figure('Position', [30 20 2000 700])
sbplt1 = [1:8];
sbplt2 = [9];
sbplt3 = [10];
sbp_no = 10;    
dff_rz = (dff_r - nanmean(dff_r')') ./ nanstd(dff_r')';  
ccs = corr(L_state', dff_rz', 'rows', 'complete');
[a,b]=sort(ccs);
zz = dff_rz(b,:);

% clean
nans_t = any(isnan(dff_rz), 1);

dff_rz(:, nans_t) = [];
MI_whisker_r(:, nans_t) = [];
L_state(:, nans_t) = [];
A_state(:, nans_t) = [];
Q_state(:, nans_t) = [];
MI_wheel_r(:, nans_t) = [];


zz = dff_rz(b,:);
sorted_zz = sort(zz(~isnan(zz)), 'ascend');
n = numel(sorted_zz);
lower_i = max(floor(n * 0.01), 1);
upper_i = max(floor(n * 0.99), 1);
lv = round(sorted_zz(lower_i), 1);
uv = round(sorted_zz(upper_i), 1);

xlm = [0,size(zz,2)];
subplot(sbp_no,1,sbplt1); hold on
imagesc(zz, [lv,uv]);

ylb = ylabel('Mossy fiber axon number');
set(ylb, 'Units', 'Normalized', 'Position', [-.04, 0.5, 0]);
yticks([1, Ntot])
ylim([0+.5,Ntot+.5])
xticklabels([]);
xlim(xlm)
cb = colorbar();
set(cb,'position',[0.1 .82 .02 .05],'Ticks', [lv, uv],'FontSize', 11)
titleStr = {'Z-scored', '\DeltaF/F'};
title(cb, titleStr, 'FontSize', 10, 'FontWeight', 'normal');

pos = get(gca, 'Position');
pos(1) = 0.15; pos(3) = 0.75;
pos(2) = .26; pos(4) = .7;
set(gca, 'Position', pos, 'LineWidth', 1, 'FontSize', 17, 'XColor', 'none', 'YDir', 'normal','TickDir', 'out')
time = (1:length(L_state)) ; % 100 = 1s

% Subplot for MI_whisker_r
subplot(sbp_no,1,sbplt3) 
hold on
wsk_r = MI_whisker_r;
ylb = ylabel({'Whisk'});
set(ylb, 'Units', 'Normalized', 'Position', [-0.04, 0.5, 0]);
pos = get(gca, 'Position');
pos(1) = 0.15; pos(3) = 0.75;
pos(2) = .04; pos(4) = .075;

set(gca, 'Position', pos,'LineWidth', 1, 'FontSize', 17, 'XColor', 'none','TickDir', 'out')
xlim(xlm)
ylim([-0.05 0.2])
AS_regions = A_state == 1;
AS_start = find(diff([0; AS_regions(:)]) == 1);
AS_end = find(diff([AS_regions(:); 0]) == -1);
for i = 1:length(AS_start)
    patch([time(AS_start(i)), time(AS_end(i)), time(AS_end(i)), time(AS_start(i))], ...
          [-0.05, -0.05, 0.8, 0.8], ...
          [255 205 255]/255, 'EdgeColor', 'none');
end
QW_regions = Q_state == 1;
QW_start = find(diff([0; QW_regions(:)]) == 1);
QW_end = find(diff([QW_regions(:); 0]) == -1);
for i = 1:length(QW_start)
    patch([time(QW_start(i)), time(QW_end(i)), time(QW_end(i)), time(QW_start(i))], ...
          [-0.05, -0.05, 0.8, 0.8], ...
          [192 255 255]/255, 'EdgeColor', 'none');
end

plot(time, wsk_r, 'LineWidth',1, 'color', [255,163,26]/255)

hold on;
x0 = max(xlim)/4.5; 
y0 = min(wsk_r) ; 
plot(x0-[0,6000], [y0-0.05,y0-0.05], 'k-', 'LineWidth',3);
text(x0-3000, y0-0.15, '60 s', 'FontSize', 15, 'HorizontalAlignment', 'center'); 

% Subplot for MI_wheel_r
subplot(sbp_no,1,sbplt2)
hold on
ylb = ylabel({'Loco'});
set(ylb, 'Units', 'Normalized', 'Position', [-.04, 0.5, 0]);

pos = get(gca, 'Position');
pos(1) = 0.15; pos(3) = 0.75;
pos(2) = .16; pos(4) = .075;
set(gca, 'Position', pos)
xlim(xlm)
ylim([-0.05 0.8])
yticks([0 0.4 0.8])
for i = 1:length(AS_start)
    patch([time(AS_start(i)), time(AS_end(i)), time(AS_end(i)), time(AS_start(i))], ...
          [-0.05, -0.05, 0.8, 0.8], ...
          [255 205 255]/255, 'EdgeColor', 'none');
end

for i = 1:length(QW_start)
    patch([time(QW_start(i)), time(QW_end(i)), time(QW_end(i)), time(QW_start(i))], ...
          [-0.05, -0.05, 0.8, 0.8], ...
          [192 255 255]/255, 'EdgeColor', 'none');
end
plot(time, MI_wheel_r, 'color', [173,210,157]/255, 'LineWidth',1);
set(gca, 'LineWidth', 1, 'FontSize', 17, 'XColor', 'none')


currentDate = datestr(now, 'yyyy-mm-dd');
savepath2 = fullfile(savepath, 'Figures', 'Figure3', currentDate);
if ~exist(savepath2, 'dir')
    mkdir(savepath2);
end
fileName = ['HeatMap__' suffix];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


%% fig 3b pc


files = dir(fullfile(concat_folder, 'concat_*.mat'));
idx = 5;
filePath = fullfile(concat_folder, files(idx).name)
load(filePath);

dff_rz = (dff_r - nanmean(dff_r')') ./ nanstd(dff_r')'; 
zz = dff_rz;
[cc_wL, zc_wL] = bootstrap_cc(zz', L_state', 100, 200);
sig_level = 2;
pm = zc_wL > sig_level;
pm_val = zc_wL(pm);
[~, pm_ids] = find(pm == 1);
nm = zc_wL < -sig_level;
nm_val = zc_wL(nm);
[~, nm_ids] = find(nm == 1);
ns = abs(zc_wL) < sig_level;
ns_val = zc_wL(ns);
[~, ns_ids] = find(ns == 1);

dt = 10;
nsh = 10;
fs = 15;
% intrepolate nans
zzf = fillmissing(zz, 'linear', 2);

if idx ==1
    zzf = zzf(:,100:47120);
    A_state = A_state(100:47120);
    Q_state = Q_state(100:47120);
end

[coeff, score, LATENT, TSQUARED, EXPLAINED] = pca(zzf', 'NumComponents',20);
dFF_L = zzf(:, A_state==1);
dFF_Q = zzf(:, Q_state==1);              

mu_L = nanmean(dFF_L,2);
mu_Q = nanmean(dFF_Q,2);

w = mu_L - mu_Q;
w = w / norm(w);

A_or_Q = w'*(zzf - nanmean(zzf,2));
A_or_Q = (A_or_Q - min(A_or_Q)) / (max(A_or_Q) - min(A_or_Q)); 
%A_or_Q = (A_or_Q - nanmean(A_or_Q)) ./ nanstd(A_or_Q);

c = nan(size(A_or_Q,2),3);
ls = 1:nsh:length(A_state);% lower sampling
score = score(ls, :);
A_or_Q = A_or_Q(ls);
A_state = A_state(ls);
Q_state = Q_state(ls);

figure('Position',[100,100,350,300]); hold on;
for k = 1:1:length(score)-1
    c_ = A_or_Q(k);
    if c_ < 0
        c_ = 0;
    elseif c_ > 1
       c_ = 1;
    end
    c(k,:) = c_*[1,0,1] + (1-c_)*[0,1,1];

    if ~isnan(c(k,:))
        plot3(score(k:k+1,1),score(k:k+1,2),score(k:k+1,3),'-','Color',c(k,:),'LineWidth', 1);
    end
end

if idx == 2
    view(60,27)
elseif idx==3  
    view(-114,51)
elseif idx==4
    view(160,-78)
elseif idx== 5 || idx==6||idx ==7
    view(74,68)
elseif idx==1
    view(-34,21)
elseif idx==8
    view(115,-58)
elseif idx==9
    view(-72,26)
elseif idx==10
    view(-73,-51)
end

set(gca,'FontSize',15,'XTick',[],'YTick',[],'ZTick',[])
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')

set(gca, 'LineWidth', 1, 'FontSize', fs, 'TickDir', 'out', 'TickLength',[.025,.025])

fileName = ['pc_activitiy_proj_' files(idx).name(1:end-4)];
fullFilePathPDF = fullfile(savepath, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

    
%% fig 3c manifold_distance
%Select
folder_names = {
   '171212_16_19_37';
   '191018_13_39_41';
   '191018_13_56_55';
   '191018_14_30_00';
   '191018_14_11_33';
   '191209_13_44_12';
   '191209_14_04_14';
   '191209_14_32_39';
   '191209_15_01_22';
   '191209_14_18_13';
   '191209_14_46_58';
   '200130_13_21_13 FunctAcq';
   '200130_13_36_14 FunctAcq';
   '200130_13_49_09 FunctAcq';
   %'200130_14_02_12 FunctAcq';
   '200130_14_15_24 FunctAcq';
   '200130_14_29_30 FunctAcq';
    };

num_PCs = 2;
angle_LQ_val = [];
angle_LQ_sh = [];
angle_LQ_sh_bl = [];
Ns = [];

dist_L = [];
dist_Q = [];
dist_LQ = [];

ijks = [];

for ii = 1:length(folder_names)

    file=char(folder_names(ii));
    quickAnalysis;
    block_size = 100;     
    dt = 10;
    
    [N,T] = size(dff_r);
    Ns = [Ns, N];
    
    dFF_L = dff_rz(:, A_state==1);
    dFF_Q = dff_rz(:, Q_state==1);  
    dsr = 100;
    xL = dFF_L(:,1:dsr:end)';
    xQ = dFF_Q(:,1:dsr:end)';
    
    if size(dFF_L,2)<100
        xL = dFF_L(:,1:end)';
    end

    x = pdist2(xL, xL);
    dist_L = [dist_L, nanmean(x(:))];
    
    x = pdist2(xQ, xQ);
    dist_Q = [dist_Q, nanmean(x(:))];
    
    x = pdist2(xL, xQ);
    dist_LQ = [dist_LQ, nanmean(x(:))];

end

figure('Position', [100,100,350,350])
avg_L = nanmean(dist_L);
sem_L = nanstd(dist_L) / sqrt(sum(~isnan(dist_L)));
avg_Q = nanmean(dist_Q);
sem_Q = nanstd(dist_Q) / sqrt(sum(~isnan(dist_Q)));
avg_LQ = nanmean(dist_LQ);
sem_LQ = nanstd(dist_LQ) / sqrt(sum(~isnan(dist_LQ)));
bar_positions = 1:3;

b = bar(bar_positions, [avg_Q, avg_L, avg_LQ], 'FaceColor', 'none', 'EdgeColor', 'black','LineWidth', 1);

hold on;

for iii = 1:length(dist_L)
    plot(bar_positions, [dist_L(iii), dist_Q(iii), dist_LQ(iii)], 'o-', 'Color', [0.7, 0.7, 0.7]);
    hold on
end

scatter(ones(size(dist_L)) * bar_positions(1), dist_L, 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'magenta','LineWidth', 1.5);
scatter(ones(size(dist_Q)) * bar_positions(2), dist_Q, 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'cyan','LineWidth', 1.5);
scatter(ones(size(dist_LQ)) * bar_positions(3), dist_LQ, 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.5,0.5,0.5],'LineWidth', 1.5);

hErrorbar = errorbar(bar_positions, [avg_L, avg_Q, avg_LQ], [sem_L, sem_Q, sem_LQ], 'k', 'linestyle', 'none');
uistack(hErrorbar, 'top');
fs = 15;
xlim([0,4])
ylim([10 30])
yticks([10 20 30])
box('off')
ylabel('Avg. distance')
xticks([1,2,3])
xticklabels({'AS', 'QW', 'AS-QW'})
xlabel({'Within/between manifolds'})
set(gca, 'LineWidth', 1, 'FontSize', fs, 'TickDir', 'out', 'TickLength',[.025,.025])

box("off")

[p_LQ, h_LQ] = signrank(dist_L, dist_Q);
[p_L_LQ, h_L_LQ] = signrank(dist_L, dist_LQ);
[p_Q_LQ, h_Q_LQ] = signrank(dist_Q, dist_LQ);

disp(['AA and QQ: p value = ', num2str(p_LQ)]);
disp(['AA and AQ: p value = ', num2str(p_L_LQ)]);
disp(['QQ and AQ: p value  = ', num2str(p_Q_LQ)]);


fileName = ['Manifold_distance_all'];
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% fig 3d Manifold Angle
num_PCs = 2;
angle_LQ_val = [];
angle_LQ_min = [];
angle_LQ_sh_bl_max = [];
angle_LQ_sh_bl_min = [];
Ns = [];
dt = 10;

for i = 1:length(folder_names)
    file = char(folder_names(i));
    quickAnalysis;
    [N, T] = size(dff_rz);
    Ns = [Ns, N];

    if all(L_state == 0)
        continue
    end

    dFF_L = dff_rz(:, A_state == 1);
    dFF_Q = dff_rz(:, Q_state == 1);
    clear L_state
    
    coeff_L = pca(dFF_L');
    coeff_Q = pca(dFF_Q');
    
    A_metric = eye(size(coeff_Q,1));
    [theta, ~, ~] = subspacea(coeff_Q(:,1:num_PCs), coeff_L(:,1:num_PCs), A_metric);
    angle_LQ_val = [angle_LQ_val, max(theta)];
    angle_LQ_min = [angle_LQ_min, min(theta)];

    ashes_max = [];
    ashes_min = [];
    N_itr = 20;
    for j = 1:N_itr
        bl_sz = 1000;
        dff_sh_bl = dff_r(:, randblock(1:floor(T / bl_sz) * bl_sz, [1, bl_sz]));
        assert(mod(size(dff_sh_bl,2), 2) == 0, 'Shuffled data must have even columns.');
        
        dff_1 = dff_sh_bl(:, 1:size(dff_sh_bl, 2)/2);
        dff_2 = dff_sh_bl(:, size(dff_sh_bl, 2)/2+1:end);
        
        coeff_1 = pca(dff_1');
        coeff_2 = pca(dff_2');
        
        [theta_sh, ~, ~] = subspacea(coeff_1(:,1:num_PCs), coeff_2(:,1:num_PCs), A_metric);
        ashes_max = [ashes_max, max(theta_sh)];
        ashes_min = [ashes_min, min(theta_sh)];
    end

    angle_LQ_sh_bl_max = [angle_LQ_sh_bl_max, mean(ashes_max, 'omitnan')];
    angle_LQ_sh_bl_min = [angle_LQ_sh_bl_min, mean(ashes_min, 'omitnan')];
end

[p_max, ~, ~] = signrank(angle_LQ_val, angle_LQ_sh_bl_max);
[p_min, ~, ~] = signrank(angle_LQ_min, angle_LQ_sh_bl_min);
disp(['Max Angle p = ', num2str(p_max)]);
disp(['Min Angle p = ', num2str(p_min)]);

figure('Position', [100, 100, 250, 300])
hold on

h1 = plot(nan, nan, 'o', 'LineWidth', 2, 'Color', [0,0.5,0.4], 'MarkerFaceColor', [0,0.5,0.4]);
h2 = plot(nan, nan, 'o', 'Color', [0.7 0.7 0.7], 'MarkerFaceColor', [0.7 0.7 0.7], 'LineWidth', 1);

for i = 1:length(angle_LQ_val)
    plot([1, 2], [angle_LQ_val(i), angle_LQ_sh_bl_max(i)], '-o', ...
        'LineWidth', 2, 'Color', [0,0.5,0.4], 'MarkerFaceColor', [0,0.5,0.4])
    plot([1, 2], [angle_LQ_min(i), angle_LQ_sh_bl_min(i)], ...
        'Color', [0.7 0.7 0.7], 'Marker', 'o', 'LineWidth', 1, ...
        'MarkerFaceColor', [0.7 0.7 0.7])
end

legend([h1, h2], {'LA', 'SA'}, 'Location', 'best','Box','off');


xlim([0.5, 2.5])
ylim([0, 2])
xticks([1, 2])
xticklabels({'QW-AS', 'Shuffle'})
xtickangle(45)
fs = 15;
set(gca, 'LineWidth', 1, 'FontSize', fs, 'TickDir', 'out', 'TickLength', [.025, .025])
ylabel('Angle (rad.)')

fileName = 'manifold_angle_all';
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% fig 3e manifold_angle_percentile PM/NM
num_PCs = 2;
Ns = [];
prcnts = 10:10:90;
angles_LQ = [];
angles_LQ_shuffle = [];

for file_i = 1:length(folder_names)
    file = char(folder_names(file_i));
    quickAnalysis;
    
    dt = 10;
    block_size = 100; 
    
    if ~all(L_state==0)
        [cc_wL, zc_wL] = bootstrap_cc(dff_rz', L_state', block_size, 100);
    else
        continue
    end

    [N,T] = size(dff_rz);
    Ns = [Ns, N];
    vals = prctile(abs(zc_wL), 100 - prcnts);

    angle_LQ_min_current = nan(1, length(prcnts));
    angle_LQ_sh = nan(1, length(prcnts));

    for i = 1:length(prcnts)
        ids = find(abs(zc_wL) < vals(i));
        
        dFF_L = dff_rz(ids, L_state==1);
        dFF_Q = dff_rz(ids, L_state==0);
        
        coeff_L = pca(dFF_L');
        coeff_Q = pca(dFF_Q');
        
        A_metric = eye(size(coeff_Q,1));
        [theta, ~, ~] = subspacea(coeff_Q(:,1:num_PCs), coeff_L(:,1:num_PCs), A_metric);
        
        angle_LQ_min_current(i) = max(theta);
        
        ashes = [];
        N_itr = 20;
        for kk = 1:N_itr
            bl_sz = 1000; % 10s
            dff_sh_bl = dff_r(:, randblock(1:floor(T/bl_sz)*bl_sz, [1, bl_sz]));

            dff_1 = dff_sh_bl(ids, 1:size(dff_sh_bl,2)/2);
            dff_2 = dff_sh_bl(ids, size(dff_sh_bl,2)/2+1:end);
            
            coeff_1 = pca(dff_1');
            coeff_2 = pca(dff_2');
            
            [theta_sh, ~, ~] = subspacea(coeff_1(:,1:num_PCs), coeff_2(:,1:num_PCs), A_metric);
            ashes = [ashes, max(theta_sh)];
        end
        angle_LQ_sh(i) = nanmean(ashes);
    end
    
    angles_LQ = [angles_LQ; angle_LQ_min_current];
    angles_LQ_shuffle = [angles_LQ_shuffle; angle_LQ_sh];
end

figure('Position', [100,100,400,300])
hold on

meanangle = nanmean(angles_LQ, 1);
meanshuffle = nanmean(angles_LQ_shuffle, 1);
SEM1 = nanstd(angles_LQ, 0, 1) / sqrt(size(angles_LQ, 1));
SEM2 = nanstd(angles_LQ_shuffle, 0, 1) / sqrt(size(angles_LQ_shuffle, 1));

errorbar(prcnts, meanangle, SEM1, 'k-', 'LineWidth', 2)
errorbar(prcnts, meanshuffle, SEM2, '-', 'LineWidth', 2, 'Color', [0.7, 0.7, 0.7])

legend('avg.', 'shuffle', 'Location', 'Best')
legend boxoff
ylim([0, 1.7])
fs = 15;
set(gca, 'LineWidth', 1, 'FontSize', fs, 'TickDir', 'out', 'TickLength', [.025, .025])

ylabel('Angle (rad.)')
xlabel('PM/NM excl. percentile (%)')

fileName = 'manifold_angle_percentiles';
fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% fig 3 fg
num_PCs = 2;
Ns = [];
prcnts = 10:10:90;
angles_LQ = [];
angles_LQ_shuffle = [];
sig_level = 2;

for md = 1:2
    for file_i = 1:length(folder_names)
        file = folder_names{file_i};
        quickAnalysis;
        dt = 10;
        block_size = 100;

        if ~all(L_state == 0)
            [cc_wL, zc_wL] = bootstrap_cc(dff_rz', L_state', block_size, 200);
        else
            continue;
        end

        [N, T] = size(dff_rz);
        Ns = [Ns, N];

        pm = zc_wL > sig_level;
        [~, pm_ids] = find(pm == 1);
        nm = zc_wL < -sig_level;
        [~, nm_ids] = find(nm == 1);
        vals_pm = prctile(zc_wL(pm_ids), 100 - prcnts);
        vals_nm = prctile(zc_wL(nm_ids), 100 - prcnts);

        angle_LQ_val = nan(1, length(prcnts));
        angle_LQ_shuffle = nan(1, length(prcnts));

        for i = 1:length(prcnts)
            if md == 1
                suffix = 'PM';
                ids = find(abs(zc_wL(pm_ids)) < abs(vals_pm(i)));
                zz = dff_rz(pm_ids, :);
            elseif md == 2
                suffix = 'NM';
                ids = find(zc_wL(nm_ids) < vals_nm(i));
                zz = dff_rz(nm_ids, :);
            end

            dFF_L = zz(ids, A_state == 1);
            dFF_Q = zz(ids, Q_state == 1);
            coeff_L = pca(dFF_L');
            coeff_Q = pca(dFF_Q');
            if sum(all(coeff_Q == 1)) == 1 || sum(all(coeff_L == 1)) == 1
                break;
            end

            A_metric = eye(size(coeff_Q, 1));
            [theta, ~, ~] = subspacea(coeff_Q(:,1:num_PCs), coeff_L(:,1:num_PCs), A_metric);
            angle_LQ_val(i) = max(theta);

            ashes = [];
            N_itr = 20;
            for kk = 1:N_itr
                bl_sz = 1000;
                dff_sh_bl = zz(:, randblock(1:floor(T/bl_sz)*bl_sz, [1, bl_sz]));

                dff_1 = dff_sh_bl(ids, 1:size(dff_sh_bl,2)/2);
                dff_2 = dff_sh_bl(ids, size(dff_sh_bl,2)/2+1:end);

                coeff_1 = pca(dff_1');
                coeff_2 = pca(dff_2');

                A_metric = eye(size(coeff_1, 1));
                [theta_sh, ~, ~] = subspacea(coeff_1(:,1:num_PCs), coeff_2(:,1:num_PCs), A_metric);
                ashes = [ashes, max(theta_sh)];
            end

            angle_LQ_shuffle(i) = nanmean(ashes);
        end

        angles_LQ = [angles_LQ; angle_LQ_val];
        angles_LQ_shuffle = [angles_LQ_shuffle; angle_LQ_shuffle];
    end
   
    figure('Position', [100, 100, 400, 300]);
    hold on;

    meanangle = nanmean(angles_LQ, 1);
    meanshuffle = nanmean(angles_LQ_shuffle, 1);
    SEM1 = nanstd(angles_LQ, 0, 1) / sqrt(size(angles_LQ, 1));
    SEM2 = nanstd(angles_LQ_shuffle, 0, 1) / sqrt(size(angles_LQ_shuffle, 1));

    errorbar(prcnts, meanangle, SEM1, 'k-', 'LineWidth', 2);
    errorbar(prcnts, meanshuffle, SEM2, '-', 'LineWidth', 2, 'Color', [0.7, 0.7, 0.7]);

    legend('avg.', 'shuffle', 'Location', 'Best');
    legend boxoff;
    ylim([0, 1.7]);
    xlim([0, 100]);

    fs = 15;
    set(gca, 'LineWidth', 1, 'FontSize', fs, 'TickDir', 'out', 'TickLength', [0.025, 0.025]);

    ylabel('Angle (rad.)');
    xlabel([suffix, ' excl. percentile (%)']);

    fileName = ['manifold_angle_percentiles_', suffix];
    fullFilePathPDF = fullfile(savepath2, [fileName,'.pdf']);
    exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');
end

