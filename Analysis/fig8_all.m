clear all; clc

folder_names = {'171212_16_19_37', 
    'superBernie', 
    'superBill', 
    'superJeremy'};

folderPath = 'X:\MFB\MFB_AH_2023\Correlation_data';
c_file = dir(fullfile(folderPath, 'concat_*.mat'));
N_GC = 300;
N_MF = 2400;
dend_L = 21.7; % um
dt = 10;
ite_l = 10; 
th = 2;

prob_PM = []; dd_PM = [];
prob_NM = []; dd_NM = [];
prob_NS = []; dd_NS = [];

results_itr = struct('dist', [], 'num_PM', [], 'num_NM', [], 'num_NSM', [], ...
    'pm_out', [], 'nm_out', [], 'nsm_out', [], 'PPPP', [], 'PPPS', [], ...
    'PPPN', [], 'PPNN', [], 'PPNS', [], 'PPSS', [], 'PNNN', [], 'PNNS', [], ...
    'PNSS', [], 'PSSS', [], 'NNNN', [], 'NNNS', [], 'NNSS', [], 'NSSS', [], ...
    'SSSS', []);


for fi = 1:length(c_file)
    filePath = fullfile(folderPath, c_file(fi).name);
    load(filePath);

    if strcmp(cell2mat(extractBetween(c_file(fi).name, 'concat_', '.mat')), '171212_16_19_37')
        file ='171212_16_19_37';
        quickAnalysis;
        dff_r = dff0_r;
    else
        dff_r = reshape(permute(ROI_dff_All, [3,2,1]), [], size(ROI_dff_All,1))';
    end

    dff_rz = (dff_r - nanmean(dff_r')') ./ nanstd(dff_r')';
    dff_rz_0 = dff_rz; 

    [cc_wL, zc_wL] = bootstrap_cc(dff_rz', L_state', 100, 300);
    sig_lev = 2;

    pmids = find(zc_wL(:) > sig_lev);
    nmids = find(zc_wL(:) < -sig_lev);
    nsmids = find(abs(zc_wL(:)) < sig_lev);

    dff_dig = dff_rz_0;
    dff_dig(pmids, :) = 1;
    dff_dig(nmids, :) = -1;
    dff_dig(nsmids, :) = 0; 

    for ite = 1:ite_l
        ccs_GC = [];
        out_mod_all = [];
        unimf =[];
        pm_out = 0;
        nm_out = 0;
        nsm_out = 0;
        dist_all = zeros(1, N_GC);
        
        PPPP_count = 0;
        PPPN_count = 0;
        PPPS_count = 0;
        PPNN_count = 0;
        PPNS_count = 0;
        PPSS_count = 0;
        PNNN_count = 0;
        PNNS_count = 0;
        PNSS_count = 0;
        PSSS_count = 0;
        NNNN_count = 0;
        NNNS_count = 0;
        NNSS_count = 0;
        NSSS_count = 0;
        SSSS_count = 0;

        Xall = xyz(:, 1);
        Yall = xyz(:, 2);
        Zall = xyz(:, 3);

        x_GC = my_rand([1, N_GC], min(Xall(:)) + dend_L, max(Xall(:)) - dend_L);
        x_MF = my_rand([1, N_MF], min(Xall(:)), max(Xall(:)));

        y_GC = my_rand([1, N_GC], min(Yall(:)) + dend_L, max(Yall(:)) - dend_L);
        y_MF = my_rand([1, N_MF], min(Xall(:)), max(Xall(:)));

        z_GC = my_rand([1, N_GC], min(Zall(:)) + dend_L, max(Zall(:)) - dend_L);
        z_MF = my_rand([1, N_MF], min(Zall(:)), max(Zall(:)));

        xyz_GC = [x_GC; y_GC; z_GC]';

%         xyz = [x_MF;y_MF;z_MF]'; % simulated MF
%         pmr = length(pmids)/(length(pmids)+length(nmids)+length(nsmids))
%         nmr = length(nmids)/(length(pmids)+length(nmids)+length(nsmids))
%         smr = length(nsmids)/(length(pmids)+length(nmids)+length(nsmids))
%         num_pm = round(pmr * N_MF);
%         num_nm = round(nmr * N_MF);
%         num_sm = N_MF - num_pm - num_nm;
%         rng shuffle
%         rand_indices = randperm(N_MF);
%         pmids = rand_indices(1:num_pm);
%         nmids = rand_indices(num_pm+1:num_pm+num_nm);
%         nsmids = rand_indices(num_pm+num_nm+1:end);
%         dff_dig = zeros(N_MF, size(dff_rz_0, 2));
%         dff_dig(pmids, :) = 1;
%         dff_dig(nmids, :) = -1;
%         dff_dig(nsmids, :) = 0; 

        for igc = 1:N_GC
            d = xyz_GC(igc,:) - xyz(:,:);
            dd_PM = [dd_PM, sqrt(sum(d(pmids,:).^2,2))'];
            dd_NM = [dd_NM, sqrt(sum(d(nmids,:).^2,2))'];
            dd_NS = [dd_NS, sqrt(sum(d(nsmids,:).^2,2))'];
            dd2 = sqrt(sum(d.^2, 2));
            [a, b] = sort(dd2);

            idx = 1:4;
            while true
                [is_independent, detected_indices] = check_MF_independence(b, idx, group_ids);
                if is_independent
                    break;
                else
                    for j = 1:length(detected_indices)
                        idx(detected_indices(j)) = max(idx) + 1;
                    end
                end
            end

            dff_GC = nansum(dff_dig(b(idx)));
            mean_dis(igc) = nanmean(a(idx));

            if dff_GC > th
                ccs_GC(igc) = 1;
            elseif dff_GC < -th
                ccs_GC(igc) = -1;
            else
                ccs_GC(igc) = 0;
            end

            out_mod_all = [out_mod_all, ccs_GC(igc)];

            GC_ids = b(idx);
            unimf = [unimf GC_ids];
            in_pids = ismember(GC_ids, pmids);
            in_nids = ismember(GC_ids, nmids);
            in_nsids = ismember(GC_ids, nsmids);

            if all(in_pids)
                PPPP_count = PPPP_count + 1;
            elseif all(in_nids)
                NNNN_count = NNNN_count + 1;
            elseif all(in_nsids)
                SSSS_count = SSSS_count + 1;
            elseif sum(in_pids)==3 && sum(in_nids)==1
                PPPN_count =PPPN_count +1;
            elseif sum(in_pids)==3 && sum(in_nsids)==1
                PPPS_count =PPPS_count +1;
            elseif sum(in_pids)==2 && sum(in_nids)==2
                PPNN_count = PPNN_count+1;
            elseif sum(in_pids)==2 && sum(in_nids)==1 && sum(in_nsids)==1
                PPNS_count = PPNS_count+1;
            elseif sum(in_pids)==2 && sum(in_nsids)==2
                PPSS_count = PPSS_count+1;
            elseif sum(in_pids)==1 && sum(in_nids)==3
                PNNN_count = PNNN_count+1;
            elseif sum(in_pids)==1 && sum(in_nids)==2 && sum(in_nsids)==1
                PNNS_count = PNNS_count+1;
            elseif sum(in_pids)==1 && sum(in_nids)==1 && sum(in_nsids)==2
                PNSS_count = PNSS_count+1;
            elseif sum(in_pids)==1 && sum(in_nsids)==3
                PSSS_count = PSSS_count+1;
            elseif sum(in_nids)==3 && sum(in_nsids)==1
                NNNS_count = NNNS_count+1;
            elseif sum(in_nids)==2 && sum(in_nsids)==2
                NNSS_count = NNSS_count+1;
            elseif sum(in_nids)==1 && sum(in_nsids)==3
                NSSS_count = NSSS_count+1;
            end

            % expect output (should = Th=0)
            if sum(in_pids) - sum(in_nids) > 0 || sum(in_pids) + sum(in_nsids) == 4 
                pm_out = pm_out + 1;
            elseif sum(in_nids) - sum(in_pids) > 0 || sum(in_nids) + sum(in_nsids) == 4 
                nm_out = nm_out + 1;
            elseif sum(in_nids) == sum(in_pids) || sum(in_nsids) == 4 
                nsm_out = nsm_out + 1;
            end
        end

      % Store iteration results
        results_itr(fi).dist(:, ite) = mean(mean_dis);
        results_itr(fi).num_PM(ite) = sum(out_mod_all == 1);
        results_itr(fi).num_NM(ite) = sum(out_mod_all == -1);
        results_itr(fi).num_NSM(ite) = sum(out_mod_all == 0);
        results_itr(fi).pm_out(ite) = pm_out;
        results_itr(fi).nm_out(ite) = nm_out;
        results_itr(fi).nsm_out(ite) = nsm_out;

        results_itr(fi).PPPP(:, ite) = PPPP_count/N_GC;
        results_itr(fi).PPPN(:, ite) = PPPN_count/N_GC;
        results_itr(fi).PPPS(:, ite) = PPPS_count/N_GC;
        results_itr(fi).PPNN(:, ite) = PPNN_count/N_GC;
        results_itr(fi).PPNS(:, ite) = PPNS_count/N_GC;
        results_itr(fi).PPSS(:, ite) = PPSS_count/N_GC;
        results_itr(fi).PNNN(:, ite) = PNNN_count/N_GC;
        results_itr(fi).PNNS(:, ite) = PNNS_count/N_GC;
        results_itr(fi).PNSS(:, ite) = PNSS_count/N_GC;
        results_itr(fi).PSSS(:, ite) = PSSS_count/N_GC;
        results_itr(fi).NNNN(:, ite) = NNNN_count/N_GC;
        results_itr(fi).NNNS(:, ite) = NNNS_count/N_GC;
        results_itr(fi).NNSS(:, ite) = NNSS_count/N_GC;
        results_itr(fi).NSSS(:, ite) = NSSS_count/N_GC;
        results_itr(fi).SSSS(:, ite) = SSSS_count/N_GC;
    end
end


%%
filePath = fullfile('X:\MFB\MFB_AH_2023\Correlation_data', ['results_itr_th2.mat']);
save(filePath, 'results_itr');


%%
all_PPPP = []; all_PPPN = []; all_PPPS = [];
all_PPNN = []; all_PPNS = []; all_PPSS = [];
all_PNNN = []; all_PNNS = []; all_PNSS = [];
all_PSSS = []; all_NNNN = []; all_NNNS = [];
all_NNSS = []; all_NSSS = []; all_SSSS = [];

for i = 1:4

    strf = string(['X:\MFB\MFB_AH_2023\Correlation_data\results_itr_th' num2str(i-1) '.mat'])
    load(strf)
    all_PPPP = [all_PPPP; results_itr(i).PPPP(:)];
    all_PPPN = [all_PPPN; results_itr(i).PPPN(:)];
    all_PPPS = [all_PPPS; results_itr(i).PPPS(:)];
    all_PPNN = [all_PPNN; results_itr(i).PPNN(:)];
    all_PPNS = [all_PPNS; results_itr(i).PPNS(:)];
    all_PPSS = [all_PPSS; results_itr(i).PPSS(:)];
    all_PNNN = [all_PNNN; results_itr(i).PNNN(:)];
    all_PNNS = [all_PNNS; results_itr(i).PNNS(:)];
    all_PNSS = [all_PNSS; results_itr(i).PNSS(:)];
    all_PSSS = [all_PSSS; results_itr(i).PSSS(:)];
    all_NNNN = [all_NNNN; results_itr(i).NNNN(:)];
    all_NNNS = [all_NNNS; results_itr(i).NNNS(:)];
    all_NNSS = [all_NNSS; results_itr(i).NNSS(:)];
    all_NSSS = [all_NSSS; results_itr(i).NSSS(:)];
    all_SSSS = [all_SSSS; results_itr(i).SSSS(:)];
end

sem_PPPP = std(all_PPPP) / sqrt(length(all_PPPP));
sem_PPPN = std(all_PPPN) / sqrt(length(all_PPPN));
sem_PPPS = std(all_PPPS) / sqrt(length(all_PPPS));
sem_PPNN = std(all_PPNN) / sqrt(length(all_PPNN));
sem_PPNS = std(all_PPNS) / sqrt(length(all_PPNS));
sem_PPSS = std(all_PPSS) / sqrt(length(all_PPSS));
sem_PNNN = std(all_PNNN) / sqrt(length(all_PNNN));
sem_PNNS = std(all_PNNS) / sqrt(length(all_PNNS));
sem_PNSS = std(all_PNSS) / sqrt(length(all_PNSS));
sem_PSSS = std(all_PSSS) / sqrt(length(all_PSSS));
sem_NNNN = std(all_NNNN) / sqrt(length(all_NNNN));
sem_NNNS = std(all_NNNS) / sqrt(length(all_NNNS));
sem_NNSS = std(all_NNSS) / sqrt(length(all_NNSS));
sem_NSSS = std(all_NSSS) / sqrt(length(all_NSSS));
sem_SSSS = std(all_SSSS) / sqrt(length(all_SSSS));

mean_PPPP = mean(mean(all_PPPP));
mean_PPPN = mean(mean(all_PPPN));
mean_PPPS = mean(mean(all_PPPS));
mean_PPNN = mean(mean(all_PPNN));
mean_PPNS = mean(mean(all_PPNS));
mean_PPSS = mean(mean(all_PPSS));
mean_PNNN = mean(mean(all_PNNN));
mean_PNNS = mean(mean(all_PNNS));
mean_PNSS = mean(mean(all_PNSS));
mean_PSSS = mean(mean(all_PSSS));
mean_NNNN = mean(mean(all_NNNN));
mean_NNNS = mean(mean(all_NNNS));
mean_NNSS = mean(mean(all_NNSS));
mean_NSSS = mean(mean(all_NSSS));
mean_SSSS = mean(mean(all_SSSS));

patterns = {'NNNN', 'NNNS', 'PNNN', 'NNSS', 'PNNS', 'NSSS', 'PNSS', 'PPNN', 'SSSS', 'PSSS', 'PPNS', 'PPSS', 'PPPN', 'PPPS', 'PPPP'};
mean_values_map = containers.Map(patterns, [mean_NNNN, mean_NNNS, mean_PNNN, mean_NNSS, mean_PNNS, ...
    mean_NSSS, mean_PNSS, mean_PPNN, mean_SSSS, mean_PSSS, mean_PPNS, mean_PPSS, mean_PPPN, mean_PPPS, mean_PPPP]);
sem_values_map = containers.Map(patterns, [sem_NNNN, sem_NNNS, sem_PNNN, sem_NNSS, sem_PNNS, ...
    sem_NSSS, sem_PNSS, sem_PPNN, sem_SSSS, sem_PSSS, sem_PPNS, sem_PPSS, sem_PPPN, sem_PPPS, sem_PPPP]);

mean_values = cellfun(@(x) mean_values_map(x), patterns);
sem_values = cellfun(@(x) sem_values_map(x), patterns);

theoretical_values =[1.1,3.78, 4.03, 4.93,10.7, 2.92, 10,7.63,0.676,...
         3.34 ,15.2,8.32,9.69 ,11.3,6.36];

figure('Position', [100, 100, 850, 350]);
b = bar(mean_values * 100);

colors = zeros(length(patterns), 3);
for i = 1:length(patterns)
    numP = sum(patterns{i} == 'P');
    numN = sum(patterns{i} == 'N');
    numS = sum(patterns{i} == 'S');
    total = numP + numN + numS;

    rR = numP / total;
    bR = numN / total;
    gR = numS / total;

    colors(i, :) = [1, 0, 0] * rR + [0, 0, 1] * bR + [0.6, 0.6, 0.6] * gR;
end
b.CData = colors;
b.FaceColor = 'flat';
hold on;

x_pos = 1:length(mean_values);
errorbar(x_pos, mean_values * 100, sem_values * 100, 'k.', 'LineWidth', 1.5);

for i = 1:length(theoretical_values)
    y_value = theoretical_values(i);
    x_coords = [i - 0.2, i + 0.2, i];
    y_coords = [y_value, y_value, y_value - 0.4];
    fill(x_coords, y_coords, [255, 178, 102]/255, 'EdgeColor', 'none');
end

ylabel('Fraction (%)');
yticks([0 10 20]);
xticklabels(patterns);
set(gca, 'LineWidth', 1, 'FontSize', 15, 'TickDir', 'out', 'Box', 'off');

mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
folderPath = fullfile(mfbFolderPath, 'Figures', 'Figure8', currentDate);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

fileName = ['GC_INPUT.png'];
fullFilePath = fullfile(folderPath, fileName);
% save
print(fullFilePath, '-dpng', '-r300');

%% fraction - output

colors = {'r', 'b', 'k'};
outdata = zeros(3,4);

for i = 0:3
    file = ['X:\MFB\MFB_AH_2023\Correlation_data\results_itr_th' num2str(i) '.mat'];
    load(file);
    
    all_num_PM  = cat(1, results_itr(:).num_PM);
    all_num_NM  = cat(1, results_itr(:).num_NM);
    all_num_NSM  = cat(1, results_itr(:).num_NSM);

    PM_out_p = all_num_PM / mean(all_num_PM + all_num_NM + all_num_NSM) * 100;
    NM_out_p = all_num_NM / mean(all_num_PM + all_num_NM + all_num_NSM) * 100;
    NSM_out_p = all_num_NSM / mean(all_num_PM + all_num_NM +  all_num_NSM) * 100;

    outdata(1,i+1) = [mean(mean(PM_out_p,2))];
    outdata(2,i+1) = [mean(mean(NM_out_p,2))];
    outdata(3,i+1) = [mean(mean(NSM_out_p,2))];
end

figure('Position', [100, 100, 400, 400]);
hold on;

b = bar(outdata, 'grouped'); 

[numGroups, numBars] = size(outdata);

groupColors = [1 0 0; 0 0 1; 0.5 0.5 0.5];

for i = 1:numGroups
    for j = 1:numBars
        b(j).FaceColor = 'flat';
        b(j).CData(i, :) = groupColors(i, :);
    end
end

xticks([1, 2, 3]);
xticklabels({'PM', 'NM', 'NSM'});
xlim([0.5, 3.5]);
ylim([0, 100]);
xlabel('Granule cell responses');
ylabel('Fraction (%)');

groupCenters = b(1).XData;
barWidth = b(1).BarWidth / numBars;

thresholdLabels = {'0', '1', '2', '3'};
yTextPos = -3; 

for i = 1:length(groupCenters)
    for j = 1:numBars
        xPos = groupCenters(i) - 1.5 * barWidth + (j - 1) * barWidth;
        text(xPos, yTextPos, thresholdLabels{j}, ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
end

ax = gca;
ax.XAxis.TickLength = [0 0];
set(gca, 'LineWidth', 1, 'FontSize', 12, 'TickDir', 'out', 'TickLength', [.025, .025]);

ylim([-6, 100]);
yticks([0 25 50 75 100])


mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
folderPath = fullfile(mfbFolderPath, 'Figures', 'Figure8', currentDate);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
fileName = ['GC_integration_output.png'];
fullFilePath = fullfile(folderPath, fileName);
% save
print(fullFilePath, '-dpng', '-r300');


%%
figure('Position',[100,100,350,300])
subplot(111)
hold on

d_max = 250;
dx = 5;
bins = 0:dx:d_max;

N_tot = N_GC;
[y,x] = histcounts(dd_PM, bins, 'Normalization','cumcount');
dx = diff(x(1:2));
xx = x(1:end-1) + dx/2;
%N_tot = length(dd_PM);
plot(xx,y/N_tot, 'LineWidth',2, 'color','r')

[y,x] = histcounts(dd_NM, bins, 'Normalization','cumcount');
dx = diff(x(1:2));
xx = x(1:end-1) + dx/2;
%N_tot = length(dd_NM);
plot(xx,y/N_tot, 'LineWidth',2, 'color','b')

[y,x] = histcounts(dd_NS, bins, 'Normalization','cumcount');
dx = diff(x(1:2));
xx = x(1:end-1) + dx/2;
%N_tot = length(dd_NS);
plot(xx,y/N_tot, 'LineWidth',2,'color', 'k')

xlim([0,250])

legend('PM', 'NM', 'NSM', 'Location', 'best')
legend boxoff

xlabel('Distance (μm)')
ylabel('MFBs (avg. per GC)')
set(gca, 'LineWidth', 1, 'FontSize', 15, 'TickDir', 'out')


mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
folderPath = fullfile(mfbFolderPath, 'Figures', 'Figure8', currentDate);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

fileName = ['GC_dist_dist.png'];
fullFilePath = fullfile(folderPath, fileName);
% save
print(fullFilePath, '-dpng', '-r300');

%%
figure('Position',[100,100,350,300])
subplot(111)
hold on

prob_PM = exp(-dd_PM/dend_L);
prob_NM = exp(-dd_NM/dend_L);
prob_NS = exp(-dd_NS/dend_L);

%d_max = 250;
dx = 5;
x = 0:dx:d_max;
pr_PM = nan*zeros(1,length(x)-1);
pr_NM = nan*zeros(1,length(x)-1);
pr_NS = nan*zeros(1,length(x)-1);

xx = x(1:end-1) + dx/2;
p_tot = nansum(prob_PM(:)) + nansum(prob_NM(:)) + nansum(prob_NS(:));
for igc = 1:length(x)-1
    ids = logical((dd_PM > x(igc)) .* (dd_PM < x(igc+1)));
    pr_PM(igc) = nansum(prob_PM(ids))/p_tot*100;
    
    ids = logical((dd_NM > x(igc)) .* (dd_NM < x(igc+1)));
    pr_NM(igc) = nansum(prob_NM(ids))/p_tot*100;

    ids = logical((dd_NS > x(igc)) .* (dd_NS < x(igc+1)));
    pr_NS(igc) = nansum(prob_NS(ids))/p_tot*100;
end

plot(xx, pr_PM, 'r-', 'LineWidth',2)
plot(xx, pr_NM, 'b-', 'LineWidth',2)
plot(xx, pr_NS, 'k-', 'LineWidth', 2);


xlim([0,250])

legend('PM', 'NM', 'NSM', 'Location', 'best')
legend boxoff

xlabel('Distance (μm)')
ylabel('Conn. prob. (%)')
set(gca, 'LineWidth', 1, 'FontSize', 15, 'TickDir', 'out')


fileName = ['GC_dist_prob.png'];
fullFilePath = fullfile(folderPath, fileName);
% save
print(fullFilePath, '-dpng', '-r300');



