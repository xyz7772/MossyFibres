clear all; close all;clc
initialize
%%
folder_names = {'Animal1', 'Animal2', 'Animal3', 'Animal4'};
concat_folder = 'X:\MFB\Concatenation data'; % path to concatenated data folder
confiles = dir(fullfile(concat_folder, 'concat_*.mat'));
N_GC = 300;
dend_L = 21.7;
dt = 10;
ite_l = 10;

for th = 0:3  % th = 0,1,2,3
    prob_PM = []; dd_PM = [];
    prob_NM = []; dd_NM = [];
    prob_NS = []; dd_NS = [];

    results_itr = struct('dist', [], 'num_PM', [], 'num_NM', [], 'num_NSM', [], ...
        'pm_out', [], 'nm_out', [], 'nsm_out', [], 'PPPP', [], 'PPPS', [], ...
        'PPPN', [], 'PPNN', [], 'PPNS', [], 'PPSS', [], 'PNNN', [], 'PNNS', [], ...
        'PNSS', [], 'PSSS', [], 'NNNN', [], 'NNNS', [], 'NNSS', [], 'NSSS', [], ...
        'SSSS', []);

    for fi = 1:length(confiles)
        filePath = fullfile(data_home, confiles(fi).name);
        load(filePath);

        if strcmp (file, 'Animal1') ==1
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
            unimf = [];
            pm_out = 0;
            nm_out = 0;
            nsm_out = 0;
            dist_all = zeros(1, N_GC);

            PPPP_count = 0; PPPN_count = 0; PPPS_count = 0; PPNN_count = 0;
            PPNS_count = 0; PPSS_count = 0; PNNN_count = 0; PNNS_count = 0;
            PNSS_count = 0; PSSS_count = 0; NNNN_count = 0; NNNS_count = 0;
            NNSS_count = 0; NSSS_count = 0; SSSS_count = 0;

            Xall = xyz(:, 1);
            Yall = xyz(:, 2);
            Zall = xyz(:, 3);

            x_GC = my_rand([1, N_GC], min(Xall(:)) + dend_L, max(Xall(:)) - dend_L);
            y_GC = my_rand([1, N_GC], min(Yall(:)) + dend_L, max(Yall(:)) - dend_L);
            z_GC = my_rand([1, N_GC], min(Zall(:)) + dend_L, max(Zall(:)) - dend_L);

            xyz_GC = [x_GC; y_GC; z_GC]';

            for igc = 1:N_GC
                d = xyz_GC(igc,:) - xyz(:,:);
                dd_PM = [dd_PM, sqrt(sum(d(pmids,:).^2,2))'];
                dd_NM = [dd_NM, sqrt(sum(d(nmids,:).^2,2))'];
                dd_NS = [dd_NS, sqrt(sum(d(nsmids,:).^2,2))'];
                dd2 = sqrt(sum(d.^2, 2));
                [a, b] = sort(dd2);

                idx = 1:4;
                % check independence
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
            end

            results_itr(fi).dist(:, ite) = mean(mean_dis);
            results_itr(fi).num_PM(ite) = sum(out_mod_all == 1);
            results_itr(fi).num_NM(ite) = sum(out_mod_all == -1);
            results_itr(fi).num_NSM(ite) = sum(out_mod_all == 0);
            results_itr(fi).pm_out(ite) = pm_out;
            results_itr(fi).nm_out(ite) = nm_out;
            results_itr(fi).nsm_out(ite) = nsm_out;
        end
    end

    filePath = fullfile([path_home '\preprocessing\GC modeling\results_itr_th' num2str(th) '.mat']);
    save(filePath, 'results_itr');
end

%%
all_PPPP = []; all_PPPN = []; all_PPPS = [];
all_PPNN = []; all_PPNS = []; all_PPSS = [];
all_PNNN = []; all_PNNS = []; all_PNSS = [];
all_PSSS = []; all_NNNN = []; all_NNNS = [];
all_NNSS = []; all_NSSS = []; all_SSSS = [];

for i = 1:4
    load([path_home '\preprocessing\GC modeling\results_itr_th' num2str(i-1) '.mat']);

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

patterns = {'NNNN', 'NNNS', 'PNNN', 'NNSS', 'PNNS', 'NSSS', 'PNSS', 'PPNN', 'SSSS', 'PSSS', 'PPNS', 'PPSS', 'PPPN', ...
    'PPPS', 'PPPP'};
mean_values_map = containers.Map(patterns, [mean_NNNN, mean_NNNS, mean_PNNN, mean_NNSS, mean_PNNS, ...
    mean_NSSS, mean_PNSS, mean_PPNN, mean_SSSS, mean_PSSS, mean_PPNS, mean_PPSS, mean_PPPN, mean_PPPS, mean_PPPP]);
sem_values_map = containers.Map(patterns, [sem_NNNN, sem_NNNS, sem_PNNN, sem_NNSS, sem_PNNS, ...
    sem_NSSS, sem_PNSS, sem_PPNN, sem_SSSS, sem_PSSS, sem_PPNS, sem_PPSS, sem_PPPN, sem_PPPS, sem_PPPP]);

mean_values = cellfun(@(x) mean_values_map(x), patterns);
sem_values = cellfun(@(x) sem_values_map(x), patterns);

theoretical_values =[1.1,3.78,4.03,4.93,10.7,2.92,10,7.63,0.676,3.34,15.2,8.32,9.69,11.3,6.36]; % pre-simulated value

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
for i = 1:length(patterns)
    data_points = eval(['all_' patterns{i}]) * 100;
    num_points = length(data_points);
    
    if num_points > 0
        scatter(repmat(x_pos(i), num_points, 1), data_points, '.','MarkerEdgeColor',[0.5,0.5,0.5]);
    end
end

for i = 1:length(theoretical_values)
    y_value = theoretical_values(i);
    x_pos = i;
    scatter(x_pos, y_value, 30, [255, 178, 102]/255, 'filled', 'd');
end

ylabel('Fraction (%)');
yticks([0 15 30]);
ylim([0 30]);
xticklabels(patterns);
set(gca, 'LineWidth', 1, 'FontSize', 15, 'TickDir', 'out', 'Box', 'off');

currentDate = datestr(now, 'yyyy-mm-dd');
savepath2 = fullfile(savepath, 'Figures', 'Figure8', currentDate);
if ~exist(savepath2, 'dir')
    mkdir(savepath2);
end

fileName = ['GC_INPUT'];
fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%%
colors = {'r', 'b', 'k'};
outdata = zeros(3, 4);
semdata = zeros(3, 4);

all_PM_out_p = cell(1, 4);
all_NM_out_p = cell(1, 4);
all_NSM_out_p = cell(1, 4);

for i = 0:3
    file = ['results_itr_th' num2str(i) '.mat'];
    load(file);
    
    all_num_PM  = cat(1, results_itr(:).num_PM);
    all_num_NM  = cat(1, results_itr(:).num_NM);
    all_num_NSM  = cat(1, results_itr(:).num_NSM);

    PM_out_p = all_num_PM / mean(all_num_PM + all_num_NM + all_num_NSM) * 100;
    NM_out_p = all_num_NM / mean(all_num_PM + all_num_NM + all_num_NSM) * 100;
    NSM_out_p = all_num_NSM / mean(all_num_PM + all_num_NM + all_num_NSM) * 100;

    outdata(1, i+1) = mean(PM_out_p);
    outdata(2, i+1) = mean(NM_out_p);
    outdata(3, i+1) = mean(NSM_out_p);

    semdata(1, i+1) = std(PM_out_p) / sqrt(length(PM_out_p));
    semdata(2, i+1) = std(NM_out_p) / sqrt(length(NM_out_p));
    semdata(3, i+1) = std(NSM_out_p) / sqrt(length(NSM_out_p));

    all_PM_out_p{i+1} = PM_out_p;
    all_NM_out_p{i+1} = NM_out_p;
    all_NSM_out_p{i+1} = NSM_out_p;
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

ax = gca;
ax.XAxis.TickLength = [0 0];
set(gca, 'LineWidth', 1, 'FontSize', 12, 'TickDir', 'out', 'TickLength', [.025, .025]);

ylim([-6, 100]);
yticks([0 25 50 75 100]);

bar_x_positions = nan(numGroups, numBars);
for j = 1:numBars
    bar_x_positions(:, j) = b(j).XEndPoints;
end
for j = 1:numBars
    errorbar(bar_x_positions(:, j), outdata(:, j), semdata(:, j), 'k.', 'LineWidth', 1.5);
end
for j = 1:numBars  
    scatter(bar_x_positions(1, j) * ones(size(all_PM_out_p{j})), all_PM_out_p{j}, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.7);
    scatter(bar_x_positions(2, j) * ones(size(all_NM_out_p{j})), all_NM_out_p{j}, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.7);
    scatter(bar_x_positions(3, j) * ones(size(all_NSM_out_p{j})), all_NSM_out_p{j}, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.7);
end

% Save
fileName = ['GC_integration_output'];
fullFilePathPDF = fullfile(savepath2, [fileName, '.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

