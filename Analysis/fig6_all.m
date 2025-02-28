%% Select
folder_names = {
   '171212_16_19_37';
   '191018_13_39_41';
   %'191018_13_56_55'; % only run at the beginning
   '191018_14_30_00'; % misaligned?
   %'191018_14_11_33'; avg speed =0
   '191209_13_44_12';
   %'191209_14_04_14'; %beh not good
   '191209_14_32_39';
   '191209_15_01_22'; %L_state is all QW, but good wshiker  
   '191209_14_18_13';
   '191209_14_46_58';
   '200130_13_21_13 FunctAcq';
   '200130_13_36_14 FunctAcq';
   '200130_13_49_09 FunctAcq';
   %'200130_14_02_12 FunctAcq';
   '200130_14_15_24 FunctAcq';
   '200130_14_29_30 FunctAcq';
    };
%% get beta for each MF

plot_decode_all_MF = 1;

if plot_decode_all_MF == 1

    beta_loco_all_files = {};
    beta_wsk_all_files = {};
    beta_state_all_files = {};
    beta_loco_sorted_indices_all = {}; 
    beta_wsk_sorted_indices_all = {}; 
    beta_state_sorted_indices_all = {};

    cvev_loco_all_files = {};
    cvev_wsk_all_files = {};
    cvev_state_all_files = {};

    for file_i = 1:length(folder_names)
        file = char(folder_names(file_i));
        quickAnalysis;
        Y_loco = MI_wheel_r';
        Y_wsk = MI_whisker_r';
        Y_state = L_state';
        X = dff_rz';

        Y_loco = Y_loco(valid_t);
        Y_wsk = Y_wsk(valid_t);
        Y_state = Y_state(valid_t);
        X = X(valid_t, :);
        
        Xtrain = X(:, :) - mean(X(:, :));
        Xtrain2 = zscore(X);
        [Ytrain_loco_z, ~, ~] = zscore(Y_loco);
        [Ytrain_wsk_z, ~, ~] = zscore(Y_wsk); 
        Ytrain_state = Y_state;

        lambda = 0.1;
        beta_loco_i = ridge(Ytrain_loco_z, Xtrain, lambda, 0);
        beta_wsk = ridge(Ytrain_wsk_z, Xtrain, lambda, 0);
        beta_state = ridge(Ytrain_state, Xtrain, lambda, 0);

        [~, sorted_indices_loco] = sort(beta_loco_i(2:end), 'descend');
        [~, sorted_indices_wsk] = sort(beta_wsk(2:end), 'descend');
        [~, sorted_indices_state] = sort(beta_state(2:end), 'descend');

        beta_loco_all_files{file_i} = beta_loco_i;
        beta_wsk_all_files{file_i} = beta_wsk;
        beta_state_all_files{file_i} = beta_state;
        beta_loco_sorted_indices_all{file_i} = sorted_indices_loco;
        beta_wsk_sorted_indices_all{file_i} = sorted_indices_wsk;
        beta_state_sorted_indices_all{file_i} = sorted_indices_state;

    end

    folderPath = 'X:\MFB\MFB_AH_2023\Correlation_data\decode';
    fileName = 'Sorted_Beta.mat';
    fullPath = fullfile(folderPath, fileName);
    save(fullPath, 'beta_loco_all_files', 'beta_wsk_all_files', 'beta_state_all_files','beta_loco_sorted_indices_all', 'beta_wsk_sorted_indices_all' ,'beta_state_sorted_indices_all');
end


%% Find optimal MFs by lasso
plot_decode_all_MF=1;

for file_i = 1:length(folder_names)
    
    file=char(folder_names(file_i));
    quickAnalysis;

    Y = MI_wheel_r';
    Y2 = MI_whisker_r';
    Y3 = L_state';
    X = dff_rz';

    Y= Y(valid_t);
    Y2= Y2(valid_t);
    Y3= Y3(valid_t);
    X= X(valid_t, :);
    L_state = L_state(valid_t);

    % may need to write a loop, from best MF to top 100, 
    X_n_loco = X;
    X_n_wsk = X;
    X_n_state = X; 

    lambda = 1e-3:1e-4:0.05;
    [B_loco, FitInfo_loco] = lasso(X_n_loco, Y, 'Lambda',lambda,'CV', 10);
    c_Min = B_loco(:, FitInfo_loco.Index1SE);
    if sum(c_Min ~= 0) >0
        num_MF_loco(file_i) = sum(c_Min ~= 0);
    else
        num_MF_loco(file_i) = nan;
    end

    [B_wsk, FitInfo_wsk] = lasso(X_n_wsk, Y2, 'Lambda',lambda,'CV', 10);
    c_Min = B_wsk(:, FitInfo_wsk.Index1SE);
    if sum(c_Min ~= 0) >0
        num_MF_wsk(file_i) = sum(c_Min ~= 0);
    else
        num_MF_wsk(file_i) = nan;
    end

    [B_state, FitInfo_state] = lasso(X_n_state, Y3, 'Lambda',lambda,'CV', 10); 
    c_Min = B_state(:, FitInfo_state.Index1SE);
    if sum(c_Min ~= 0) >0
        num_MF_state(file_i) = sum(c_Min ~= 0);
    else
        num_MF_state(file_i) = nan;
    end

end

save('X:\MFB\MFB_AH_2023\Correlation_data\decode\Best_MFs_num.mat', "num_MF_state","num_MF_loco","num_MF_wsk");


%% plot traces
k = 4;
set(figure, 'Position', [400, 100, 500, 290]);
% Plot Loco Actual
cell1 = Actual_loco_all{k}{end};
minL = min(cellfun(@(x) numel(x), cell1));
Data = zeros(minL, numel(cell1));
for i = 1:numel(cell1)
    Data(:, i) = cell1{i}(1:minL);
end
average_test = mean(Data, 2);

subplot(6, 1, 1);
h1 = plot(average_test, 'color',[0.6,0.6,0.6],'LineWidth',1);
box('off')
axis off;

disp(Cvev_loco{k}(1:5))  

for plot_i=1:5  
cell2 = Decode_loco_all{k}{plot_i};
minL = min(cellfun(@(x) numel(x), cell2));
Data = zeros(minL, numel(cell2));
for i = 1:numel(cell2)
    Data(:, i) = cell2{i}(1:minL);
end
    average_pred = mean(Data, 2);
    
    subplot(6, 1, plot_i+1);
    h2=plot(average_pred,'color',[173,210,157]/255,'LineWidth',1); % Predicted in blue
    box('off')
    axis off;
end


xlabel('Time');
ylabel('Corr')
ax = gca;
xlims = ax.XLim;
ylims = ax.YLim;
x0 = xlims(2) - 1000; % 10s
y0 = ylims(1); %

hold on
plot(ax, [x0, x0+1000], [y0, y0], 'k-', 'LineWidth', 3);
text(ax, x0, y0 - diff(ylims) * 0.2, '10 s', 'FontSize', 13, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
folderPath = fullfile(mfbFolderPath, 'Figures', 'Figure5', currentDate);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
fileName = [ 'Ridge_Loco_' folder_names{k} '_MF#'];
fullFilePathPDF = fullfile(folderPath, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


% Plot Wsk
set(figure, 'Position', [400, 100, 500, 290]);

% Plot whisk Actual
cell1=Actual_whisk_all{k}{end};
minL = min(cellfun(@(x) numel(x), cell1));
Data = zeros(minL, numel(cell1));
for i = 1:numel(cell1)
    Data(:, i) = cell1{i}(1:minL);
end
average_test = mean(Data, 2);

subplot(6, 1, 1);
h1=plot(average_test, 'color',[0.6,0.6,0.6],'LineWidth',1);
box('off')
axis off;

disp(Cvev_wsk{k}(1:5))  
for plot_i=1:5
cell2=Decode_whisk_all{k}{plot_i};
minL = min(cellfun(@(x) numel(x), cell2));
Data = zeros(minL, numel(cell2));
for i = 1:numel(cell2)
    Data(:, i) = cell2{i}(1:minL);
end
average_pred = mean(Data, 2);
subplot(6, 1, plot_i+1);
h2=plot(average_pred,'color',[255,163,26]/255,'LineWidth',1);
box('off')
axis off;
end

xlabel('Time');
ylabel('Corr')
ax = gca;
xlims = ax.XLim;
ylims = ax.YLim;
x0 = xlims(2) - 1000; % 10s
y0 = ylims(1);

hold on
plot(ax, [x0, x0+1000], [y0, y0], 'k-', 'LineWidth', 3);
text(ax, x0, y0 - diff(ylims) * 0.2, '10 s', 'FontSize', 13, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

fileName = [ 'Ridge_whisk_' folder_names{k} '_MF#'];
fullFilePathPDF = fullfile(folderPath, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

% Plot States
set(figure, 'Position', [400, 100, 500, 290]);

% Plot state Actual
cell1=Actual_state_all{k}{end};
minL = min(cellfun(@(x) numel(x), cell1));
Data = zeros(minL, numel(cell1));
for i = 1:numel(cell1)
    Data(:, i) = cell1{i}(1:minL);
end
average_test = cell1{i};

subplot(6, 1, 1);
h1=plot(average_test, 'color',[0.6,0.6,0.6],'LineWidth',1);
box('off')
axis off;
disp(Cvev_state{k}(1:5))
for plot_i=1:5
cell2=Decode_state_all{k}{plot_i};

minL = min(cellfun(@(x) numel(x), cell2));
Data = zeros(minL, numel(cell2));
for i = 1:numel(cell2)
    Data(:, i) = cell2{i}(1:minL);
end
average_pred = mean(Data, 2);

subplot(6, 1, plot_i+1);
h2=plot(average_pred,'color',[255,53,255]/255,'LineWidth',1);
box('off')
axis off;
end

xlabel('Time');
ylabel('Corr')
ax = gca;
xlims = ax.XLim;
ylims = ax.YLim;
x0 = xlims(2) - 1000; % 10s
y0 = ylims(1);

hold on
plot(ax, [x0, x0+1000], [y0, y0], 'k-', 'LineWidth', 3);
text(ax, x0, y0 - diff(ylims) * 0.2, '10 s', 'FontSize', 13, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

fileName = [ 'Ridge_state_' folder_names{k} '_MF#'];
fullFilePathPDF = fullfile(folderPath, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%% fig. 6b

load('X:\MFB\MFB_AH_2023\Correlation_data\decode\results_all_MF_sorted_beta.mat')
figure('Position', [400,100,450,350]);
for i=1:length(Cvev_loco)
    if i==6
        continue
    end
    cell_data = Cvev_loco{i};

    modes = 1:length(cell_data);

    hold on;
    plot((modes-1)*15+1, cell_data, '-', 'Color', [51, 204, 204]/255, 'LineWidth', 1.5)

end

allData = [];
for i=1:length(Cvev_loco)
    cell_data = Cvev_loco{i};
    allData = [allData; cell_data];
end

overallMean = nanmean(allData);
overallSEM = nanstd(allData) / sqrt(length(allData));
hold on;
x = (modes-1)*15+1;
fill([x, fliplr(x)], [overallMean + overallSEM, fliplr(overallMean - overallSEM)], ...
     [0.6, 0.6, 0.6], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
plot(x, overallMean, 'k', 'LineWidth', 1.5);

title('Loco');
xlabel('MFA#');
ylabel('CVEV');

xlim([1 240]);
xticks([1 240])
ylim([0 0.85]);
yticks([0 0.85])

set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out')

fileName = ['Decode_allMF#_loco'];
fullFilePathPDF = fullfile(folderPath, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


figure('Position', [400, 100, 450, 350]);
for i=1:length(Cvev_wsk)
    cell_data = Cvev_wsk{i};
    modes = 1:length(cell_data);
    hold on;
    plot((modes-1)*15+1, cell_data, '-', 'Color', [51, 204, 204]/255, 'LineWidth', 1.5)
end


allData = [];
for i=1:length(Cvev_wsk)
    cell_data = Cvev_wsk{i};
    allData = [allData; cell_data];
end

overallMean = nanmean(allData);
overallSEM = nanstd(allData) / sqrt(length(allData));
hold on;
x = (modes-1)*15+1;
fill([x, fliplr(x)], [overallMean + overallSEM, fliplr(overallMean - overallSEM)], ...
     [0.6, 0.6, 0.6], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
plot(x, overallMean, 'k', 'LineWidth', 1.5);

title('Whisk');
xlabel('MFA#');
ylabel('CVEV');

xlim([1 240]);
xticks([1 240])
ylim([0 0.85]);
yticks([0 0.85])

drawnow;
set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out')

fileName = ['Decode_allMF#_whisk'];
fullFilePathPDF = fullfile(folderPath, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


figure('Position', [400, 100, 450, 350]);
for i=1:length(Cvev_loco)
    cell_data = Cvev_state{i};
    modes = 1:length(cell_data);
    hold on;
    plot((modes-1)*15+1, cell_data, '-', 'Color', [51, 204, 204]/255, 'LineWidth', 1.5)
end

allData = [];
for i=1:length(Cvev_state)
    cell_data = Cvev_state{i};
    allData = [allData; cell_data];
end

overallMean = nanmean(allData);
overallSEM = nanstd(allData) / sqrt(length(allData));
hold on;
x = (modes-1)*15+1;
fill([x, fliplr(x)], [overallMean + overallSEM, fliplr(overallMean - overallSEM)], ...
     [0.6, 0.6, 0.6], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
plot(x, overallMean, 'k', 'LineWidth', 1.5);

title('States');
xlabel('MFA#');
ylabel('CVEV');

xlim([1 240]);
xticks([1 240])
ylim([0 1]);
yticks([0 1])

drawnow;

set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out')

fileName = ['Decode_allMF#_state'];
fullFilePathPDF = fullfile(folderPath, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');


%% plot by first/optimal MF based on corr with behaviour

load ('X:\MFB\MFB_AH_2023\Correlation_data\decode\Best_MFs_num.mat')
load('X:\MFB\MFB_AH_2023\Correlation_data\decode\Sorted_Beta.mat')

Cvev_loco=[];
Cvev_wsk=[];
Cvev_state=[];
itr=1;

for file_i = 1:length(folder_names)
    file = char(folder_names(file_i));
    quickAnalysis;

    Y = MI_wheel_r';
    Y2 = MI_whisker_r';
    Y3 = L_state';
    X = dff_rz';


    Y = Y(valid_t);
    Y2 = Y2(valid_t);
    Y3 = Y3(valid_t);
    X = X(valid_t, :);
    L_state = L_state(valid_t);

    cvev_loco_results = zeros(itr,1);
    cvev_wsk_results = zeros(itr,1);
    cvev_state_results = zeros(itr,1);

    for j = 1:itr
       % 1. Best MFs in loco (Y)
        N = num_MF_loco(file_i);
        N=1
        selected_loco = randperm(size(X,2), N);
        
        beta_loco_i = beta_loco_sorted_indices_all{file_i};
        beta_loco_arr= zeros(size(X, 2), 1);
        beta_loco_arr(beta_loco_i) = 1:length(beta_loco_i);
        
        [~, sorted_indices_corr] = sort(corr(X,Y), 'descend');
        corr_loco_rank = zeros(size(X, 2), 1);
        corr_loco_rank(sorted_indices_corr) = 1:length(sorted_indices_corr);
        
        beta_norm = (beta_loco_arr - min(beta_loco_arr)) / (max(beta_loco_arr) - min(beta_loco_arr));
        corr_norm = (corr_loco_rank - min(corr_loco_rank)) / (max(corr_loco_rank) - min(corr_loco_rank));
        
        combined_rank_loco = 0.5 * beta_norm + 0.5 * corr_norm;
        [~, final_rank_loco] = sort(combined_rank_loco, 'ascend');
        selected_loco_sorted = final_rank_loco(1:N); % 排序选取
        
        % 2. Best MFs in wsk (Y2)
        N = num_MF_wsk(file_i);
        N =1
        selected_wsk = randperm(size(X,2), N);
        
        beta_wsk = beta_wsk_sorted_indices_all{file_i};
        beta_wsk_array = zeros(size(X, 2), 1);
        beta_wsk_array(beta_wsk) = 1:length(beta_wsk);
        
        [~, sorted_indices_corr_wsk] = sort(corr(X,Y2), 'descend');
        corr_wsk_rank = zeros(size(X, 2), 1);
        corr_wsk_rank(sorted_indices_corr_wsk) = 1:length(sorted_indices_corr_wsk);
        
        % 归一化排名
        beta_wsk_norm = (beta_wsk_array - min(beta_wsk_array)) / (max(beta_wsk_array) - min(beta_wsk_array));
        corr_wsk_norm = (corr_wsk_rank - min(corr_wsk_rank)) / (max(corr_wsk_rank) - min(corr_wsk_rank));
        
        % 计算综合排名
        combined_rank_wsk = 0.5 * beta_wsk_norm + 0.5 * corr_wsk_norm;
        [~, final_rank_wsk] = sort(combined_rank_wsk, 'ascend');
        selected_wsk_sorted = final_rank_wsk(1:N); % 排序选取
        
        % 3. Best MFs in states (Y3)
        N = num_MF_state(file_i);
        N=1
        if isnan(N)
            nancheck = 1;
        else
            nancheck = 0;
            selected_state = randperm(size(X,2), N);
        end
        
        beta_state = beta_state_sorted_indices_all{file_i};
        beta_state_array = zeros(size(X, 2), 1);
        beta_state_array(beta_state) = 1:length(beta_state);
        
        [~, sorted_indices_corr_state] = sort(corr(X,Y3), 'descend');
        corr_state_rank = zeros(size(X, 2), 1);
        corr_state_rank(sorted_indices_corr_state) = 1:length(sorted_indices_corr_state);
        
        % 归一化排名
        beta_state_norm = (beta_state_array - min(beta_state_array)) / (max(beta_state_array) - min(beta_state_array));
        corr_state_norm = (corr_state_rank - min(corr_state_rank)) / (max(corr_state_rank) - min(corr_state_rank));
        
        % 计算综合排名
        combined_rank_state = 0.5 * beta_state_norm + 0.5 * corr_state_norm;
        [~, final_rank_state] = sort(combined_rank_state, 'ascend');
        
        % 选择最优神经元
        if nancheck == 0
            selected_state_sorted = final_rank_state(1:N); 
        else 
            selected_state_sorted = 1; 
        end

        X_n_loco = X(:, selected_loco_sorted);
        X_n_wsk = X(:, selected_wsk_sorted);
        X_n_state = X(:, selected_state_sorted); 

        cvp = cvpartition(size(X, 1), 'KFold', 5);
        cvev_loco = zeros(cvp.NumTestSets, 1);
        cvev_wsk = zeros(cvp.NumTestSets, 1);
        cvev_state = zeros(cvp.NumTestSets, 1);

        [B_loco_all, FitInfo_loco_all] = lasso(X_n_loco, Y, 'CV', 10);
        lambda_opt_loco = FitInfo_loco_all.LambdaMinMSE;
        intercept_loco = FitInfo_loco_all.Intercept(FitInfo_loco_all.IndexMinMSE);
        
        [B_wsk_all, FitInfo_wsk_all] = lasso(X_n_wsk, Y2, 'CV', 10);
        lambda_opt_wsk = FitInfo_wsk_all.LambdaMinMSE;
        intercept_wsk = FitInfo_wsk_all.Intercept(FitInfo_wsk_all.IndexMinMSE);
        
        [B_state_all, FitInfo_state_all] = lasso(X_n_state, Y3, 'CV', 10);
        lambda_opt_state = FitInfo_state_all.LambdaMinMSE;
        intercept_state = FitInfo_state_all.Intercept(FitInfo_state_all.IndexMinMSE);
        
        for i = 1:cvp.NumTestSets
            trainIdx = cvp.training(i);
            testIdx = cvp.test(i);
        
            B_loco = lasso(X_n_loco(trainIdx, :), Y(trainIdx), 'Lambda', lambda_opt_loco);
            B_wsk = lasso(X_n_wsk(trainIdx, :), Y2(trainIdx), 'Lambda', lambda_opt_wsk);
            B_state = lasso(X_n_state(trainIdx, :), Y3(trainIdx), 'Lambda', lambda_opt_state);
        
            Y_pred_loco = X_n_loco(testIdx, :) * B_loco + intercept_loco;
            Y_pred_wsk = X_n_wsk(testIdx, :) * B_wsk + intercept_wsk;
            Y_pred_state = X_n_state(testIdx, :) * B_state + intercept_state;
        
            cvev_loco(i) = cal_cvev(Y(testIdx), Y_pred_loco);
            cvev_wsk(i) = cal_cvev(Y2(testIdx), Y_pred_wsk);
            cvev_state(i) = cal_cvev(Y3(testIdx), Y_pred_state);
        end
    
            % Store results for each iteration
            cvev_loco_results(j) = mean(cvev_loco);
            cvev_wsk_results(j) = mean(cvev_wsk);
            if nancheck ~= 1
                cvev_state_results(j) = mean(cvev_state);
            end
        end

    % Average CVEV across all iterations
    meanCVEV_loco = mean(cvev_loco_results);
    meanCVEV_wsk = mean(cvev_wsk_results);
    if nancheck ~= 1
        meanCVEV_state = mean(cvev_state_results);
    else
        meanCVEV_state = nan;
    end

    % Store or output the results
    Cvev_loco{file_i} = meanCVEV_loco
    Cvev_wsk{file_i} = meanCVEV_wsk
    Cvev_state{file_i} = meanCVEV_state
end


%% Save the results
save('X:\MFB\MFB_AH_2023\Correlation_data\decode\results_lasso_1_MF_sorted_beta.mat','Cvev_loco','Cvev_wsk','Cvev_state');

%% plot for fig 6c

load('X:\MFB\MFB_AH_2023\Correlation_data\decode\results_lasso_MFs_sorted_beta.mat');
result_wsk_lasso = cellfun(@(x) x(1), Cvev_wsk);
result_loco_lasso = cellfun(@(x) x(1), Cvev_loco);
result_state_lasso = cellfun(@(x) x(1), Cvev_state);

load('X:\MFB\MFB_AH_2023\Correlation_data\decode\results_lasso_1_MF_sorted_beta.mat');
result_wsk_lasso_1 = cellfun(@(x) x(1), Cvev_wsk);
result_loco_lasso_1 = cellfun(@(x) x(1), Cvev_loco);
result_state_lasso_1 = cellfun(@(x) x(1), Cvev_state);

% WSK
set(figure, 'Position', [400, 100, 350, 350]);
hold on;

bar_data_wsk = [mean(result_wsk_lasso), mean(result_wsk_lasso_1)];
b1 = bar([1, 2], bar_data_wsk, 'FaceColor', 'k', 'FaceAlpha', 0.1);
scatter(ones(size(result_wsk_lasso)), result_wsk_lasso, 50, 'k');
scatter(ones(size(result_wsk_lasso_1)) + 1, result_wsk_lasso_1, 50, 'k');
for i = 1:length(result_wsk_lasso)
    plot([1, 2], [result_wsk_lasso(i), result_wsk_lasso_1(i)], 'k--', 'Color', [0.5 0.5 0.5]);
end

set(gca, 'XTick', [1, 2], 'XTickLabel', {'Lasso', 'Best MF'},'FontSize',15,'LineWidth',1);
ylabel('CVEV');
title('Whisk');
yticks([0 0.5 1])
ylim([0 1])
[p_wsk, h_wsk] = signrank(result_wsk_lasso, result_wsk_lasso_1)


fileName = ['Lasso_MF1vsall_whisk'];
fullFilePathPDF = fullfile(folderPath, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

%Loco
set(figure, 'Position', [400, 100, 350, 350]);
hold on;

bar_data_loco = [mean(result_loco_lasso), mean(result_loco_lasso_1)];
b2 = bar([1, 2], bar_data_loco, 'FaceColor', 'k', 'FaceAlpha', 0.1);
scatter(ones(size(result_loco_lasso)), result_loco_lasso, 50, 'k');
scatter(ones(size(result_loco_lasso_1)) + 1, result_loco_lasso_1, 50, 'k');

for i = 1:length(result_loco_lasso)
    plot([1, 2], [result_loco_lasso(i), result_loco_lasso_1(i)], 'k--', 'Color', [0.5 0.5 0.5]);
end

set(gca, 'XTick', [1, 2], 'XTickLabel', {'Lasso', 'Best MF'},'FontSize',15,'LineWidth',1);
ylabel('CVEV');
title('Loco');
yticks([0 0.5 1])
ylim([0 1])


mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
folderPath = fullfile(mfbFolderPath, 'Figures', 'Figure6', currentDate);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

fileName = ['Lasso_MF1vsall_loco'];
fullFilePathPDF = fullfile(folderPath, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

[p_loco, h_loco] = signrank(result_loco_lasso, result_loco_lasso_1)


% state
set(figure, 'Position', [400, 100, 350, 350]);
hold on

bar_data_state = [nanmean(result_state_lasso), nanmean(result_state_lasso_1)]; 
b3 = bar([1, 2], bar_data_state, 'FaceColor', 'k', 'FaceAlpha', 0.1); 
scatter(ones(size(result_state_lasso)), result_state_lasso, 50, 'k');
scatter(ones(size(result_state_lasso_1)) + 1, result_state_lasso_1, 50, 'k'); 
for i = 1:length(result_state_lasso)
    plot([1, 2], [result_state_lasso(i), result_state_lasso_1(i)], 'k--', 'Color', [0.5 0.5 0.5]);
end

set(gca, 'XTick', [1, 2], 'XTickLabel', {'Lasso', 'Best MF'},'FontSize',15,'LineWidth',1);
ylabel('CVEV');
title('State');
yticks([0 0.5 1])
ylim([0 1])
[p_state, h_state] = signrank(result_state_lasso, result_state_lasso_1)

fileName = ['Lasso_MF1vsall_state'];
fullFilePathPDF = fullfile(folderPath, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');

