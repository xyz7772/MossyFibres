clear all; close all; clc
%% all folders
folder_names = {
   '171212_16_19_37';
   '191018_13_39_41';
   %'191018_13_56_55'; % sessions are removed due to no beh
   '191018_14_30_00';
   %'191018_14_11_33';
   '191209_13_44_12';
   %'191209_14_04_14';
   '191209_14_32_39';
   '191209_15_01_22'; % No loco,state, but good whisker  
   '191209_14_18_13';
   '191209_14_46_58';
   '200130_13_21_13 FunctAcq';
   '200130_13_36_14 FunctAcq';
   '200130_13_49_09 FunctAcq';
   %'200130_14_02_12 FunctAcq';
   '200130_14_15_24 FunctAcq';
   '200130_14_29_30 FunctAcq';
    };
savepath = 'X:\MFB';
%% Plot by 1-100 modes
plot_decode_all_modes = 1;
if plot_decode_all_modes == 1
    Cvev_loco = cell(length(folder_names), 1);
    Cvev_wsk = cell(length(folder_names), 1);
    Cvev_state = cell(length(folder_names), 1);

    for file_i = 1:length(folder_names)
        file = char(folder_names(file_i));
        quickAnalysis;

        Y = MI_wheel_r';
        Y2 = MI_whisker_r';
        X = dff_rz';
        Y3 = L_state';

        X = X(valid_t,:);
        Y = Y(valid_t);
        Y2 = Y2(valid_t);
        Y3 = Y3(valid_t);

        % PCA
        [coeff, score, ~, ~, ~] = pca(X);

        NC = 1:100;
        mean_cvev_loco_all = zeros(1, length(NC));
        mean_cvev_wsk_all = zeros(1, length(NC));
        mean_cvev_state_all = zeros(1, length(NC));

        for sj = 1:length(NC)
            N = NC(sj);
            X_n = score(:, 1:N);

            cvp = cvpartition(size(X_n, 1), 'KFold', 5);
            lambda = 10;

            cvev_loco = zeros(cvp.NumTestSets, 1);
            cvev_wsk = zeros(cvp.NumTestSets, 1);
            cvev_state = zeros(cvp.NumTestSets, 1);

            for f = 1:cvp.NumTestSets
                trainIdx = cvp.training(f);
                testIdx = cvp.test(f);

                X_train = X_n(trainIdx, :);
                Y_train_all = {Y(trainIdx), Y2(trainIdx), Y3(trainIdx)};
                X_test = X_n(testIdx, :);
                Y_test_all = {Y(testIdx), Y2(testIdx), Y3(testIdx)};

                [CVEV, ~, ~] = myRidge({X_train, X_train, X_train}, Y_train_all, ...
                                       {X_test, X_test, X_test}, Y_test_all, lambda);

                cvev_loco(f) = CVEV(1);
                cvev_wsk(f) = CVEV(2);
                cvev_state(f) = CVEV(3);
            end

            mean_cvev_loco_all(sj) = mean(cvev_loco);
            mean_cvev_wsk_all(sj) = mean(cvev_wsk);
            mean_cvev_state_all(sj) = mean(cvev_state);
        end

        Cvev_loco{file_i} = mean_cvev_loco_all;
        Cvev_wsk{file_i} = mean_cvev_wsk_all;
        Cvev_state{file_i} = mean_cvev_state_all;
    end
end

%% Save the results
save([savepath '\decode\results_1-100.mat'], 'Cvev_loco','Cvev_wsk','Cvev_state');

%% Decode by optimal modes/selected modes

for file_i = 1:length(folder_names)

    file = char(folder_names(file_i));
    quickAnalysis;

    Y = MI_wheel_r';
    Y2 = MI_whisker_r';
    Y3 = L_state';

    X = dff_rz';
    [coeff, score, ~, ~, ~] = pca(X);

    N1 = 10;
    N2 = 10;
    N3 = 10;

    X_n1 = score(:, 1:N1);
    X_n2 = score(:, 1:N2);
    X_n3 = score(:, 1:N3);

    cvp = cvpartition(size(X_n1, 1), 'KFold', 5);
    lambda = 10;

    cvev_loco = zeros(cvp.NumTestSets, 1);
    cvev_wsk = zeros(cvp.NumTestSets, 1);
    cvev_state = zeros(cvp.NumTestSets, 1);

    for f = 1:cvp.NumTestSets

        trainIdx = cvp.training(f);
        testIdx = cvp.test(f);

        Y_train_all = {Y(trainIdx), Y2(trainIdx), Y3(trainIdx)};
        Y_test_all = {Y(testIdx), Y2(testIdx), Y3(testIdx)};

        [CVEV, ~, Y_pred] = myRidge({X_n1(trainIdx, :), X_n2(trainIdx, :), X_n3(trainIdx, :)}, Y_train_all, ...
                               {X_n1(testIdx, :), X_n2(testIdx, :), X_n3(testIdx, :)}, Y_test_all, lambda);

        cvev_loco(f) = CVEV(1);
        cvev_wsk(f) = CVEV(2);
        cvev_state(f) = CVEV(3);
    end


    Cvev_loco{file_i} = mean(cvev_loco);
    Cvev_wsk{file_i} = mean(cvev_wsk);
    Cvev_state{file_i} = mean(cvev_state);

    Actual_loco_all{file_i} = Y(testIdx);
    Actual_whisk_all{file_i} = Y2(testIdx);
    Actual_state_all{file_i} = Y3(testIdx);
    Decode_loco_all{file_i} = Y_pred(1);
    Decode_whisk_all{file_i} = Y_pred(2);
    Decode_state_all{file_i} = Y_pred(3);
end



%% Save the results
save([savepath '\decode\results_optimal_modes.mat'], 'Cvev_loco','Cvev_wsk','Cvev_state');

%% plot test vs pred
 
k = 13; % choose a datset

set(figure, 'Position', [400, 100, 500, 190]);

% Plot Loco
cell1 = Actual_loco_all{k};
minL = min(cellfun(@(x) numel(x), cell1));
Data = zeros(minL, numel(cell1));
for i = 1:numel(cell1)
    Data(:, i) = cell1{i}(1:minL);
end
average_test = mean(Data, 2);

cell2 = Decode_loco_all{k};
minL = min(cellfun(@(x) numel(x), cell2));
Data = zeros(minL, numel(cell2));
for i = 1:numel(cell2)
    Data(:, i) = cell2{i}(1:minL);
end
average_pred = mean(Data, 2);

subplot(3, 1, 1);
h1 = plot(average_test, 'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;
h2 = plot(average_pred,'color',[173,210,157]/255,'LineWidth',1);
disp(['Loco Cvev=' num2str(Cvev_loco{k}, '%.2f')])

xlabel('Time');
ylabel('Corr')
box('off')
axis off;

% Plot Wsk
cell3 = Actual_whisk_all{k};
minL = min(cellfun(@(x) numel(x), cell3));
Data = zeros(minL, numel(cell3));
for i = 1:numel(cell3)
    Data(:, i) = cell3{i}(1:minL);
end
average_test = mean(Data, 2);

cell4 = Decode_whisk_all{k};
minL = min(cellfun(@(x) numel(x), cell4));
Data = zeros(minL, numel(cell4));
for i = 1:numel(cell4)
    Data(:, i) = cell4{i}(1:minL);
end
average_pred = mean(Data, 2);
subplot(3, 1, 2);
h3 = plot(average_test, 'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;
h4 = plot(average_pred,'color',[255,163,26]/255,'LineWidth',1);

disp(['Whisk Cvev=' num2str(Cvev_wsk{k}, '%.2f')])
box('off')
axis off;

% Plot states
cell5 = Actual_state_all{k}{end};
average_test = cell5;

cell6 = Decode_state_all{k};
minL = min(cellfun(@(x) numel(x), cell6));
Data = zeros(minL, numel(cell6));
for i = 1:numel(cell6)
    Data(:, i) = cell6{i}(1:minL);
end
average_pred = mean(Data, 2);
subplot(3, 1, 3);
h5 = plot(average_test, 'color',[0.6,0.6,0.6],'LineWidth',1); 
hold on;
h = plot(average_pred,'color',[255,53,255]/255,'LineWidth',1);

ax = gca;
xlims = ax.XLim;
ylims = ax.YLim;
x0 = xlims(2) - 1000; % 10s
y0 = ylims(1);
plot(ax, [x0, x0+1000], [y0, y0]*2, 'k-', 'LineWidth', 3);
text(ax, x0, y0 - diff(ylims) * 0.25, '10 s', 'FontSize', 13, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

disp(['state Cvev=' num2str(Cvev_state{k}, '%.2f')])
box('off')
axis off;

currentDate = datestr(now, 'yyyy-mm-dd');
savepath2 = fullfile(savepath, 'Figures', 'Figure5', currentDate);
if ~exist(savepath2, 'dir')
    mkdir(savepath2);
end

fileName = [ 'Ridge_' file '_mode# ' num2str(N) '.png'];
fullFilePath = fullfile(savepath2, fileName);
print(fullFilePath, '-dpng', '-r300');

