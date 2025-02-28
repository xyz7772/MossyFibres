%% Plot extrapolation for different subsample populations
%Jeremy
Jeremy_names = {
    '191018_13_39_41';
    '191018_13_56_55';
    '191018_14_30_00';
    '191018_14_11_33';
    };

%Bernie
Bernie_names = {
   '191209_13_44_12';
   '191209_14_04_14';
   '191209_14_32_39';
   '191209_15_01_22';
   '191209_14_18_13';
   '191209_14_46_58';
    };

%Bill
Bill_names = {
    '200130_13_21_13 FunctAcq';
    '200130_13_36_14 FunctAcq';
    '200130_13_49_09 FunctAcq';
    '200130_14_02_12 FunctAcq';
    '200130_14_15_24 FunctAcq';
    '200130_14_29_30 FunctAcq';
    };

%Nigel
Nigel_names = {
    '171212_16_19_37';
    };


folder_names = {
   '171212_16_19_37'; %1
   '191018_13_39_41'; %2
   '191018_13_56_55'; %3 inactive
   '191018_14_30_00'; %4
   '191018_14_11_33'; %5 ianctive
   '191209_13_44_12'; %6
   '191209_14_04_14'; %7 inactive
   '191209_14_32_39'; %8
   '191209_15_01_22'; %9
   '191209_14_18_13'; %10
   '191209_14_46_58'; %11
   '200130_13_21_13 FunctAcq'; %12
   '200130_13_36_14 FunctAcq'; %13
   '200130_13_49_09 FunctAcq'; %14
   %'200130_14_02_12 FunctAcq'; 
   '200130_14_15_24 FunctAcq'; %15
   '200130_14_29_30 FunctAcq'; %16
    };



%% Calculate dimensionality for chosen subpopulation size (~10 min for each subpopulation size)

basedir="X:\MFB";

for sub_i= [75 100 125 150 175]

N_sub = sub_i
N_file= length(folder_names);
acquisition_rate=100;
dimmax = nan(N_file,1);
varexp = cell(N_file,1);

for dataset_i = 1:N_file
    % Load data
    file=char(folder_names(dataset_i));

    quickAnalysis;
    zz = dff_rz;
    zz(isnan(zz)) = 0;
    if size(zz,1) >= N_sub
        % Calculate dimensionality for grouped axons
        [varexp{dataset_i},dimmax(dataset_i),~] = get_dim(zz,N_sub,N_sub,acquisition_rate);
    end

end

save([char(basedir),'/Processed/dimensionality_N',num2str(N_sub)],'dimmax','varexp')
end


%% fig 4e

valid_dataset = [1,2,4,6,8:16];
Jeremy_ids = [2,3,4,5];
Bernie_ids = [6:11];
Bill_ids = [12:16];

c_bill=[.1,.5,.8];
c_jeremy=[1,.7,0];
c_bernie=[.7,.1,.7];
c_nigel=[1,.4,.4];
mice_colors = containers.Map({'Bernie', 'Bill', 'Jeremy', 'Nigel'}, ...
                                 {c_bernie, c_bill, c_jeremy, c_nigel});

N_sub = [75:25:175];
basedir="X:\MFB";
slope = zeros(size(N_sub));
sem=[];
N_file= length(folder_names);
for k = 1:length(N_sub)
    % Recalculate max variance and dimensionality
    varmax = nan(N_file,1);
    dimmax = nan(N_file,1);
    slope_values = nan(N_file,1);
    for dataset_i = valid_dataset
        load(fullfile(basedir, 'Processed/MFA subpopulation dim/', ['dimensionality_N',num2str(N_sub(k))]))
        if ~isempty(varexp{dataset_i})
            [varmax(dataset_i),dimmax(dataset_i)] = max(nanmean(varexp{dataset_i},1));
        end
        slope_values(dataset_i) = (varmax(dataset_i)'*varmax(dataset_i))\(varmax(dataset_i)'*dimmax(dataset_i));
    end
    slope_ind = N_sub(k)./slope_values;
    ix = find(~isnan(varmax));
    slope(k) = (varmax(ix)'*varmax(ix))\(varmax(ix)'*dimmax(ix));
    sem(k) = nanstd(slope_ind) / sqrt(sum(~isnan(slope_ind)));
end


figure('Position', [30 20 500 250])
bar(N_sub', N_sub./slope, 'FaceColor', [.6, .6, .6], 'EdgeColor', 'w', 'LineWidth', 1)
hold on;
errorbar(N_sub', N_sub./slope, sem, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
set(gca, 'FontSize', 13, 'Box', 'off')
xlim([min(N_sub)-10, max(N_sub)+10])

xlabel('Number of MFA')
ylabel('MFA per dimension')

% print
mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
folderPath = fullfile(mfbFolderPath, 'Figures', 'Figure4', currentDate);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

fileName = ['subgroups_' char("errorbar") '_corr'];
fullFilePathPDF = fullfile(folderPath, [fileName,'.pdf']);
exportgraphics(gcf, fullFilePathPDF, 'ContentType', 'vector');



%% color set

for dataset_i = 1:length(folder_names)
    mice = getMice(folder_names{dataset_i}, Jeremy_names, Bernie_names, Bill_names, Nigel_names);
    mice_colors_all{dataset_i} = mice_colors(mice);
end

%% Extrapolate to estimate full dimensionality

load("X:\MFB\MFB_AH_2023\varexp.mat")
varmax_all = []; dimmax_all = []; col_all = [];
for dataset_ix = 1:length(varmax)
    if ~isempty(varexp_all{dataset_ix})
    for k = 1:10
        [varmax_,dimmax_] = max(varexp_all{dataset_ix}(k,:));
        varmax_all = [varmax_all; varmax_];
        dimmax_all = [dimmax_all; dimmax_];
    end
    end
end
figure, 
ix = find(~isnan(varmax_all));
plot(varmax_all,dimmax_all,'vk','MarkerFaceColor',[.6,.6,.6],'MarkerEdgeColor',[.6,.6,.6])
hold on, plot([0,1],[0,1]*(varmax_all(ix)'*varmax_all(ix))\(varmax_all(ix)'*dimmax_all(ix)),'k')

for dataset_ix = length(varmax)
    if ~isnan(varmax(dataset_ix))
        plot(varmax(dataset_ix),dimmax(dataset_ix),'vk','MarkerFaceColor',mice_colors_all{dataset_ix},'MarkerEdgeColor',mice_colors_all{dataset_ix})
        hold on;
    end
end
xlabel('Max variance explained')
ylabel('Number of components')
set(gca,'FontSize',18)

%ylim([0 80])

fileName = ['subgroups_' char("all") '_linear.png'];
fullFilePath = fullfile(folderPath, fileName);
% save
print(fullFilePath, '-dpng', '-r300');
%% another way to plot 

c_bill=[.1,.5,.8];
c_jeremy=[1,.7,0];
c_bernie=[.7,.1,.7];
c_nigel=[1,.4,.4];
dataset_colors = {c_nigel,c_jeremy,c_jeremy,c_bernie,c_bernie,c_bernie,c_bernie,c_bernie,c_bill,c_bill,c_bill,c_bill,c_bill};

varmax_all = []; dimmax_all = []; col_all = [];
for dataset_ix = 1:13
    if ~isempty(varexp{dataset_ix})
    for k = 1:10
        [varmax_,dimmax_] = max(varexp_all{dataset_ix}(k,:));
        varmax_all = [varmax_all; varmax_];
        dimmax_all = [dimmax_all; dimmax_];
    end
    end
end

figure('Position', [30 20 650 300])
ix = find(~isnan(varmax_all));
plot(varmax_all,dimmax_all,'vk','MarkerFaceColor',[.6,.6,.6],'MarkerEdgeColor',[.6,.6,.6])
hold on, plot([0,1],[0,1]*(varmax_all(ix)'*varmax_all(ix))\(varmax_all(ix)'*dimmax_all(ix)),'k')


xlabel('Max variance explained')
ylabel('Number of components')
set(gca,'FontSize',15)
box('off')
%yticks([0 40 80])
xticks([0 0.5 1])

for dataset_ix = 1:13
    if ~isnan(varmax(dataset_ix))
        plot(varmax(dataset_ix),dimmax(dataset_ix),'vk','MarkerFaceColor',dataset_colors{dataset_ix},'MarkerEdgeColor',dataset_colors{dataset_ix})
    end
end


fileName = ['subgroups_' char("Nigel") '_linear.png'];
fullFilePath = fullfile(folderPath, fileName);
% save
print(fullFilePath, '-dpng', '-r300');



%% Effective dimensionality
% Select
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

D_eff_all=[];

for file_i  =1:length(folder_names)
    
 file = char(folder_names(file_i));
        quickAnalysis_Yizhou;

cov_matrix = cov(dff_rz','partialrows'); 

[eigenvectors, eigenvalues_matrix] = eig(cov_matrix);
eigenvalues = diag(eigenvalues_matrix);

lambda_hat_sq = (eigenvalues.^2) ./ (sum(eigenvalues).^2);

D_eff = 1 / sum(lambda_hat_sq);

disp(['Effective Dimensionality (D_eff): ', num2str(D_eff)]);

D_eff_all = [D_eff_all D_eff];
end

%% Effective dimensionality with smoothing
% Select
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

D_eff_all = [];

window_size = 5;

for file_i = 1:length(folder_names)
    
    file = char(folder_names(file_i));
    quickAnalysis;
    
    dff_rz_smoothed = smoothdata(dff_rz, 2, 'movmean', window_size);
    
    cov_matrix = cov(dff_rz_smoothed','partialrows'); 

    [eigenvectors, eigenvalues_matrix] = eig(cov_matrix);
    eigenvalues = diag(eigenvalues_matrix);

    lambda_hat_sq = (eigenvalues.^2) ./ (sum(eigenvalues).^2);
    D_eff = 1 / sum(lambda_hat_sq);

    disp(['Effective Dimensionality (D_eff): ', num2str(D_eff)]);
    D_eff_all = [D_eff_all D_eff];
end



