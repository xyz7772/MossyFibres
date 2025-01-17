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

%Select
folder_names = {
   '171212_16_19_37';
   '191018_13_39_41';
   %'191018_13_56_55';
   '191018_14_30_00';
   %'191018_14_11_33';
   '191209_13_44_12';
   %'191209_14_04_14';
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

% for no puff
folder_names = {
   %'171212_16_19_37';
   %'191018_13_39_41';
   %'191018_13_56_55';
   '191018_14_30_00';
   %'191018_14_11_33';
   %'191209_13_44_12';
   %'191209_14_04_14';
   '191209_14_32_39';
   '191209_15_01_22';
   '191209_14_18_13';
   '191209_14_46_58';
   %'200130_13_21_13 FunctAcq';
   '200130_13_36_14 FunctAcq';
   '200130_13_49_09 FunctAcq';
   %'200130_14_02_12 FunctAcq';
   '200130_14_15_24 FunctAcq';
   %'200130_14_29_30 FunctAcq';
    };

varexp=cell(length(folder_names),1);
numcomp=zeros(length(folder_names),1);
c_bill=[.1,.5,.8];
c_jeremy=[1,.7,0];
c_bernie=[.7,.1,.7];
c_nigel=[1,.4,.4];
dataset_colors = {c_nigel,c_jeremy,c_jeremy,c_bernie,c_bernie,c_bernie,c_bernie,c_bernie,c_bill,c_bill,c_bill,c_bill,c_bill};
mice_colors = containers.Map({'Bernie', 'Bill', 'Jeremy', 'Nigel'}, ...
                                 {c_bernie, c_bill, c_jeremy, c_nigel});

%% fig. 4C cross validated explained variance 
% This section will calculate CVEV of all datasets and save varexp.mat to
% a directory
ve = cell(length(folder_names),1);
varmax =[];
dimmax =[];
figure;

for dataset_i = 1:length(folder_names)

    file=char(folder_names(dataset_i));
    quickAnalysis_Yizhou;

    mice = getMice(folder_names{dataset_i}, Jeremy_names, Bernie_names, Bill_names, Nigel_names);
    mice_colors_all{dataset_i} = mice_colors(mice);
    zz = dff_rz;
    zz(isnan(zz)) = 0;
    N_sub = size(zz,1);
    nSamples = size(dff_rz,1);    
    Nsub = 100;

    % calc cvev
    [ve{dataset_i}, dm, vm] = get_dim(zz, Nsub,  Nsub, 100);
    
    ve_mean = nanmean(ve{dataset_i}, 1);
    ve_sem = nanstd(ve{dataset_i}, 1) / sqrt(size(ve{dataset_i},1));
    
    varexp{dataset_i}=ve_mean; 
    numcomp(dataset_i)=dm;
    varmax = [varmax vm];
    dimmax = [dimmax dm];
    varexp_all{dataset_i}=ve{dataset_i};

    % plot
        
    x=1:nSamples;
    y_mean = varexp{dataset_i};
    y_sem = ve_sem;
    face_color = mice_colors_all{dataset_i};
    y_top = y_mean + y_sem;
    y_bot = y_mean - y_sem;
    
    ix = find(~isnan(y_top) & ~isnan(y_bot));
    x = x(ix);
    y_mean = y_mean(ix);
    y_top = y_top(ix);
    y_bot = y_bot(ix);
    
    fill([x, fliplr(x)], [y_top, fliplr(y_bot)], face_color, 'LineStyle', 'none', 'FaceAlpha', 0.3);
    hold on;
    plot(x, y_mean, 'Color', 'k', 'LineWidth', 2);
    
    set(gca, 'Box', 'off', 'FontSize', 18);
        [varmax(dataset_i),dm(dataset_i)] = max(nanmean(varexp{dataset_i},1));
    plot(dm(dataset_i),varmax(dataset_i)+.03,'v','MarkerFaceColor',face_color,'MarkerEdgeColor',face_color);
    
    if strcmp(mice, 'Bernie')
        h(1) = plot(nan, nan, 'Color', face_color,'LineWidth', 2);
    elseif strcmp(mice, 'Bill')
        h(2) = plot(nan, nan, 'Color', face_color,'LineWidth', 2);
    elseif strcmp(mice, 'Jeremy')
        h(3) = plot(nan, nan, 'Color', face_color,'LineWidth', 2);
    elseif strcmp(mice, 'Nigel')
        h(4) = plot(nan, nan, 'Color', face_color,'LineWidth', 2);
    end

end

xlabel('Number of components')
ylabel('Variance explained (cross-val)')
ylim([0, .5])
yticks([0 0.25 0.5])

mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
folderPath = fullfile(mfbFolderPath, 'Figures', 'Figure4', currentDate);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

fileName = ['PCA_cvev'];
fullFilePath = fullfile(folderPath, fileName);
print(fullFilePath, '-dpng', '-r300');

save("X:\MFB\MFB_AH_2023\varexp.mat", 'varexp','varexp_all','varmax','dimmax');


%% fig. 4ab
% Calculate the PCA of datasets with clear behaviour to obtain the first 20 normalized eigenvalues
% and the loadings of the first 3 PCs.

folder_names = {
   '171212_16_19_37';
   '191018_13_39_41';
   '191018_14_30_00';
   '191209_13_44_12';
   '191209_14_32_39';
   '191209_15_01_22';
   '191209_14_18_13';
   '191209_14_46_58';
   '200130_13_21_13 FunctAcq';
   '200130_13_36_14 FunctAcq';
   '200130_13_49_09 FunctAcq';
   '200130_14_15_24 FunctAcq';
   '200130_14_29_30 FunctAcq';
};

sig_level = 2;
all_loadings = [];
eigen_20 = [];


for dataset_ix = 1:length(folder_names)
    file = char(folder_names(dataset_ix));
    quickAnalysis;

    zz = dff0_rz;
    nan_times = any(isnan(zz), 1);  
    zz = zz(:, ~nan_times);
    L_state = L_state(~nan_times);
    dff_rz = dff_rz(:, ~nan_times);

    [coeff, score, latent, tsquared, explained] = pca(zz);

    if size(coeff, 2) < 3
        continue;
    end

    total_variance = sum(latent);
    ne = latent / total_variance;
    eigen_20 = [eigen_20, ne(1:20)];

    loadings = coeff(:, 1:3) .* sqrt(latent(1:3))';
    all_loadings = [loadings;all_loadings];
end


figure(Position=[200 100 880 300]); hold on;
for ix = 1:size(eigen_20,2)
    face_color = dataset_colors{ix};    
    plot(eigen_20(:,ix),'color',face_color,'linewidth',2)
    
    mean_eigen = nanmean(eigen_20,2);
    sem_mean_eigen = nanstd(eigen_20,0,2) / sqrt(size(eigen_20,2));
    errorbar(1:20, mean_eigen, sem_mean_eigen, 'k', 'linestyle', 'none', 'LineWidth', 2);
    bar(1:20, mean_eigen, 'FaceAlpha',0, 'LineWidth',2)
    xlabel('Num. components')
    ylabel('Normalized eigenvalue')
    set(gca, 'FontSize', 15, 'XTick', 1:20)
    ylim([0, .2])
    yticks([0 0.1 0.2])
end


fileName = ['PCA_normalised_eigenvalue'];
fullFilePath = fullfile(folderPath, fileName);
print(fullFilePath, '-dpng', '-r300');


figure;
% plot loading
for i = 1:3
    subplot(3, 1, i);
    histogram(all_loadings(:, i),'Binwidth',0.05,'Normalization', 'probability', 'FaceColor', 'b');
    xlabel(sprintf('PC%d', i));
    ylabel('Prob.');
    box("off")
    set(gca, 'FontSize', 15)
    set(gcf, 'Position', [100, 100, 250, 600]);
    xlim([-1,1])
end


fileName = ['PCA_loadings'];
fullFilePath = fullfile(folderPath, fileName);
print(fullFilePath, '-dpng', '-r300');


%% 4d
load("X:\MFB\MFB_AH_2023\varexp.mat")

c_bill=[.1,.5,.8];
c_jeremy=[1,.7,0];
c_bernie=[.7,.1,.7];
c_nigel=[1,.4,.4];
dataset_colors = {c_nigel,c_jeremy,c_jeremy,c_jeremy,c_bernie,c_bernie,c_bernie,c_bernie,c_bernie,c_bernie,c_bill,c_bill,c_bill,c_bill,c_bill,c_bill};


varmax_all = []; dimmax_all = []; col_all = [];
for dataset_ix = 1:length(varexp)
    if ~isempty(varexp{dataset_ix})
        for k = 1:10
            [varmax_,dimmax_] = max(varexp_all{dataset_ix}(k,:));
            varmax_all = [varmax_all; varmax_];
            dimmax_all = [dimmax_all; dimmax_];
        end
    end
end

figure('Position', [30 20 650 350])
ix = find(~isnan(varmax_all));
plot(varmax_all,dimmax_all,'vk','MarkerFaceColor',[.6,.6,.6],'MarkerEdgeColor',[.6,.6,.6])
hold on, plot([0,1],[0,1]*(varmax_all(ix)'*varmax_all(ix))\(varmax_all(ix)'*dimmax_all(ix)),'k')

xlabel('Max variance explained')
ylabel('Lower bound of dimensionality')
set(gca,'FontSize',15)
box('off')
yticks([0 40 80])
xticks([0 0.5 1])

for dataset_ix = 1:length(varexp)
    if ~isnan(varmax(dataset_ix))
        plot(varmax(dataset_ix),dimmax(dataset_ix),'vk','MarkerFaceColor',dataset_colors{dataset_ix},'MarkerEdgeColor',dataset_colors{dataset_ix})
    end
end
mfbFolderPath = 'X:\MFB';
currentDate = datestr(now, 'yyyy-mm-dd');
folderPath = fullfile(mfbFolderPath, 'Figures', 'Figure4', currentDate);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

fileName = ['PCA_low_bound'];
fullFilePath = fullfile(folderPath, fileName);
print(fullFilePath, '-dpng', '-r300');

%% 4e

subpopulation_plot


