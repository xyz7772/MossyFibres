
addpath([path_home '/Analysis/Analysis_pipeline/Utilities/matlab-networks-toolbox-master']);

%% Parameters

% threshold on the radius of ROIs belonging to each cluster
%cluster_radius = 10; % Melania's map
%cluster_radius = 5; % other
cluster_radius = 5;
% except for Melania
switch folder_name(1:6) 
    case '161215'; cluster_radius = 10;
    case '161221'; cluster_radius = 10;
    case '161222'; cluster_radius = 10;
end
        
cluster_diameter = 2*cluster_radius;

% the permissiveness of the initial clustering (don't change if unsure)
exp_sharpness = 5;%cluster_diameter;

% minimum no of points to be counted as a cluster
min_roi_cluster = 3;

%% Clustering
      
Nroi = size(X0,1);

% pairwise distance of ROIs wrt each other
DD0 = zeros(Nroi, Nroi);
for i = 1:Nroi
    for j = i:Nroi
        d = sqrt((X0(i)-X0(j)).^2 + (Y0(i)-Y0(j)).^2 + ((Z0(i)-Z0(j))).^2);
        DD0(i,j) = d; DD0(j,i) = d;
    end
end

% the distance matrix to infer the clusters from
DD = exp(-DD0 /exp_sharpness); 

% infer the initial clusters / modules 
clstrs = newmanEigenvectorMethod(DD);

% how many have you inferred?
clst_N0 = size(clstrs,2);

%% Refining

%clear ROI_clusters;
ROI_clusters = [];

centers = [];
clst_cnt = 0;
for i = 1:clst_N0
    
    % get the ids of ROIs belonging to the i-th cluster
    ids = clstrs{1,i};
    
    % compute the center of the cluster and the distance of ROIs to it
    Xm0 = mean(X0(ids)); Ym0 = mean(Y0(ids)); Zm0 = mean(Z0(ids));
    dx = X0(ids)-Xm0; dy = Y0(ids)-Ym0; dz = Z0(ids)-Zm0;
    dd0 = sqrt(dx.^2+dy.^2+dz.^2);
    
    % exclude ROIs farther than the cluster raidus
    ids(find(dd0 >= cluster_radius)) = [];
    
    % compute the new center of the cluster
    Xm = mean(X0(ids)); 
    Ym = mean(Y0(ids)); 
    Zm = mean(Z0(ids));
    
    % a criterion to make sure there is a minimum distance between centers
    if i == 1
        crit = 0;
    else
        crit = sum(sqrt(sum((([Xm Ym Zm] - centers).^2),2)) < cluster_diameter);
    end
    centers = [centers; [Xm Ym Zm]];
    
    if length(ids) > 0
        % compute the distance of all ROIs to the center
        Dall = sqrt((X0-Xm).^2 + (Y0-Ym).^2 + (Z0-Zm).^2);
        % is there any ROI closeby who was missed? include it!
        nids = find(Dall <= 1*cluster_radius);
        ids = unique([ids transpose(nids)]);
    end;
        
    % include the clusters if the size is larger than one 
    % and there is no other too-close-cluster found so far 
    if length(ids) >= min_roi_cluster && crit == 0
        clst_cnt = clst_cnt + 1;
        ROI_clusters{clst_cnt} = ids;  
    end
end

clusters_no = length(ROI_clusters);

%% Sorting the clusters spatially

clear ROI_clusters_sorted;

ds = [];
for i = 1:clusters_no
    ids = ROI_clusters{i};
    % the center of the cluster
    ds = [ds (1+mean(X0(ids))) + (1-mean(Y0(ids))) * 10000 + (1-mean(Z0(ids)))*100];
end

[dss, sids] = sort(ds);

for i = 1:clusters_no
    ROI_clusters_sorted{i} = ROI_clusters{sids(i)};
end

ROI_clusters = ROI_clusters_sorted;

for i = 1:clusters_no
    ids = ROI_clusters{i};
    Xm = nanmean(X0(ids)); Ym = nanmean(Y0(ids)); Zm = nanmean(Z0(ids));

    clusterCenters{i} = [Xm; Ym; Zm];
end
