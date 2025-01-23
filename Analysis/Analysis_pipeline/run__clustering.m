
run '/Users/sadra/Desktop/MFB/Init/Initialize';

folder_names = {
'170712_15_35_44'
};

data_path0 = '/Volumes/silverlab-general/Hana/';
fig_path0 = '/Volumes/silverlab-general/Hana/all group clustering/';

calculate_clusters = 1;

plot_rawZstacks = 0;
plot_overlaid = 0;
plot_clusters = 1;
plot_clusters_overlaid = 1;
zoomed = 1;
zoomL = 50;

for ij = 1:length(folder_names)
    folder_name = char(folder_names(ij))
    
    if ~calculate_clusters
        load([path_home '/Data_processed/processedData_' folder_name '.mat']);
        ROI_clusters = prepData.clusters;
    end
    
    data_path = [data_path0, char(folder_name), ' FunctAcq'];
    path_ZstackImages = [data_path, '/Zstack Images'];
    fig_path = [fig_path0, char(folder_name)];
    
    path_header = data_path;
    read__HeaderInfo_pointExp;
    
    Xim = Xn * frame_size/2 + frame_size/2;
    Yim = Yn * frame_size/2 + frame_size/2;
    Dim = sqrt(Xim.^2 + Yim.^2);
    
    dd = abs(nanmean(diff(zz_norm)));
    zcl = (Z0 - min(Z0))/(max(Z0)-min(Z0));
    
%% Plotting
if plot_overlaid
    path_plot = strcat(fig_path,'/Overlaid');
    if ~exist(path_plot, 'dir'); mkdir(path_plot); end;
    cd(path_plot);

    for j = 1:length(zz_norm)
        %zsid = find(abs(Zn-zz_norm(j)) <= dd/2);
        zsid = find(Zn == zz_norm(j));
        %length(zsid);
        
        if length(zsid) > 0
            gimdata = imread(strcat(path_ZstackImages,'/GreenChannel_',num2str(j,'%04i'),'.tif'));
            
            fig = figure('Position', [50, 200, 1000, 600], 'visible','off');          
            title(strcat('Z: @', num2str(zz_norm(j))), 'FontSize', 10); hold on;

            imagesc(gimdata); 
            caxis([0 500]); colormap(flipud(gray));
            xlim([0, frame_size]); ylim([0, frame_size]);
         
            set(gca,'DataAspectRatio',[1 1 1])
            set(gca, 'Ydir', 'reverse');

            for kk = 1:length(zsid)
                jj = zsid(kk);
                cl = [Xim(jj)/512 Yim(jj)/512 zcl(jj)];

                plot(Xim(jj), Yim(jj), 'o', 'MarkerSize',15, 'LineWidth',1, 'MarkerEdgeColor', cl);
            end

            if ~isempty(zsid)
                h = legend(num2str(zsid), 'Location', 'northeastoutside'); 
                set(h,'FontSize',7); 
            end

            fig_name = strcat('zStackOverlaid_',num2str(j));
            print(fig_name, '-dpng');
        end
    end
    close('all');
end

%% compute clusters and plot the overall map
if plot_clusters
    path_plot = strcat(fig_path,'/Clusters');
    if ~exist(path_plot, 'dir'); mkdir(path_plot); end;
    cd(path_plot);
        
    % -- compute clusters based on spatial information
    if calculate_clusters
        compute__Clusters;
    end
    
    clusters_no = length(ROI_clusters);
    
    %clusters_info.ROIs = ROI_clusters;
    %clusters_info.centers = clusterCenters;
    %clusters_info.dist_threshold = dist_threshold;
    %clusters_info.exp_sharpness = exp_sharpness;

    save(strcat('ROI_clusters',folder_name,'.mat'), 'ROI_clusters')
    
    % -
    
    figure('Position', [100, 100, 800, 600], 'visible', 'off'); 
    hold on; %rotate3d on;
    pbaspect([1 1 1]);

    scatter3(X0, Y0, Z0, 100, 'LineWidth', 1, 'MarkerEdgeColor', 'k');
    
    for i = 1:clusters_no
        ids = ROI_clusters{i};
        Xm = mean(X0(ids)); Ym = mean(Y0(ids)); Zm = mean(Z0(ids));

        X1 = (-min(X0) + Xm)/(-min(X0)+max(X0));
        Y1 = (-min(Y0) + Ym)/(-min(Y0)+max(Y0));
        Z1 = (-min(Z0) + Zm)/(-min(Z0)+max(Z0));

        plot3(Xm, Ym, Zm, 's', 'MarkerEdgeColor', [X1 Y1 Z1], 'MarkerSize', 40, 'LineWidth',2);

        lbl = i;
        text(Xm+5, Ym+2, 0.1, num2str(lbl), 'Color', [X1 Y1 Z1], 'fontsize', 15);
        
%         if i < clusters_no/2
%             text((i-clusters_no/4)/clusters_no*4, 1.1, 0, num2str(lbl), ...
%                 'Color', [X1 Y1 Z1], 'fontsize', 15);
%         else
%             text(1.1, (i-3*clusters_no/4)/clusters_no*4, num2str(lbl), ...
%                 'Color', [X1 Y1 Z1], 'fontsize', 15);
%         end
        
        scatter3(X0(ids), Y0(ids), Z0(ids), 100, 'LineWidth', 1,'MarkerEdgeColor', [X1 Y1 Z1]);
    end

    xlabel('X'); ylabel('Y'); zlabel('Z');
    set(gca, 'LineWidth', 1, 'FontSize', 15);
    set(gca, 'Ydir', 'reverse');
    xlim([0,250]);
    ylim([0,250]);
    
    fig_name = strcat('Clusters__', folder_name);
    print(fig_name, '-dpng');   
end

%% Plotting clusters over the closest zstack image
if plot_clusters_overlaid
    
    path_plot = strcat(fig_path,'/Clusters');
    if ~exist(path_plot, 'dir'); mkdir(path_plot); end
    cd(path_plot);
    
    clusters_no = length(ROI_clusters);

    for i = 1:clusters_no
        strcat('Cluster #', num2str(i));
        ids = ROI_clusters{i};
        
        zm = nanmean(Zn(ids));
        [zs, zsid] = sort(abs(zz_norm - zm));

        zt = imread(strcat(path_ZstackImages,'/GreenChannel_',num2str(zsid(1),'%04i'),'.tif'));
    
        fig = figure('Position', [50, 200, 1000, 600], 'visible','off'); 
        title(strcat('Cluster #', num2str(i), '-- zstack:', num2str(zsid(1)), ' @',num2str(zz_norm(zsid(1)))) );
        hold on;
        
        imagesc(zt);
        caxis([0 500]); colormap(flipud(gray));
        xlim([0, frame_size]); ylim([0, frame_size]);
    
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca, 'Ydir', 'reverse');
        
        %
        Xm = mean(Xim(ids)); Ym = mean(Yim(ids)); Zm = mean(zcl(ids));
        
        for kk = 1:length(ids)
            jj = ids(kk);
            plot(Xim(jj), Yim(jj), 'o', 'MarkerSize',15, 'LineWidth',2);%, ...
            %'MarkerEdgeColor', [Xim(jj)/512 Yim(jj)/512 -min(Z0)+Z0(jj)]);
        end
            
        h = legend(num2str(transpose(ids)), 'Location', 'northeastoutside'); 
        set(h,'FontSize',12); 
        
        plot(Xm, Ym, 's', 'MarkerEdgeColor', [Xm/512 Ym/512 Zm], 'MarkerSize', 50, 'LineWidth',2);
        
        if zoomed
            Xmin = Xm-zoomL; 
            if Xmin < 0; Xmin = 0; end
            Xmax = Xm+zoomL; 
            if Xmax > frame_size; Xmax = frame_size; end
            xlim([Xmin, Xmax]);
            
            Ymin = Ym-zoomL; 
            if Ymin < 0; Ymin = 0; end
            Ymax = Ym+zoomL; 
            if Ymax > frame_size; Ymax = frame_size; end
            ylim([Ymin, Ymax]);
        else
            xlim([0, frame_size]); ylim([0, frame_size]);
        end
        
        fig_name = strcat('zStackCluster_',num2str(i));
        print(fig_name, '-dpng');
    end  
    close('all');
end

end