
run '/Users/sadra/Desktop/MFB/Init/Initialize';

path_MI = [path_home '/MI_Data'];
path_save = [path_home '/Data_processed'];
dlc_res_path = [path_home '/DeepLC/TrackingData/'];

goodFolders__wheelMI;

folder_names = folder_names_point;
    
folder_names = {'171003_14_05_44'};

%% which info to load and preprocess
load_header_info = 0;
load_speed_data = 0;

load_functional_data = 0;
compute_spatial_clusters = 0;
clean_for_MFs = 0;
add_loc_MFclust = 0;

load_MI_wheel = 1;
load_MI_whisking = 1;
load_MI_tail = 0;
load_MI_forelimb = 0;

load_DLC = 1;

load_previous_prep_data = 1;
save_prep_data = 1;
    
%% --

for ij = 1:length(folder_names) 
    folder_name = folder_names{ij}(1:15)
    
    prep_name = ['processedData_' folder_name '.mat'];
    
    path_data = [path_home '/Data_raw/' folder_name];
    path_header = path_data;
    speed_path = [path_data '/Speed_Data'];
    
    Tshift = 0;
    laser_times;
    Tshift
    
    %% load previous
    clear prepData;
    if load_previous_prep_data
        cd(path_save);
        if exist(prep_name) ~= 0; load(prep_name); end
    end
    
    %% header Info
    if load_header_info
        %'header info ...'
        
        read__HeaderInfo_pointExp;

        % prepData.header
        prepData.header.dt = dt;
        prepData.header.T = T;
        prepData.header.Ntr = Ntr;
        prepData.header.Nroi = Nroi;
        prepData.header.Npt = Npt;
        prepData.header.Tcycl = Tcycl;
        prepData.header.Ncycl = Ncycl;
        prepData.header.Troi = Troi;
        prepData.header.pixel_size = pixel_size;
        prepData.header.XYZ = [X0 Y0 Z0];
        prepData.header.XYZ_norm = [Xn Yn Zn];

    end
    
    dt = prepData.header.dt;
    T = prepData.header.T;
    Ntr = prepData.header.Ntr;
    Nroi = prepData.header.Nroi;
    Npt = prepData.header.Npt;
    Tcycl = prepData.header.Tcycl;
    Ncycl = prepData.header.Ncycl;
    Troi = prepData.header.Troi;
    XYZ = prepData.header.XYZ;
    %XYZ_norm = prepData.header.XYZ_norm;
        
    %% 
    if compute_spatial_clusters
        
        if ~load_header_info
        X0 = XYZ(:,1);Y0 = XYZ(:,2);Z0 = XYZ(:,3); 
        end
        compute__Clusters;
                
        prepData.clusters = ROI_clusters;
    end
    
    %% 
    if load_speed_data
        %'speed data ...'
        
        read__SpeedData;

        % prepData.speed
        prepData.speed.spDall = speed;
        prepData.speed.spdraw = spdraw;
        prepData.speed.Tb_tr = Tb_tr;
        prepData.speed.Tf_tr = Tf_tr;
        prepData.speed.Tint = Tint;
        prepData.speed.dt_encoder = dt_encoder;
    
    end
    
    %%
    
    if load_MI_whisking || load_MI_tail || load_MI_forelimb || load_MI_wheel || load_DLC
       Tb_tr = prepData.speed.Tb_tr;
       Tf_tr = prepData.speed.Tf_tr;
    end
    
    if load_MI_wheel  
        %'wheel MI ...'
        file_name = [path_MI, '/MI_wheel', folder_name, ' FunctAcq.mat'];
        prepData.MIs.wheel = read_MIdata(file_name, Ntr, T, Tb_tr, Tf_tr, Tshift);
    end
    
    if load_MI_whisking  
        %'whisking MI ...'
        file_name = [path_MI, '/MI_', folder_name, ' FunctAcq.mat'];
        prepData.MIs.whisking = read_MIdata(file_name, Ntr, T, Tb_tr, Tf_tr, Tshift);
    end

    % -
    if load_MI_tail
        %'tail MI ...'
        file_name = [path_MI, '/MI_Tail', folder_name, ' FunctAcq.mat'];
        prepData.MIs.body = read_MIdata(file_name, Ntr, T, Tb_tr, Tf_tr, Tshift);
    end
    
    % -
    if load_MI_forelimb
        %'forelimb MI ...'
        file_name = [path_MI, '/MI_RightForelimb', folder_name, ' FunctAcq.mat'];
        prepData.MIs.forelimb = read_MIdata(file_name, Ntr, T, Tb_tr, Tf_tr, Tshift);
    end
    
    %% 
    if load_DLC
        % - video times
        dlc_data_path = [dlc_res_path 'EyeCam-relative times_' folder_name '.txt'];
        if exist(dlc_data_path, 'file') == 2
            zt = importdata(dlc_data_path);
            tt = (zt(2:end,2)-Tshift)/1e3;
            % - forelimb positions from DeepLabCut
            zd = importdata([dlc_res_path 'EyeCam_' folder_name 'DeepCut_resnet50_HanaTestApr24shuffle1_50000.csv']);
            xy_L = zd.data(:,2:3);
            dis_L = fastsmooth(sqrt(nansum(diff(xy_L).^2,2)), 10,3,1);
            xy_R = zd.data(:,5:6);
            dis_R = fastsmooth(sqrt(nansum(diff(xy_R).^2,2)), 10,3,1);

            fl_mv0 = (dis_L+dis_R)/2;
            fl_mv = trial_resolved_trans(tt*1e3, fl_mv0, Tb_tr, Tf_tr, Ntr, T);

            fl_con0 = fastsmooth((zd.data(:,4)+ zd.data(:,7))/2, 10,3,1);
            fl_con = trial_resolved_trans(tt*1e3, fl_con0, Tb_tr, Tf_tr, Ntr, T);

            prepData.DLC.fl_mv = fl_mv;
            prepData.DLC.fl_con = fl_con;
        end
    end
    
    %%
    if load_functional_data
        %'functional data ...'
        
        read__FunctionalData_pointExp;

        % artefact to be removed
        if strcmp(folder_name, '171212_13_36_10'); Gdff(1,:,T < 3e3) = nan; end

        prepData.functional.Gdff = Gdff;
    end
    
    %% MF df/f
    if clean_for_MFs
        cc_threshold = 0.5;

        ROI_clusters = prepData.clusters;
        cl_no = length(ROI_clusters);

        dff = prepData.functional.Gdff; 
        
        % -- clean clusters for points with low functional correlation
        for i = 1:cl_no
            zz = dff(:,ROI_clusters{i},:);

            zp = permute(zz, [2,3,1]); 
            zr = reshape(zp, size(zp,1),[]);

            ccw = corrcoef(zr', 'rows', 'complete'); 
            ccw(eye(length(ccw))==1) = nan; 
            ccw(ccw < cc_threshold) = nan;

            exid = find(isnan(nanmean(ccw)));
            ROI_clusters{i}(exid)=[];
        end
        prepData.clusters_functional = ROI_clusters;
        
        % -- add df/f for MFs (average activity of MFs within clusters)
        dff_mf = [];
        for i = 1:cl_no
            if ~isempty(ROI_clusters{i})
                dff_mf = cat(2, dff_mf, nanmean(dff(:,ROI_clusters{i},:),2));
            end
        end
        prepData.functional.Gdff_MFs = dff_mf; 
    end
    
    if add_loc_MFclust
        % - location of clustered MFs
        MF_ids = prepData.clusters_functional;
        XYZ = prepData.header.XYZ;
        Nmf = size(prepData.functional.Gdff_MFs,2);
        MF_xyz = nan*ones([Nmf,3]);
        cnt = 0;
        for i = 1:length(MF_ids)
            if ~isempty(MF_ids{i})
                cnt = cnt+1;
                MF_xyz(cnt,:) = nanmean(XYZ(MF_ids{i},:),1);
            end
        end
        prepData.header.XYZ_MF_clust = MF_xyz;
    end
    
    %%
    if save_prep_data
        cd(path_save);

        save(prep_name, 'prepData');
    end

end
