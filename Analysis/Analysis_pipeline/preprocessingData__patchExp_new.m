
run '/Users/sadra/Desktop/MFB/Init/Initialize';

path_MI = [path_home '/MI_Data'];
path_save = [path_home '/Data_processed'];
dlc_res_path = [path_home '/DeepLC/TrackingData/'];

goodFolders__wheelMI;
folder_names = folder_names_patch;

path_data0 = [path_home '/Data_raw/'];

%  folder_names = {
%     '171212_16_19_37', % Nigel patch
% %     '171212_16_42_31', % Nigel plane
% %     '180130_14_21_41', %Francois patch
% %     '180130_15_11_06',
% %     '180130_15_28_20',
% %    '180131_14_58_58',
% %    '180131_15_21_48',
% %    '180131_15_37_10'
%     }
    
% - Jeremy
% path_data0 = [path_home '/Data_raw/Jeremy/'];
% folder_names = {%'191018_13_04_19',
%                 %'191018_14_11_33',
%                 %'191018_14_30_00'
%                 
%                 '191018_13_56_55',
%                 
%                 %'191018_13_39_41',
%                        
%                 };

% bernie
path_data0 = [path_home '/Data_raw/Bernie/'];
folder_names = {
%     '191209_13_44_12',
%'191209_14_46_58',
%'191209_14_18_13',
%'191209_14_32_39',
%'191209_15_01_22',
'191209_14_04_14'
     };

%% which info to load and preprocess
load_header_info = 1;
load_speed_data = 1;

load_functional_data = 1;

load_MI_wheel = 1;
load_MI_whisking = 1;
load_MI_tail = 0;
load_MI_forelimb = 0;

load_DLC = 0;

load_previous_prep_data = 1;
save_prep_data = 1;
    
%% --

for ij = 1:length(folder_names)
    folder_name = folder_names{ij}(1:15)
    
    prep_name = ['processedData_' folder_name '.mat'];
    
    path_data = [path_data0 folder_name];

    path_header = path_data;
    path_speed = [path_data '/Speed_Data'];
    path_functional = [path_data '/Functional_data'];
    
    Tshift = 0;
    laser_times;
    %Tshift

    prefix = 'patch';
    bogus_frames = [];
    if strcmp(folder_name, '171212_16_19_37')
        % -- bogus frames to be cut (too much movement)
        %bogus_frames=[];%1553:1602,1624:1640,2798:2875,2889:2893,3384:3495,3527:3538];
        bogus_frames = [2807:2814, 3399:3401];
        
        prefix = 'clean_patch';
    elseif strcmp(folder_name, '171212_16_42_31')
        prefix = 'clean_plane';  
    end
    
    %% load previous
    clear prepData;
    if load_previous_prep_data
        cd(path_save);
        if exist(prep_name) ~= 0; load(prep_name); end
    end
    
    %% header Info
    if load_header_info
        %'header info ...'
        
        %read__HeaderInfo_patchExp;
        read2__HeaderInfo_patchExp;
        
        % Bernie - has 13 trials instead of 20 in the header; not sure why.
        if strcmp(folder_name, '191209_13_44_12')
        Ntr = 13;
        end
        
        % prepData.header
        prepData.header.dt = dt;
        prepData.header.T = T;
        prepData.header.Ntr = Ntr;
        prepData.header.Nroi = Nroi;
        %prepData.header.Npt = Npt;
        prepData.header.Tcycl = Tcycl;
        prepData.header.Ncycl = Ncycl;
        prepData.header.Troi = Troi;
        prepData.header.pixel_size = pixel_size;
        prepData.header.XYZ_patches = [X0 Y0 Z0];
        %prepData.header.XYZ_norm = [Xn Yn Zn];
        prepData.header.Npatches = Npatches;

    end
    
    dt = prepData.header.dt;
    T = prepData.header.T;
    Ntr = prepData.header.Ntr;
    Nroi = prepData.header.Nroi;
    %Npt = prepData.header.Npt;
    Tcycl = prepData.header.Tcycl;
    Ncycl = prepData.header.Ncycl;
    Troi = prepData.header.Troi;
    pixel_size = prepData.header.pixel_size;
    %XYZ = prepData.header.XYZ;
    %XYZ_norm = prepData.header.XYZ_norm;
    Npatches = prepData.header.Npatches;
    
    % -- patch numbers to analyze
    patch_numbers = 1:Npatches;
    if strcmp(folder_name, '171212_16_19_37')
        patch_numbers = [1,2, 4:12, 14:15];
        
    elseif strcmp(folder_name, '180131_14_58_58')
        patch_numbers = [1:17,19:20];
    elseif strcmp(folder_name, '180131_15_21_48')
        patch_numbers = [1:17,19:20];
    elseif strcmp(folder_name, '180131_15_37_10')
        patch_numbers = [1:17,19:20];
    
    elseif strcmp(folder_name, '191018_13_04_19')
        patch_numbers = [1:6];
    elseif strcmp(folder_name, '191018_13_39_41')
        patch_numbers = [1:13];
    elseif strcmp(folder_name, '191018_13_56_55')
        patch_numbers = [1:13];
    elseif strcmp(folder_name, '191018_14_11_33')
        patch_numbers = [1:13];
    elseif strcmp(folder_name, '191018_14_30_00')
        patch_numbers = [1:13];
        
    elseif strcmp(folder_name, '191209_13_44_12')
        patch_numbers = [1:8];
    elseif strcmp(folder_name, '191209_14_46_58')
        patch_numbers = [1:8];
    elseif strcmp(folder_name, '191209_14_18_13')
        patch_numbers = [1:8];
    elseif strcmp(folder_name, '191209_14_32_39')
        patch_numbers = [1:8];
    elseif strcmp(folder_name, '191209_15_01_22')
        patch_numbers = [1:8];
    elseif strcmp(folder_name, '191209_14_04_14')
        patch_numbers = [1:8];    
        
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
    
    if load_MI_whisking || load_MI_tail || load_MI_forelimb || load_MI_wheel
       Tb_tr = prepData.speed.Tb_tr;
       Tf_tr = prepData.speed.Tf_tr;
    end
    
    if load_MI_wheel  
        %'wheel speed ...'
        file_name = [path_MI, '/MI_wheel', folder_name, ' FunctAcq.mat'];
        prepData.MIs.wheel = read_MIdata(file_name, Ntr, T, Tb_tr, Tf_tr, Tshift);
    end
    
    if load_MI_whisking  
        %'whisking data ...'
        file_name = [path_MI, '/MI_', folder_name, ' FunctAcq.mat'];
        prepData.MIs.whisking = read_MIdata(file_name, Ntr, T, Tb_tr, Tf_tr, Tshift);
    end

    % -
    if load_MI_tail
        %'tail data ...'
        file_name = [path_MI, '/MI_Tail', folder_name, ' FunctAcq.mat'];
        prepData.MIs.body = read_MIdata(file_name, Ntr, T, Tb_tr, Tf_tr, Tshift);
    end
    
    % -
    if load_MI_forelimb
        %'forelimb data ...'
        file_name = [path_MI, '/MI_RightForelimb', folder_name, ' FunctAcq.mat'];
        prepData.MIs.forelimb = read_MIdata(file_name, Ntr, T, Tb_tr, Tf_tr, Tshift);
    end
    
    %% 
    if load_DLC
        % - video times
        zt = importdata([dlc_res_path 'EyeCam-relative times_' folder_name '.txt']);
        tt = (zt(2:end,2)-Tshift)/1e3;
        % - forelimb positions from DeepLabCut
        %zd = importdata([dlc_res_path 'EyeCam_' folder_name 'DeepCut_resnet50_HanaTestApr24shuffle1_50000.csv']);
        zd = importdata([dlc_res_path 'EyeCam-1_' folder_name 'DeepCut_resnet50_HanaTestMay3shuffle1_100000.csv']);
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
    
    %%
    if load_functional_data
        %'functional data ...'
        
        Nrois = [];
        ROI_dff_all = []; ROI_fmb_all = [];
        f_baselines = [];
        Xall=[]; Yall=[]; Zall=[];
        for patch_no = patch_numbers
            patch_no

            read__FunctionalData_patchExp;
            Nrois = [Nrois N_miniroi];

            ROI_dff_all = cat(1, ROI_dff_all, ROI_dff);
            ROI_fmb_all = cat(1, ROI_fmb_all, ROI_fmb);
            f_baselines = [f_baselines; f_baseline];

            Xall = cat(1, Xall, Xc+X0(patch_no));
            Yall = cat(1, Yall, Yc+Y0(patch_no));
            Zall = cat(1, Zall, ones(length(Xc),1)*Z0(patch_no));
        end

        prepData.functional.Gdff = permute(ROI_dff_all, [2,1,3]);
        
        prepData.header.Nrois = Nrois;
        prepData.header.XYZ = [Xall, Yall, Zall];

    end
    
    %%
    if save_prep_data
        cd(path_save);

        save(prep_name, 'prepData');
    end

end
