
folder_names = {

% - Nigel
'171212_16_19_37', % Nigel patch

};

for ij = 1:length(folder_names)
    clearvars -except folder_names ij
    
    prefix = 'clean_patch';

    folder_name = folder_names{ij}
    patch_numbers = [1,2, 4:12, 14:15];
    
    path_data = ['/Volumes/Data/Data_raw/', folder_name];
    path_functional = [path_data, '/functional_patch'];
    path_speed = [path_data '/Speed_Data'];
    path_save = path_data;

    % - reading header info
    cd(path_data);
    headerInfo = ini2struct('Experiment Header.ini');
    Ntr = str2num(headerInfo.globalParameters.numberoftrials);
    Nroi = str2num(headerInfo.globalParameters.numberofpoi);
    Npt = str2num(headerInfo.globalParameters.numberofcycles);

    % - reading functional data 
    if exist(path_functional)
        cd(path_functional);
                
        zz  = importdata(strcat(prefix,num2str(patch_no),'_ROImeanArea.txt')); 

        zd = zz.data;
        ROI_area = zd(1,2:2:end-2); 
        N_miniroi = length(ROI_area); 

        ROI_mean0 = zd(:,3:2:end-2);
        ROI_mean0(bogus_frames,:) = nan;
        ROI_mean0r = reshape(ROI_mean0, [Ncycl, Ntr, N_miniroi]);

        background = zd(:,end);

        Tr = nanmean(Troi((patch_no-1)*N_lines+1: patch_no*N_lines));

        allF_green = ones(Ntr, N_miniroi, size(T,2))*nan;
        
        for i = 1:Ntr
            if i == 1; Tf0 = (0:1:Ncycl-1)*Tcycl + Tr;
            else; Tf0 = (0:1:Ncycl-1)*Tcycl + Tr + TintCS(i-1);
            end

            for j = 1:N_miniroi
                allF_green(i,j,:) = interp1(Tf0, squeeze(ROI_mean0r(:,i,j)), T, 'spline', nan);
            end
        end

        cd(path_save);
        save('allF_green.mat', 'allF_green');
    end
end