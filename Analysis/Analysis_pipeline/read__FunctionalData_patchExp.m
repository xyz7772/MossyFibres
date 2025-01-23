%% reading functional data 
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

ROI_f = ones(N_miniroi, Ntr, size(T,2))*nan;
for i = 1:Ntr
    if i == 1; Tf0 = (0:1:Ncycl-1)*Tcycl + Tr;
        else; Tf0 = (0:1:Ncycl-1)*Tcycl + Tr + TintCS(i-1);
    end

    for j = 1:N_miniroi
        ROI_f(j,i,:) = interp1(Tf0, squeeze(ROI_mean0r(:,i,j)), T, 'spline', nan);
    end
end

ROI_fmb = ROI_f;% - nanmean(background(:));
f_baseline = prctile(reshape(ROI_f, [N_miniroi, Ntr*length(T)]),10,2);

ROI_dff = (ROI_fmb - f_baseline)./f_baseline;
zz  = importdata(strcat(prefix,num2str(patch_no),'_ROIcentroids.txt')); 
zd = zz.data;
Xc = zd(1:N_miniroi,2) *pixel_size; 
Yc = zd(1:N_miniroi,3) *pixel_size;