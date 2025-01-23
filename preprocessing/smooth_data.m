%  Select
folder_names = {
   '171212_16_19_37';
   '191018_13_39_41';
   '191018_13_56_55';
   '191018_14_30_00';
   '191018_14_11_33';
   '191209_13_44_12';
   '191209_14_04_14';
   '191209_14_32_39';
   '191209_15_01_22';
   '191209_14_18_13';
   '191209_14_46_58';
   '200130_13_21_13 FunctAcq';
   '200130_13_36_14 FunctAcq';
   '200130_13_49_09 FunctAcq';
   '200130_14_02_12 FunctAcq';
   '200130_14_15_24 FunctAcq';
   '200130_14_29_30 FunctAcq';
    };

for ii=1:length(folder_names)
    file=char(folder_names(ii));
    quickAnalysis_Yizhou_fast; % load data
    dt = 10;
    Nmf0 = size(dff0,1);
    dff0_r = reshape(permute(dff0, [1,3,2]), Nmf0, []);
    smooth_window = 200/dt;
    dff0_r_smooth = zeros(size(dff0_r));
    for i = 1:Nmf0
        i
        dff0_r_smooth(i,:) = fastsmooth(dff0_r(i,:), smooth_window, 3, 1);
    end
    
    save(['X:\MFB\Smooth_data\', file, '.mat'], 'dff0_r_smooth');
end