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
   %'200130_14_02_12 FunctAcq';
   '200130_14_15_24 FunctAcq';
   '200130_14_29_30 FunctAcq';
    };


for ii=1:length(folder_names)

file = char(folder_names(ii));

load(['X:\MFB\Smooth_data\', file, '.mat']);

DFF = dff0_r_smooth;
zz0 = DFF;
zz0(isnan(zz0)) = 0;

Nmf0 = size(dff0_r_smooth,1);

ldr = zeros(Nmf0)*nan;
for i = 1:Nmf0
    i
    for j = 1:Nmf0
        if i >= j
            ldr(i,j) = get_var_ratio(zz0(i,:), zz0(j,:));
            ldr(j,i) = ldr(i,j);
        end
    end
end

save(['X:\MFB\LDR\', file, '.mat'], 'ldr');
end
    cc0 = corr(dff0_r_smooth', 'rows', 'complete');
    cc0(eye(length(cc0))==1)=nan;








