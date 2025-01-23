
cd('/Users/sadra/Google Drive/Hana2017_Analysis/Data/forSadra /ORIGINAL_141120_16_41_00/Functional_Data');

Trial_no = 9; 
for ROI_no = [10, 13, 36, 156]
    z0 = load(strcat('ROI-',num2str(ROI_no,'%03i'),'_POI_Dwell-4.0us_ChB_Trial-',num2str(Trial_no,'%02i'),'.dat'));

    zb = median(z0);
    dff = (z0 - zb) / zb;

    [c,s] = deconvolveCa(dff, 'ar1');

    figure(); hold on;

    plot(dff);
    
    sp_ids = find(s>0);
    plot(sp_ids, ones(length(sp_ids)), 'marker','+');
end

