
%% reading MI data
function MI_data = read_MIdata(file_name, Ntr, T, Tb_tr, Tf_tr, Tshift)
    if nargin == 5; Tshift = 0; end  
    if ~exist(file_name)
        %'warning: MI does not exist for this folder'
        MI_data = [];
        disp('MI does not exist for this folder')
    else
        load(file_name, 'MI');
        % time shift (laser on)
        my_t = MI(:,2)-Tshift;
        my_z = MI(:,1);
        MI_data = trial_resolved_trans(my_t, my_z, Tb_tr, Tf_tr, Ntr, T);
    end
end