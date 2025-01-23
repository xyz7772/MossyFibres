%% 
function [z_al] = trial_resolved_trans(my_t, my_z, Tb_tr, Tf_tr, Ntr, T)
    % align to trials
    z_al = nan*ones(Ntr, length(T));
    for i = 1:Ntr
        mids = find((my_t>=Tb_tr(i)) .* (my_t< Tf_tr(i)));
        if ~isempty(mids)   
            tt = my_t(mids); 
            zz = my_z(mids);
            [tx, index] = unique(tt);
            z_al(i,:) = interp1(tx-Tb_tr(i), zz(index), T);
        end
    end
end