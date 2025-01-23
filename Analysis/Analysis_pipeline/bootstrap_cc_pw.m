
% -- bootstrap analysis; pairwise cc
function [cc0, zc] = bootstrap_cc_pw(x, block_size, N_itr)
    
    cc0 = corr(x', 'rows', 'complete')';

    % - shuffle / cross validate
    Nmf = size(x,1);
    cc_sh = ones([N_itr, Nmf, Nmf]);
    ll1 = floor(size(x,2)/block_size) * block_size;
    x_sh = x;
    for itr = 1:N_itr
        for ish = 1:Nmf
            x_sh(ish,1:ll1) = randblock(x(ish,1:ll1), block_size);
        end
        cc_sh(itr, :,:) = corr(x_sh', 'rows', 'complete');
    end
    zc = (cc0 - squeeze(nanmean(cc_sh)) ) ./ squeeze(nanstd(cc_sh));
end
