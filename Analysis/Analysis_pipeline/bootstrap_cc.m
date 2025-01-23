
% -- bootstrap analysis; cc with behvaiour
function [cc0, zc] = bootstrap_cc(x1, x2, block_size, N_itr)
    
    cc0 = corr(x1, x2, 'rows', 'complete')';

    % - shuffle / cross validate
    cc_sh = ones([N_itr,size(x1,2)]);
    ll1 = floor(length(x2)/block_size) * block_size;
    x2_sh = x2;
    for itr = 1:N_itr
        x2_sh(1:ll1) = randblock(x2(1:ll1), block_size);
        %x2_sh = randblock(x2(1:ll1), block_size);
        %x2_sh = [x2_sh; x2(ll1+1:end)];
        cc_sh(itr,:) = corr(x1, x2_sh, 'rows', 'complete');
    end
    zc = (cc0 - nanmean(cc_sh) ) ./ nanstd(cc_sh);
end