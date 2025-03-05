function [peak_time, onset_time, rel_inc, maxima, minima] = peak_latency(x, tis_before, tis_after)

    % tpuff = 5e3;
    % time_wind_before = 0.5*1e3;
    % time_wind_after = 2*1e3;

    % Tb1 = tpuff-t_before; Tb2 = tpuff;
    % Ta1 = tpuff; Ta2 = tpuff+t_after;
    % tis_before = find( (T>=Tb1) .* (T<Tb2) ); 
    % tis_after = find( (T>=Ta1) .* (T<=Ta2) ); 
    
    %tt = T(tis_after)/1e3;
    std_th = 3;

    dt = 10;

    Ntot = size(x,1);

    peak_time = ones([1, Ntot])*nan;
    onset_time = ones([1, Ntot])*nan;
    rel_inc = ones([1,Ntot])*nan;
    maxima = ones([1,Ntot])*nan;
    minima = ones([1,Ntot])*nan;

    for ii = 1:Ntot    
        x_s = fastsmooth(x(ii,:), 10,3,1);

        zz_before = squeeze(x_s(tis_before));

        zz_after = squeeze(x_s(tis_after));

        [amx,bmx] = nanmax(zz_after);
        sig_inc = (amx > (nanmean(zz_before) + std_th*nanstd(zz_before)));
        
        [amn,bmn] = nanmin(zz_after);
        sig_dec = (amn < (nanmean(zz_before) - std_th*nanstd(zz_before)));

        if sig_inc
            peak_time(ii) = bmx*dt;
            
            sig_lev = find(zz_after > (nanmean(zz_before) + std_th*nanstd(zz_before)));
            if ~isempty(sig_lev)
                onset_time(ii) = sig_lev(1)*dt;
            end
            rel_inc(ii) = (amx - nanmean(zz_before)) / nanstd(zz_before);
            maxima(ii) = amx;
%             
%         if sig_dec
%             peak_time(ii) = bmn*dt;
% 
%             sig_lev = find(zz_after < (nanmean(zz_before) - std_th*nanstd(zz_before)));
%             if ~isempty(sig_lev)
%                 onset_time(ii) = sig_lev(1)*dt;
%             end     
%             rel_inc(ii) = (amn - nanmean(zz_before)) / nanstd(zz_before);
%             minima(ii) = amn;
% 
%         end


    end

end