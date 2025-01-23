
run '/Users/sadra/Desktop/MFB/Init/Initialize';

goodFolders__wheelMI;
            
% -
N_itr = 500; 
sig_lev = 3;
T_state = 3; %seconds

states_analysis = 1;
all_states_analysis = 0;

% for all
corr_whisking_analysis = 0;
corr_fl_analysis = 0;
pairwise_cc_analysis = 0;

% excluding after puff
corr_noPuff = 0;
corr_noPuff_locomotion = 1;
corr_noPuff_whisking = 1;
corr_noPuff_forelimb = 0;
corr_noPuff_pairwise = 1;

regress_analysis = 0;

NN_analysis = 0;

save_results = 1;

%folder_names = {'191209_13_44_12'};
%folder_names = {'191018_13_39_41'};

%  folder_names = {'191018_13_04_19', %with wheel MI
%                  '191018_13_56_55', %with wheel MI
%                  
%                  '191018_14_11_33', %with wheel MI
%                  '191018_14_30_00', %with wheel MI
%                  
% %                 %'191018_13_39_41',                
%                 };


folder_names = {
    %'191209_14_04_14'
% '191209_14_46_58',
%'191209_14_18_13',
% '191209_14_32_39',
'191209_15_01_22',
     };
 
 %'191018_13_56_55'

%%
for ijk = 1:length(folder_names)
    
    folder_name = folder_names{ijk}(1:15)
    prep_name = [path_home '/Data_processed/processedData_' folder_name '.mat'];
    load(prep_name);

    dt = prepData.header.dt; 
    T = prepData.header.T;
    Ntr = prepData.header.Ntr;
    
    tids_noPuff = logical((T<5e3) + (T>10e3));

    block_size = round(200/dt); %200 ms

    Tb_tr = prepData.speed.Tb_tr;
    Tf_tr = prepData.speed.Tf_tr;
    
    if isfield(prepData.functional, 'Gdff_MFs') 
        dff = prepData.functional.Gdff_MFs;
        xyz = prepData.header.XYZ_MF_clust;    
    else
        dff = prepData.functional.Gdff;
        xyz = prepData.header.XYZ; 
    end
    
    Nmf = size(dff,2);
    if Nmf>0; dff = dff(:,:,1:length(T)); end
    
    dff_r = reshape(permute(dff, [2,3,1]), Nmf, []);
    tr = (0:size(dff_r,2)-1)*dt/1e3;
     
    % as videos are not available after these times
    nnids = [];
    if strcmp(folder_name,'170712_21_31_26'); nnids = (tr>300); 
    elseif strcmp(folder_name,'170711_17_47_08'); nnids = (tr>280); 
    elseif strcmp(folder_name,'171212_16_19_37'); nnids = (tr>390); 
    end
    %dff_r(:,nnids) = nan;
    
    dff_rz = ((dff_r' - nanmean(dff_r') ) ./ nanstd(dff_r'))';
    
    spd = prepData.speed.spDall;
    spd_r = fastsmooth(reshape(spd', 1,[]),15,3,1);
    
    if isfield(prepData, 'states') 
        run_rr = prepData.states.run_rr;
    end
    
    %% regression analysis (locomotion and whisking coefficients)
    if regress_analysis
        
        reg_loc = [];
        reg_wsk = [];
        if ~isempty(prepData.MIs.whisking) && ~isempty(dff_r)
                           
             wsk = prepData.MIs.whisking(:,:);
             wsk_r = fastsmooth(reshape(wsk', 1,[]),15,3,1);
             b2 = (wsk_r - nanmean(wsk_r)) / nanstd(wsk_r); 
             
             run_rr = prepData.states.run_rr;
             b1 = (run_rr - nanmean(run_rr)) / nanstd(run_rr);
             
            for i = 1:Nmf
                
                y = dff_r(i,:);
                nnids = logical(~isnan(b1) .* ~isnan(b2) .* ~isnan(y));
                
                Y = y(nnids);
                X = [ones(1,length(Y)); b1(nnids); b2(nnids)];
                
                r = X'\Y';
                
                reg_loc = [reg_loc, r(2)];
                reg_wsk = [reg_wsk, r(3)];
                
            end
        
        end
        
        prepData.regression.reg_loc = reg_loc;
        prepData.regression.reg_wsk = reg_wsk;
            
    end
    
    %% states based on MI wheel
    if all_states_analysis
        
        if ~isempty(prepData.MIs.whisking) && ~isempty(dff)
            
            whl = prepData.MIs.wheel;
            whl_rs = fastsmooth(reshape(whl', 1,[]),T_state*1e3/dt,3,1);
            L_th = prctile(whl_rs(~isnan(run_rr)),10)*1.2;
            L_state = zeros(size(whl_rs));
            L_state(whl_rs>L_th) = 1;
            L_state(nnids) = nan;
            
            wsk = prepData.MIs.whisking;
            wsk_rs = fastsmooth(reshape(wsk', 1,[]),T_state*1e3/dt,3,1);
            W_th = prctile(wsk_rs(~isnan(run_rr)),10)*1.2;
            W_state = zeros(size(wsk_rs));
            W_state(wsk_rs>W_th) = 1;
            W_state(nnids) = nan;
            
            %
            wsk_r = fastsmooth(reshape(wsk', 1,[]),10,3,1);
            
            mids = logical((W_state==1));
            [cc_ww_W, zc_ww_W] = bootstrap_cc(dff_r(:,mids)', wsk_r(mids)', block_size, N_itr);
            
            mids = logical((W_state==1).*(L_state==0));
            [cc_ww_WL0, zc_ww_WL0] = bootstrap_cc(dff_r(:,mids)', wsk_r(mids)', block_size, N_itr);
            
            %
            prepData.all_states.W = W_state;
            prepData.all_states.L = L_state;
            
            prepData.all_states.cc_ww_W = cc_ww_W;
            prepData.all_states.cc_ww_WL0 = cc_ww_WL0;
            
            prepData.all_states.zc_ww_W = zc_ww_W;
            prepData.all_states.zc_ww_WL0 = zc_ww_WL0;
        end
         
    end
    
    %% states based on MI wheel
    if states_analysis
         
         whl = prepData.MIs.wheel;
         whl_r = fastsmooth(reshape(whl', 1,[]),15,3,1);
    
         whl_b = prctile(whl_r,10);
         whl_th = whl_b+.1;
         whl_s = fastsmooth(reshape(whl', 1,[]),100,3,1);         
         run_r = zeros(size(whl_s));
         run_r(whl_s <= whl_th) = 0;
         run_r(whl_s >= whl_th) = 1;
         %run_r(isnan(whl_r)) = nan;

         run_rr = run_r;
         ids_trans = [1, find(abs(diff(run_r)) == 1)];
         ids_long = ((diff(ids_trans)*dt/1e3) >= T_state);
         for ii = 1:length(ids_long)
             if ~ids_long(ii) 
                 %run_rr(ids_trans(ii):ids_trans(ii+1)) = nan;
                 run_rr(ids_trans(ii):ids_trans(ii+1)) = 0;
             end
         end 
         run_rr(nnids) = nan;
         
         if Nmf > 0
             
            % - modulation with locomotion
            [cc_wL, zc_wL] = bootstrap_cc(dff_r(:,~isnan(run_rr))', run_rr(~isnan(run_rr))', block_size, N_itr);

            run_rr_aligned = trial_resolved_trans(tr*1e3, run_rr, Tb_tr, Tf_tr, Ntr, T);

            prepData.states.run_rr = run_rr;
            prepData.states.run_rr_aligned = run_rr_aligned;
            prepData.states.T_state = T_state;
            prepData.states.whl_th = whl_th;
            %prepData.states.loc_mod = z_mod;

            prepData.CC.cc_wL = cc_wL;
            prepData.CC.zc_wL = zc_wL;
   
         end
    end
    
    %% correlations excluding after-puff
    if corr_noPuff
        
        if isfield(prepData, 'states') && ~isempty(dff)
        
            run_rr_aligned = prepData.states.run_rr_aligned;
            
            dff_r_noPuff = reshape(permute(dff(:,:,tids_noPuff), [2,3,1]), Nmf, []);
            run_rr_noPuff = reshape(run_rr_aligned(:,tids_noPuff)', 1, []);
            
            % - locomotion
            if corr_noPuff_locomotion
                
            [cc_wL_noPuff, zc_wL_noPuff] = bootstrap_cc(dff_r_noPuff(:,~isnan(run_rr_noPuff))', ...
                                          run_rr_noPuff(~isnan(run_rr_noPuff))', block_size, N_itr);
            
            prepData.CC.cc_wL_noPuff = cc_wL_noPuff;
            prepData.CC.zc_wL_noPuff = zc_wL_noPuff;
        
            end
            
            % - whisking 
            if corr_noPuff_whisking
                if ~isempty(prepData.MIs.whisking)
                    wsk = prepData.MIs.whisking(:,:);
                    wsk_r_noPuff = fastsmooth(reshape(wsk(:,tids_noPuff)', 1,[]),3,3,1);

                    [cc_ww_noPuff, zc_ww_noPuff] = bootstrap_cc(dff_r_noPuff(:,~isnan(run_rr_noPuff))', ...
                                                                wsk_r_noPuff(~isnan(run_rr_noPuff))', block_size, N_itr);

                    [cc_ww_Quiet_noPuff, zc_ww_Quiet_noPuff] = bootstrap_cc(dff_r_noPuff(:,run_rr_noPuff==0)', ...
                                                                wsk_r_noPuff(run_rr_noPuff==0)', block_size, N_itr);
                    
                    if sum(run_rr_noPuff==1)>0
                        [cc_ww_Loc_noPuff, zc_ww_Loc_noPuff] = ...
                        bootstrap_cc(dff_r_noPuff(:,run_rr_noPuff==1)', wsk_r_noPuff(run_rr_noPuff==1)', block_size, N_itr);
                    else
                     cc_ww_Loc_noPuff = nan*ones(size(cc_ww_noPuff));
                     zc_ww_Loc_noPuff = nan*ones(size(zc_ww_noPuff));
                    end 
             
                    prepData.CC.cc_ww_noPuff = cc_ww_noPuff;
                    prepData.CC.zc_ww_noPuff = zc_ww_noPuff;

                    prepData.CC.cc_ww_Quiet_noPuff = cc_ww_Quiet_noPuff;
                    prepData.CC.zc_ww_Quiet_noPuff = zc_ww_Quiet_noPuff;
                    
                    prepData.CC.cc_ww_Loc_noPuff = cc_ww_Loc_noPuff;
                    prepData.CC.zc_ww_Loc_noPuff = zc_ww_Loc_noPuff;
                end
            end
        
            % - pairwise
            if corr_noPuff_pairwise
                if Nmf>1
                   [cc_pw_noPuff, zc_pw_noPuff] = bootstrap_cc_pw(dff_r_noPuff(:,~isnan(run_rr_noPuff)), block_size, N_itr);
                   [cc_pw_Quiet_noPuff, zc_pw_Quiet_noPuff] = bootstrap_cc_pw(dff_r_noPuff(:,run_rr_noPuff==1), block_size, N_itr);
                end
                
                prepData.CC.cc_pw_noPuff = cc_pw_noPuff;
                prepData.CC.zc_pw_noPuff = zc_pw_noPuff;
                
                prepData.CC.cc_pw_Quiet_noPuff = cc_pw_Quiet_noPuff;
                prepData.CC.zc_pw_Quiet_noPuff = zc_pw_Quiet_noPuff;
       
%                prepData.CC.cc_pw_Q = cc_pw_Q;
%                prepData.CC.zc_pw_Q = zc_pw_Q;
%                prepData.CC.cc_pw_L = cc_pw_L;
%                prepData.CC.zc_pw_L = zc_pw_L;
% 
%                 prepData.CC.cc_ww_noPuff = cc_ww_noPuff;
%                 prepData.CC.zc_ww_noPuff = zc_ww_noPuff;
% 
%                 prepData.CC.cc_ww_Quiet_noPuff = cc_ww_Quiet_noPuff;
%                 prepData.CC.zc_ww_Quiet_noPuff = zc_ww_Quiet_noPuff;
            end
            
        end
    end
    
     %% corr with whisking
     if corr_whisking_analysis
         
        cc_ww = []; zc_ww = [];
        cc_ww_Q = []; zc_ww_Q = [];
        cc_ww_L = []; zc_ww_L = [];
        if ~isempty(prepData.MIs.whisking) && ~isempty(dff_r)
                           
             wsk = prepData.MIs.whisking(:,:);
             wsk_r = fastsmooth(reshape(wsk', 1,[]),15,3,1);
         
             [cc_ww, zc_ww] = bootstrap_cc(dff_r(:,~isnan(run_rr))', wsk_r(:,~isnan(run_rr))', block_size, N_itr);
             %cc_ww = corr(dff_r', wsk_r', 'rows', 'complete')';

             if sum(run_rr == 0)>0
                %cc_ww_Q = corr(dff_r(:,run_rr==0)', wsk_r(run_rr==0)', 'rows', 'complete')';
                [cc_ww_Q, zc_ww_Q] = bootstrap_cc(dff_r(:,run_rr==0)', wsk_r(run_rr==0)', block_size, N_itr);
             else 
                 cc_ww_Q = nan*ones(size(cc_ww));
                 zc_ww_Q = nan*ones(size(zc_ww));
             end

             if sum(run_rr==1)>0
                %cc_ww_L = corr(dff_r(:,run_rr==1)', wsk_r(run_rr==1)', 'rows', 'complete')';
                [cc_ww_L, zc_ww_L] = bootstrap_cc(dff_r(:,run_rr==1)', wsk_r(run_rr==1)', block_size, N_itr);
             else
                 cc_ww_L = nan*ones(size(cc_ww));
                 zc_ww_L = nan*ones(size(zc_ww));
             end 
        end
        
     prepData.CC.cc_whisk_all = cc_ww;
     prepData.CC.zc_whisk_all = zc_ww;

     prepData.CC.cc_whisk_quiet = cc_ww_Q;
     prepData.CC.zc_whisk_quiet = zc_ww_Q;

     prepData.CC.cc_whisk_locomotion = cc_ww_L;
     prepData.CC.zc_whisk_locomotion = zc_ww_L;
    end
    
   %% pairwise correlations
   if pairwise_cc_analysis
        
       cc_pw_Q = nan*ones([Nmf,Nmf]);
       cc_pw_L = nan*ones([Nmf,Nmf]);
       zc_pw_Q = nan*ones([Nmf,Nmf]);
       zc_pw_L = nan*ones([Nmf,Nmf]);
       
       if Nmf>1
           if nansum(run_rr==0) ~= 0
                [cc_pw_Q, zc_pw_Q] = bootstrap_cc_pw(dff_r(:,run_rr==0), block_size, N_itr);
           end
           if nansum(run_rr==1) ~= 0
                [cc_pw_L, zc_pw_L] = bootstrap_cc_pw(dff_r(:,run_rr==1), block_size, N_itr);
           end
       end
       
       prepData.CC.cc_pw_Q = cc_pw_Q;
       prepData.CC.zc_pw_Q = zc_pw_Q;
       prepData.CC.cc_pw_L = cc_pw_L;
       prepData.CC.zc_pw_L = zc_pw_L;
   end
   
   %% correlation with forelimb movement
   if corr_fl_analysis
   
      if isfield(prepData, 'DLC')
        flm = prepData.DLC.fl_mv;
        flm_r = reshape(flm', 1,[]);

        [cc_flm, zc_flm] = bootstrap_cc(dff_r(:,~isnan(run_rr))', flm_r(~isnan(run_rr))', block_size, N_itr);
        
        prepData.CC.cc_flm = cc_flm;
        prepData.CC.zc_flm = zc_flm;
      end
    
   end
   
   %%
   if NN_analysis
    
       how_many_NN = 4; 
       neib_dist = 40; %um
       
       for ij = 1:2
           cc = nan*ones([1,Nmf]);
           zc = nan*ones([1,Nmf]);
           if ij == 1
               if isfield(prepData.CC, 'cc_flm')
                   cc = prepData.CC.cc_flm;
                   zc = prepData.CC.zc_flm;
               end
           elseif ij == 2
               if isfield(prepData.CC, 'cc_wL')
                   cc = prepData.CC.cc_wL;
                   zc = prepData.CC.zc_wL;
               end
           end         
          
       N = size(xyz,1);
       dd = nan*ones(N,N);
       for i = 1:N
           dd(i,:) = sqrt(nansum((xyz(i,:) - xyz).^2,2)); 
           %dd(i,1:i)=nan;
           dd(i,i)=nan;
       end
       
       dd_sim_all = []; dd_diff_all = [];
       dd_sim_NN = []; dd_diff_NN = [];
       nn_ids_pall = {}; nn_ids_nall = {}; nn_ids_all = {};
       if Nmf > 4
           for i = 1:N
               [nn_dd,nn_ids0] = sort(dd(i,:));
               nn_ids0_p = nn_ids0(zc(nn_ids0)>sig_lev);
               nn_ids0_n = nn_ids0(zc(nn_ids0)<-sig_lev);
               
               nn_ids = nn_ids0(dd(i,nn_ids0)<neib_dist);  
               if length(nn_ids) > how_many_NN
                   nn_ids = nn_ids(1:how_many_NN);
               end
               nn_ids_p = nn_ids(zc(nn_ids)>sig_lev);
               nn_ids_n = nn_ids(zc(nn_ids)<-sig_lev);
               
               nn_ids_all{i} = nn_ids;
               nn_ids_pall{i} = nn_ids_p;
               nn_ids_nall{i} = nn_ids_n;
                   
               if zc(i) > sig_lev
                   dd_sim_all = [dd_sim_all, dd(i,nn_ids0_p)];
                   dd_diff_all = [dd_diff_all, dd(i,nn_ids0_n)];

                   dd_sim_NN = [dd_sim_NN, dd(i,nn_ids_p)];
                   dd_diff_NN = [dd_diff_NN, dd(i,nn_ids_n)];

               elseif zc(i) < -sig_lev
                   dd_sim_all = [dd_sim_all, dd(i,nn_ids0_n)];
                   dd_diff_all = [dd_diff_all, dd(i,nn_ids0_p)];

                   dd_sim_NN = [dd_sim_NN, dd(i,nn_ids_n)];
                   dd_diff_NN = [dd_diff_NN, dd(i,nn_ids_p)];
               end
           end
       end
       
       prepData.spatial.NN_ids = nn_ids_all;
       prepData.spatial.NN_ids_P = nn_ids_pall;
       prepData.spatial.NN_ids_N = nn_ids_nall;
       
       if ij == 1
       prepData.spatial.dd_sim_flm = dd_sim_all;
       prepData.spatial.dd_diff_flm = dd_diff_all;
       
       prepData.spatial.dd_sim_flm_NN = dd_sim_NN;
       prepData.spatial.dd_diff_flm_NN = dd_diff_NN;
       
       elseif ij == 2
       prepData.spatial.dd_sim_wL = dd_sim_all;
       prepData.spatial.dd_diff_wL = dd_diff_all;
       
       prepData.spatial.dd_sim_wL_NN = dd_sim_NN;
       prepData.spatial.dd_diff_wL_NN = dd_diff_NN;
       
       end
       
       end
   end
   
   %%
   
   if save_results
      save(prep_name, 'prepData');
   end
end

