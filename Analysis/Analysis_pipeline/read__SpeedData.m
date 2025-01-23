
%% reading speed data
cd(path_speed);

spddata = importdata('Speed data 001.txt');

% -- translating timestamps to real time
%if size(spddata.textdata, 2) == 2
    spd_timestamps = cell2mat(spddata.textdata(:,2));
% elseif size(spddata.textdata, 2) == 1
%     spd_timestamps = cell2mat(spddata.textdata(:,1));
%     spd_timestamps = spd_timestamps(:,12:end);
% end

%spdTreal0 = zeros([1,length(spd_timestamps)]);
spdTreal0 = (60^2 *str2num(spd_timestamps(:,1:2))+ ...
             60*str2num(spd_timestamps(:,4:5)) + ...
             str2num(spd_timestamps(:,7:8)) ) *1e3 + ...
             str2num(spd_timestamps(:,10:end))/1e3;
         
spdTreal = spdTreal0 - spdTreal0(1) + spddata.data(1,1)/1e3;

% -- finding the end / beginning points for trials
spdraw = spddata.data;
spdraw(:,1) = spdraw(:,1) * 1e-3; % in ms

spdff = diff(spdraw(:,1));
if spdff(length(spdff)) > 0; spdff(length(spdff)) = -1e3; end
sids = find(spdff < 0);
sids = vertcat(1, sids);

Ts = spdTreal(sids);
dTs = diff(Ts);
Tint = dTs(2:2:end);

sids_begin = sids(1:2:end) +1;
if length(sids_begin) > Ntr; sids_begin = sids_begin(1:Ntr); end
sids_end = sids(2:2:end);
if length(sids_end) > Ntr; sids_end = sids_end(1:Ntr); end
Ntr = length(sids_end);

Tb_tr = spdTreal(sids_begin);
%Tb_tr(1) = Tb_tr(1) - spdraw(1,1);
Tb_tr = Tb_tr - spdraw(sids_begin,1);

% fudge factor to account for the timing issue of the encoder 
%dt_encoder = mode(diff(spdraw(:,1)));
% zde = diff(spdraw(:,1)); 
% zde(zde<0) = nan;
% dt_encoder = nanmean(zde);
dt_encoder = 2.4; %
for i = 1:Ntr; Tb_tr(i,:) = Tb_tr(i,:) + dt_encoder*(i-1);
end

%Tf_tr = Tb_tr + spdTreal(sids_end);
Tf_tr = Tb_tr + Ttrial_real*1e3;

%Tint = Tb_tr(2:end) - Tf_tr(1:end-1);

TintCS = cumsum(Tint);

% -- interpolating speed to the universal time
speed = ones(Ntr, length(T))*nan;
for i = 1:Ntr
    i0 = sids_begin(i); %sids(2*i-1)+1;
    i1 = sids_end(i); %sids(2*i);
    %speed(i,:) = interp1((spdraw(i0:i1,1) - spdraw(i0,1)), spdraw(i0:i1,2), T, 'linear', 'extrap');
    speed(i,:) = interp1(spdraw(i0:i1,1), spdraw(i0:i1,2), T);
end
speed = -speed;
