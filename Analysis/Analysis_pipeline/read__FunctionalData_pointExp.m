cd(path_data);
load('allF_green.mat');

allFb = zeros(1, Nroi);
Gdff_raw  = zeros(Ntr, Nroi, Npt);
Gdff = zeros(Ntr, Nroi, size(T,2));
for j = 1:Nroi
    Fb = prctile(reshape(allF_green(:, j, :), [1, Ntr*Npt]),10);
    allFb(j) = Fb;
    Gdff_raw(:,j,:) = (allF_green(:,j,:) - Fb) / Fb;
    for i = 1:Ntr
        Tf0 = (0:1:Ncycl-1)*Tcycl + Troi(j);
        Gdff(i,j,:) = interp1(Tf0, squeeze(Gdff_raw(i,j,:)), T);
    end
end