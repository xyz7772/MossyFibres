
%% reading header info
cd(path_header);

headerInfo = ini2struct('Experiment Header.ini');

Ntr = str2num(headerInfo.globalParameters.numberoftrials);
Nroi = str2num(headerInfo.globalParameters.numberofpoi);
Npt = str2num(headerInfo.globalParameters.numberofcycles);
Ttrial0 = str2num(headerInfo.globalParameters.acquisitionlength0x28sec0x29);
Ncycl = str2num(headerInfo.globalParameters.numberofcycles);
DwellT = str2num(headerInfo.globalParameters.dwelltime0x28us0x29);

fov = str2num(headerInfo.globalParameters.fieldofview);
frame_size = str2num(headerInfo.globalParameters.framesize);

pixel_size = fov / frame_size;

SCrt = load('Single cycle relative times.txt');
Troi = SCrt(:,1) /1e3;
% Tcycl_before = SCrt(1,2)/1e3; %- that was wrong! was it?
% Ttrial_real_before = Ncycl*Tcycl_before /1e3;
Tcycl = SCrt(end,1)/1e3 / (length(Troi)/Nroi); 
Ttrial_real = Ncycl*Tcycl /1e3;

MC_on = headerInfo.movementCorrection.movcorenabled0x3f;
if MC_on(1:4) == 'TRUE'; MC_on = 1;
else; MC_on = 0;
end

% universal time reference frame
dt = 10; %(ms)
T = 0:dt:round(Ttrial_real*1e3);

cd(path_header);

NIn = importdata('Normalised_ROI.dat');
Xn = NIn.data(:,5); 
Yn = NIn.data(:,6); 
Zn = NIn.data(:,7);

NI = importdata('ROI.dat');

X0 = NI.data(:,5) *pixel_size; 
Y0 = NI.data(:,6) *pixel_size; 
Z0 = NI.data(:,7) ;%*pixel_size;
Xf = NI.data(:,8) *pixel_size;
Yf = NI.data(:,9) *pixel_size; 
Zf = NI.data(:,10);%*pixel_size;

zdata = importdata('Zplane_Pockels_Values.dat');
zz_norm = zdata.data(:,2);
Z_aol = zdata.data(:,1);

