
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

rH = str2num(headerInfo.globalParameters.rectangleroiheight0x28pixels0x29);
rL = str2num(headerInfo.globalParameters.rectangleroilength0x28pixels0x29);

%--
zdata = importdata('Zplane_Pockels_Values.dat');
AOL_Znorm = zdata.data(:,2);
AOL_Z0 = zdata.data(:,1);

%--
SCrt = load('Single cycle relative times.txt');
Troi = SCrt(:,1)/1e3; 
%Tcycl = SCrt(1,2)/1e3; - that was wrong?
Tcycl = SCrt(end,1)/1e3; 
Ttrial_real = Ncycl*Tcycl/1e3;

dt = 10; %(ms)
T = 0:dt:round(Ttrial_real*1e3);

%--
%NI = importdata('Normalized_ROI.dat');
NI = importdata('ROI.dat');
N_lines = NI.data(1,13);

X0 = NI.data(:,5) *pixel_size; 
Y0 = NI.data(:,6) *pixel_size; 
Z0 = NI.data(:,7); 
Xf = NI.data(:,8) *pixel_size;
Yf = NI.data(:,9) *pixel_size; 
Zf = NI.data(:,10); 

Npatches = length(X0);

dX = Xf(1)-X0(1); dY = Yf(1)-Y0(1);

