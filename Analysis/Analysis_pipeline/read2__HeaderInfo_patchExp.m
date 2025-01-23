
%% reading header info
cd(path_header);

%headerInfo = ini2struct('Experiment Header.ini');

File = 'Experiment Header.ini';
headerInfo = INI('File',File);
headerInfo.read();

Ntr = headerInfo.FUNCTIONAL_IMAGING.Number_of_trials;
Nroi = headerInfo.FUNCTIONAL_IMAGING.Number_of_pois;
Ncycl = headerInfo.FUNCTIONAL_IMAGING.Number_of_cycles;
Ttrial0 = headerInfo.FUNCTIONAL_IMAGING.Acquisition_Length;
%DwellT = str2num(headerInfo.FUNCTIONAL_IMAGING.Dw);

fov = headerInfo.FUNCTIONAL_IMAGING.field_of_view;
frame_size = headerInfo.FUNCTIONAL_IMAGING.Frame_Size;

pixel_size = fov / frame_size;

rH = headerInfo.FUNCTIONAL_IMAGING.Rectangle_ROI_height;
rL = headerInfo.FUNCTIONAL_IMAGING.Rectangle_ROI_length;

%--
% zdata = importdata('Zplane_Pockels_Values.dat');
% AOL_Znorm = zdata.data(:,2);
% AOL_Z0 = zdata.data(:,1);

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

%--
SCrt = importdata('Single cycle relative times_HW.txt');
% Troi = SCrt(:,1)/1e3; 
% %Tcycl = SCrt(1,2)/1e3; - that was wrong?
% Tcycl = SCrt(end,1)/1e3; 
% Ttrial_real = Ncycl*Tcycl/1e3;

% if struct?
if isstruct(SCrt)
    SCrt = SCrt.data;
else
end
% -
zdt = find(diff(SCrt(:,1))<0);

ijk = 5;
my_scrt = SCrt(zdt(ijk)+1:zdt(ijk+1),1);

Ttrial_real = (my_scrt(end)-my_scrt(1))/1e3;
Tcycl =  Ttrial_real/ Ncycl;

Troi = my_scrt(1:Npatches*N_lines)/1e3;
%Troi = Troi_all(1:N_lines:end);

% -

dt = 10; %(ms)
T = 0:dt:round(Ttrial_real);
