% main_script Diabetic Foot by RSWE
clear;
format shortg; clock
name_path = '\\NASLIM\PieDiabetico\Pacientes\Adquiridos\8760';   %where data is located
cd(name_path)
Read_files=dir;
addpath('C:\Users\u_imagenes\Documents\MATLAB\Codigo RSWE-20200222T180013Z-001\Codigo RSWE');
addpath('C:\Users\u_imagenes\Documents\MATLAB\Codigo RSWE-20200222T180013Z-001\Codigo RSWE\imoverlay')
set(0,'DefaultFigureColormap',jet);
freq = 400:50:600;
% Load data patient: mm- index of the file in the folder 
for mm=1:18% Processing mm=1:18//Verification mm=1
cd(name_path)
load(Read_files(mm+2).name)
path1 = [name_path,'\File_',num2str(mm)];
mkdir(path1) 
cd(path1)
offset_y = 0;
RFsamples = 2048;
dinf.fmin = 80; 
dinf.fmax = 500;
dinf.fc = Trans.frequency*1e6;
dinf.c0 = 1540;
dinf.wl = dinf.c0/dinf.fc;
dinf.dz = PData(2).PDelta(3)*dinf.wl;
dinf.dx = PData(2).PDelta(1)*dinf.wl;
dinf.samplesRF = RFsamples;
dinf.PRFe = 10*dinf.fmax;
dinf.Tprf = 1/dinf.PRFe;
dinf.Troundtrip = (dinf.samplesRF/2)*(1/dinf.fc);
dinf.maxDepth = dinf.Troundtrip*dinf.c0/2;
dinf.fs = 2*dinf.fc;
dinf.offset_y = offset_y; 
offset_x = dinf.dx*size(IQData,2)/2;
dinf.offset_x = offset_x;
dinf.Bmodedz = PData(1).PDelta(3)*dinf.wl;
dinf.Bmodedx = PData(1).PDelta(1)*dinf.wl;
dinf.num_angles = 4;
dinf.offset_z = dinf.dz*size(IQData,1)/2;
[u,dinf] = pv_cal(IQData,dinf,dinf.num_angles);% particle velocity estimates using Loupas'method
save('pv_cal.mat','u','dinf','freq','IQBmodeData');
end
clock