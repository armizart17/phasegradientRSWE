% SCRIP TO GENERATE SWS MAPS OF HOMOGENEOUS REVERBERANT FIELDS
% Creation: 26/03/2024 (EMZ)

clc, clear all, close all;
addpath(genpath(pwd));
%%
%% GENERATE DATA SWS INCLUSION ONE FRAME
nWaves = 10e3; % number of waves
v_freq = [300, 350, 400, 450,500,550, 600, 650, 700]; % [Hz]
v_freq = [300, 400,500,600,700];
v_freq = 700;
nFields = 1;

% Resoluc 0.1mm/pixel
% lambda = cs/freq
% 2.5m/s y 400Hz -> 6.25mm 
% 2.5m/s y 700Hz -> 3.57mm 
% 2.5m/s y 400Hz -> 11mm 
% 2.5m/s y 400Hz -> 6.25mm 

% 700Hz - 2.5m/s 2l = 7.14mm = 71pix 4.5m/s 2l=12.8mm)129px

% 800Hz 2.5m/s 1l = 3.125mm = 33pix 4.5m/s 1l = 5.625mm = 57pix
% 900Hz 2.5m/s 1l = 2.77mm = 29pix 4.5m/s 1l = 5mm = 51pix
% 1000Hz 2.5m/s 1l = 2.5mmm = 27pix 4.5m/s 1l = 4.5mm = 47pix

v_freq = [800, 900, 1000];
window_v = 1*[57 51 47] + 0; %1 lamda pixels
window_v = 2*[57 51 47] + 1; %2 lamda pixels

nFields = 1;
cont = 1;
tic;
for freq = v_freq
    window = window_v(cont);
    path1 = 'G:\Mi unidad\PUCP_PosGraduate\Med_Imag\Proyecto_Imag_Medic\Codigos_R-SWE\CamposReverberantes';
    path2 = [path1,'\','Data', num2str(freq),'Hz-',num2str(nWaves),'ondas'];

    for field = 1:nFields
        name = ['R-FIELD_inc_',num2str(field),'.mat'];
        R_Field = load([path2, '\', name]);
        dinf.dx = min(diff(R_Field.x));
        dinf.dz = min(diff(R_Field.z));
        frame = R_Field.pv_complexZ(:,:,1); % number of frame
        
        frame = real(frame'); %transpose for Z (vertical-axial) X(horizontal-lateral)
        
        og_size = size(frame);
        mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');
    
        % FUNCION SWS MAP (Con ajuste de curva)
        [k_z,R_ax,k_x,R_lat,k,sws_matrix] = theoretical_fitting(mirror_frame,[window window],freq,dinf,og_size);
    
        % EMPAQUETAR RESULTADOS
        theoricalFitting.k_z = k_z;
        theoricalFitting.R_ax = R_ax;
        theoricalFitting.k_x = k_x;
        theoricalFitting.R_lat = R_lat;
        theoricalFitting.k = k;
        theoricalFitting.sws_matrix = sws_matrix;
    
    % Save
    mkdir([path2, '\CurveFitting'])
    save([path2, '\CurveFitting\SWS_re_inc_',num2str(field),'.mat'],'theoricalFitting');

    end
    cont = cont+1;
end
toc


%% GENERATE DATA SWS INCLUSION ONE FRAME WITH PHASE ESTIMATOR QR SOLVER

addpath('G:\Mi unidad\PUCP_PosGraduate\Med_Imag\Proyecto_Imag_Medic\Codigos_R-SWE\PhaseEstimatorFinalVersion');

% 500Hz 2.5m/s 1l = 5mm = 51pix 4.5m/s 1l = 9mm = 91pix
% 600Hz 2.5m/s 1l = 4.17mm = 41pix 4.5m/s 1l = 7.5mm = 75pix
% 700Hz 2.5m/s 1l = 3.57mm = 35pix 4.5m/s 1l = 6.42mm = 65pix

% 800Hz 2.5m/s 1l = 3.125mm = 33pix 4.5m/s 1l = 5.625mm = 57pix
% 900Hz 2.5m/s 1l = 2.77mm = 29pix 4.5m/s 1l = 5mm = 51pix
% 1000Hz 2.5m/s 1l = 2.5mmm = 27pix 4.5m/s 1l = 4.5mm = 47pix

nWaves = 10e3; % number of waves
v_freq = [500, 600, 700, 800, 900, 1000];
v_freq = [500];
nFields = 1;

window = 15; %11 pixels as described in paper
nFields = 1;

w_kernel = [window, window];
tic;
for freq = v_freq
%     path1 = 'G:\Mi unidad\PUCP_PosGraduate\Med_Imag\Proyecto_Imag_Medic\Codigos_R-SWE\CamposReverberantes';
    path1 = pwd;
    path2 = [path1,'\','Data', num2str(freq),'Hz-',num2str(nWaves),'ondas'];

    for field = 1:nFields
        name = ['R-FIELD_inc_',num2str(field),'.mat'];
        R_Field = load([path2, '\', name]);
        dinf.dx = min(diff(R_Field.x));
        dinf.dz = min(diff(R_Field.z));
        frame = R_Field.pv_complexZ(:,:,1); % number of frame
        
        %frame = (frame'); %transpose for Z (vertical-axial) X(horizontal-lateral)
        
        og_size = size(frame);
        mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');
    
        % FUNCION SWS MAP (Con linearizacion)
        [grad_z,grad_x,k,sws_matrix] = phase_estimator_QRVersion(mirror_frame, w_kernel, freq, dinf, og_size, 0.33);
    
        % EMPAQUETAR RESULTADOS
        
        pg_QR.grad_z = grad_z;
        pg_QR.grad_x = grad_x;
        pg_QR.grad_k = k;
        pg_QR.sws_matrix = sws_matrix;
       
        
    % Save
%     mkdir([path2, '\PhaseEstimator'])
%     save([path2, '\PhaseEstimator\SWS_PG_QR_inc_',num2str(field),'.mat'],'pg_QR');

    end
end
toc