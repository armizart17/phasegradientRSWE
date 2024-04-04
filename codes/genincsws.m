% SCRIP TO GENERATE SWS MAPS OF HOMOGENEOUS REVERBERANT FIELDS
% Creation: 26/03/2024 (EMZ)

clc, clear all, close all;
addpath(genpath(pwd));
%%
%% GENERATE DATA SWS INCLUSION ONE FRAME CURVE FITTING
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
%     path1 = 'G:\Mi unidad\PUCP_PosGraduate\Med_Imag\Proyecto_Imag_Medic\Codigos_R-SWE\CamposReverberantes';
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


%% GENERATE DATA SWS INCLUSION ONE FRAME WITH PHASE ESTIMATOR QR SOLVER v1.0

addpath(genpath(pwd));
fprintf('-------QR solver 1.0-------\n')

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
stride = 1;
nFields = 1;

w_kernel = [window, window];
tic;


pathdata = './dataold/';
pathout = './out/old/';

if ~exist("pathout","dir"); mkdir(pathout); end

for freq = v_freq
   
    pathfreq_in = [pathdata,'Data', num2str(freq),'Hz-',num2str(nWaves),'ondas/'];
    pathfreq_out = [pathout, 'Out', num2str(freq),'Hz/'];

    if ~exist(pathfreq_out,"dir"); mkdir(pathfreq_out); end

    for field = 1:nFields
        name = ['R-FIELD_inc_',num2str(field),'.mat'];
        R_Field = load([pathfreq_in, name]);
        dinf.dx = min(diff(R_Field.x));
        dinf.dz = min(diff(R_Field.z));
        frame = R_Field.pv_complexZ(:,:,1); % number of frame
        
        %frame = (frame'); %transpose for Z (vertical-axial) X(horizontal-lateral)
        
        og_size = size(frame);
        mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');
    
        % FUNCION SWS MAP (Con linearizacion)
        [grad_z,grad_x,k,sws_matrix] = phase_estimator_QR_kernel(mirror_frame, w_kernel, freq, dinf, og_size, 0.33, stride);
    
        % EMPAQUETAR RESULTADOS
        
        pg_QRv1.grad_z = grad_z;
        pg_QRv1.grad_x = grad_x;
        pg_QRv1.grad_k = k;
        pg_QRv1.sws_matrix = sws_matrix;
       
        
    % Save
    pathQR = [pathfreq_out, 'PhaseGradientQRv1/'];
    if ~exist(pathQR,"dir"); mkdir(pathQR); end
%     save([pathQR, 'SWS_PG_QRv1_homo_',num2str(field),'.mat'],'pg_QRv1');

    end
end
toc
fprintf('---------------------------\n')

%%
%% GENERATE DATA SWS INCLUSION ONE FRAME WITH PHASE ESTIMATOR QR SOLVER v2.0 BIG MATRIX

addpath(genpath(pwd));

fprintf('-------QR solver 2.0-------\n')

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
stride = 1;
tic;

pathdata = './dataold/';
pathout = './out/old';

if ~exist("pathout","dir"); mkdir(pathout); end

for freq = v_freq
   
    pathfreq_in = [pathdata,'Data', num2str(freq),'Hz-',num2str(nWaves),'ondas/'];
    pathfreq_out = [pathout, 'Out', num2str(freq),'Hz/'];

    if ~exist(pathfreq_out,"dir"); mkdir(pathfreq_out); end

    for field = 1:nFields
        name = ['R-FIELD_inc_',num2str(field),'.mat'];
        R_Field = load([pathfreq_in, name]);
        dinf.dx = min(diff(R_Field.x));
        dinf.dz = min(diff(R_Field.z));
        frame = R_Field.pv_complexZ(:,:,1); % number of frame
            
        %frame = (frame'); %transpose for Z (vertical-axial) X(horizontal-lateral)
        
        og_size = size(frame);
        mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');
    
        % FUNCION SWS MAP (Con linearizacion)
        [grad_z,grad_x,k,sws_matrix] = phase_estimator_QR_bigmat(mirror_frame, w_kernel, freq, dinf, og_size, 0.33, stride);
    
        % EMPAQUETAR RESULTADOS
        
        pg_QRv2.grad_z = grad_z;
        pg_QRv2.grad_x = grad_x;
        pg_QRv2.grad_k = k;
        pg_QRv2.sws_matrix = sws_matrix;
       
        
    % Save
    pathQR = [pathfreq_out, 'PhaseGradientQRv2/'];
    if ~exist(pathQR,"dir"); mkdir(pathQR); end
%     save([pathQR, 'SWS_PG_QRv2_homo_',num2str(field),'.mat'],'pg_QRv2');

    end
end
toc
fprintf('---------------------------\n')

%%




%%
tol = 1e-3;
mask = ones(M*N,1);
isotropic = true;
colMajor = false;


%% GENERATE DATA SWS INCLUSION ONE FRAME WITH PHASE ESTIMATOR WITH TNV (multi frequency)

addpath(genpath(pwd));

fprintf('-------QR solver 2.0-------\n')

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
stride = 1;
tic;

pathdata = './dataold/';
pathout = './out/old';

if ~exist("pathout","dir"); mkdir(pathout); end

for freq = v_freq
   
    pathfreq_in = [pathdata,'Data', num2str(freq),'Hz-',num2str(nWaves),'ondas/'];
    pathfreq_out = [pathout, 'Out', num2str(freq),'Hz/'];

    if ~exist(pathfreq_out,"dir"); mkdir(pathfreq_out); end

    for field = 1:nFields
        name = ['R-FIELD_inc_',num2str(field),'.mat'];
        R_Field = load([pathfreq_in, name]);
        dinf.dx = min(diff(R_Field.x));
        dinf.dz = min(diff(R_Field.z));
        frame = R_Field.pv_complexZ(:,:,1); % number of frame
            
        %frame = (frame'); %transpose for Z (vertical-axial) X(horizontal-lateral)
        
        og_size = size(frame);
        mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');
    
        % FUNCION SWS MAP (Con linearizacion)
        [grad_z,grad_x,k,sws_matrix] = phase_estimator_QR_bigmat(mirror_frame, w_kernel, freq, dinf, og_size, 0.33, stride);
    
        % EMPAQUETAR RESULTADOS
        
        pg_QRv2.grad_z = grad_z;
        pg_QRv2.grad_x = grad_x;
        pg_QRv2.grad_k = k;
        pg_QRv2.sws_matrix = sws_matrix;
       
        
    % Save
    pathQR = [pathfreq_out, 'PhaseGradientQRv2/'];
    if ~exist(pathQR,"dir"); mkdir(pathQR); end
%     save([pathQR, 'SWS_PG_QRv2_homo_',num2str(field),'.mat'],'pg_QRv2');

    end
end
toc
fprintf('---------------------------\n')


%% SIMPLE SWS PLOT
cm = 1e2;
mm = 1e3;
figure, 
sgtitle('\bf Results ESTIMATOR QRv2')
set(gcf, 'Position',[100 200 1500 500]);

 Create first subplot for Particle Velocity
subplot(121),
imagesc(R_Field.x*cm, R_Field.z*cm, real(R_Field.pv_complexZ(:,:,1)));
xlabel('Lateral'), ylabel('Axial'), title('Part. Veloc.');
axis equal; axis tight;
colormap(gca); colorbar% Apply pink colormap to the current axes

% Create second subplot for SWS
subplot(122),
imagesc(R_Field.x*cm, R_Field.z*cm, sws_matrix, [0 5]);
colormap(gca, jet); % Apply jet colormap to the current axes
colorbar;
axis equal; axis tight;
xlabel('Lateral'), ylabel('Axial'), 

uu = round(mean(sws_matrix(:)),3);
ss = round(std(sws_matrix(:)),3);
title(['SWS: ', num2str(uu), '\pm', num2str(ss) ]);


%% Plot comparing QRv1 vs QRv2

cm = 1e2;
mm = 1e3;
figure, 
sgtitle('\bf Comparison Estimator QR')
set(gcf, 'Position',[100 200 1500 500]);


% Create first subplot for Particle Velocity
subplot(131),
imagesc(R_Field.x*cm, R_Field.z*cm, real(R_Field.pv_complexZ(:,:,1)));
xlabel('Lateral'), ylabel('Axial'), title('Part. Veloc.');
axis equal; axis tight;
colormap(gca); colorbar% Apply pink colormap to the current axes


subplot(132),
imagesc(R_Field.x*cm, R_Field.z*cm, pg_QRv1.sws_matrix, [0 5]);
colormap(gca, jet); % Apply jet colormap to the current axes
colorbar;
axis equal; axis tight;
xlabel('Lateral'), ylabel('Axial'), 

uu = round(mean(pg_QRv1.sws_matrix(:)),3);
ss = round(std(pg_QRv1.sws_matrix(:)),3);
title(['QRv1 SWS: ', num2str(uu), '\pm', num2str(ss) ]);


pg_QRv2.sws_matrix_big = bigImg(pg_QRv2.sws_matrix, pg_QRv1.sws_matrix);
uu = round(mean(pg_QRv2.sws_matrix_big(:)),3);
ss = round(std(pg_QRv2.sws_matrix_big(:)),3);


subplot(133),
imagesc(R_Field.x*cm, R_Field.z*cm, pg_QRv2.sws_matrix_big, [0 5]);
colormap(gca, jet); % Apply jet colormap to the current axes
colorbar;
axis equal; axis tight;
xlabel('Lateral'), ylabel('Axial'), 

% uu = round(mean(pg_QRv2.sws_matrix(:)),3);
% ss = round(std(pg_QRv2.sws_matrix(:)),3);
title(['QRv2, st=',num2str(stride),' SWS: ', num2str(uu), '\pm', num2str(ss) ]);


%% GENERATE DATA SWS INCLUSION ONE FRAME WITH PHASE ESTIMATOR TIKHONOV v1.0 BIG MATRIX

addpath(genpath(pwd));

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
stride = 1;
tic;

% FOR REGULARIZATION parsTK
parsTK.lambda = 1e-7;
eps = 1e-100;
parsTK.version = 6;
parsTK.k = 1;
parsTK.beta = 1/(parsTK.lambda + eps); % 1e5;
parsTK.tolerance = 1e-3;
parsTK.operator = 'G';
parsTK.alpha = 2;
parsTK.maxIter = 100;

% parsTK.version = 11;
% parsTK.lambda = 1/1e8;

% parsTK.version = 3;
% parsTK.lambda = 0.000000001;
% parsTK.tolerance = 1e-3;

pathdata = './dataold/';
pathout = './out/old';

if ~exist("pathout","dir"); mkdir(pathout); end

for freq = v_freq
   
    pathfreq_in = [pathdata,'Data', num2str(freq),'Hz-',num2str(nWaves),'ondas/'];
    pathfreq_out = [pathout, 'Out', num2str(freq),'Hz/'];

    if ~exist(pathfreq_out,"dir"); mkdir(pathfreq_out); end

    for field = 1:nFields
        name = ['R-FIELD_inc_',num2str(field),'.mat'];

        R_Field = load([pathfreq_in, name]);
        dinf.dx = min(diff(R_Field.x));
        dinf.dz = min(diff(R_Field.z));
        frame = R_Field.pv_complexZ(:,:,1); % number of frame
        
        %frame = (frame'); %transpose for Z (vertical-axial) X(horizontal-lateral)
        
        og_size = size(frame);
        mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');
    
        % FUNCION SWS MAP (Con linearizacion)
%         [grad_z,grad_x,k,sws_matrix] = phase_estimator_QR_bigmat(mirror_frame, w_kernel, freq, dinf, og_size, 1.6, stride);
        [grad_z,grad_x,k,sws_matrix] = phase_estimator_TK_bigmat(mirror_frame, w_kernel, freq, dinf, og_size, 0.33, stride, parsTK);
    
        % EMPAQUETAR RESULTADOS
        
        pg_TKv1.grad_z = grad_z;
        pg_TKv1.grad_x = grad_x;
        pg_TKv1.grad_k = k;
        pg_TKv1.sws_matrix = sws_matrix;
       
        
    % Save
    pathTK = [pathfreq_out, 'PhaseGradientTikhov1/'];
    if ~exist(pathTK,"dir"); mkdir(pathTK); end
%     save([pathTK, 'SWS_PG_TKv1_homo_',num2str(field),'.mat'],'pg_TKv1');

    end
end
t = toc;
fprintf('Tikhonov version %d time : %f seconds\n', parsTK.version, t);


%% Plot comparing SWS MAPS QRv1(kernel) vs QRv2(bigmat) vs Tkv1(bigmat)

cm = 1e2;
mm = 1e3;
figure, 
sgtitle('\bf Comparison Estimators')
set(gcf, 'Position',[100 200 1500 500]);

% Create first subplot for Particle Velocity
subplot(141),
imagesc(R_Field.x*cm, R_Field.z*cm, real(R_Field.pv_complexZ(:,:,1)));
xlabel('Lateral'), ylabel('Axial'), title('Part. Veloc.');
axis equal; axis tight;
colormap(gca); colorbar% Apply pink colormap to the current axes


subplot(142),
imagesc(R_Field.x*cm, R_Field.z*cm, pg_QRv1.sws_matrix, [0 5]);
colormap(gca, jet); % Apply jet colormap to the current axes
colorbar;
axis equal; axis tight;
xlabel('Lateral'), ylabel('Axial'), 

uu = round(mean(pg_QRv1.sws_matrix(:)),3);
ss = round(std(pg_QRv1.sws_matrix(:)),3);
title(['QRv1 SWS: ', num2str(uu), '\pm', num2str(ss) ]);

pg_QRv2.sws_matrix_big = bigImg(pg_QRv2.sws_matrix, pg_QRv1.sws_matrix);
uu = round(mean(pg_QRv2.sws_matrix_big(:)),3);
ss = round(std(pg_QRv2.sws_matrix_big(:)),3);

subplot(143),
imagesc(R_Field.x*cm, R_Field.z*cm, pg_QRv2.sws_matrix_big, [0 5]);
colormap(gca, jet); % Apply jet colormap to the current axes
colorbar;
axis equal; axis tight;
xlabel('Lateral'), ylabel('Axial'), 

% uu = round(mean(pg_QRv2.sws_matrix(:)),3);
% ss = round(std(pg_QRv2.sws_matrix(:)),3);
title(['QRv2, st=',num2str(stride),' SWS: ', num2str(uu), '\pm', num2str(ss) ]);


pg_TKv1.sws_matrix_big = bigImg(pg_TKv1.sws_matrix, pg_QRv1.sws_matrix);
uu = round(mean(pg_TKv1.sws_matrix_big(:)),3);
ss = round(std(pg_TKv1.sws_matrix_big(:)),3);


subplot(144),
% imagesc(R_Field.x*cm, R_Field.z*cm, pg_TKv1.sws_matrix_big, [0 5]);
imagesc(R_Field.x*cm, R_Field.z*cm, pg_TKv1.sws_matrix_big, [0 5]);
colormap(gca, jet); % Apply jet colormap to the current axes
colorbar;
axis equal; axis tight;
xlabel('Lateral'), ylabel('Axial'), 

% uu = round(mean(pg_QRv2.sws_matrix(:)),3);
% ss = round(std(pg_QRv2.sws_matrix(:)),3);
title(['TKv1, st=',num2str(stride),' SWS: ', num2str(uu), '\pm', num2str(ss) ]);