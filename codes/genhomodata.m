% SCRIP TO GENERATE HOMOGENEOUS REVERBERANT FIELDS
% Creation: 26/03/2024 (EMZ)

clc, clear all, close all;
addpath(genpath(pwd));

%% R-FIELD CHARACTERISTICS %%
nWaves = 10e3; % number of waves
v_freq = [300,350,400,450,500,550,600,650,700]; % [Hz]
v_freq = [800,900, 1000];

v_freq = [500]; % first test 500Hz
c_back = 2.5; % background SWS [m/s] 
c_inc = 4.5; % inclusion SWS [m/s] NOT APPLICABLE HERE
nFields = 1; % number of fields to generate at each freq
% vid_duration = 25e-3; % time of field generation [s]

%% ITERATION HOMO
pathData = './data/';
if ~exist("pathData","dir"); mkdir(pathData); end


for freq = v_freq
    for field = 1:nFields
        tic;

        
        path2 = [pathData, 'Data', num2str(freq),'Hz-',num2str(nWaves),'wvs'];
        mkdir(path2); % create path
        cd(path2)
        
        [pv_complexZ,x,z] = simulate_reverberant_complexZ_vhomo(nWaves, freq, c_back);
        
        name = ['R-FIELD_homo_',num2str(field),'.mat'];
        save(name,'pv_complexZ','x','z','freq','c_back');
        
        toc
    end
end
cd('../../'); % home directory

%% ITERATION INC UP 7.5mm

pathData = './data/';
if ~exist("pathData","dir"); mkdir(pathData); end


for freq = v_freq
    for field = 1:nFields
        tic;

        
        path2 = [pathData, 'Data', num2str(freq),'Hz-',num2str(nWaves),'wvs'];
        mkdir(path2); % create path
        cd(path2)
        
        [pv_complexZ,x,z] = simulate_reverberant_complexZ_vIncup(nWaves, freq, c_back, c_inc);
        
        name = ['R-FIELD_incup_',num2str(field),'.mat'];
        save(name,'pv_complexZ','x','z','freq','c_back');
        
        toc
    end
end
cd('../../'); % home directory