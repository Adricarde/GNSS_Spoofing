clc
clear all
close all
addpath('core') 

%PARAMETROS TIPICOS
time_params.dt = 1/100;%Paso de tiempo, no disminuir mucho ni aumentar mas de 20-30 
time_params.dtGNSS = 0.5;%Frecuencia de acualizacion de medidas GPS, entre dt y 50-60
time_params.alineamiento = 25000;%Tiempo de alineamiento inicial, entre 100 y 2000
time_params.alineamiento = floor(150/time_params.dt);%Tiempo de alineamiento inicial, entre 100 y 2000
%PARAMETROS TRAYECTORIA
tray.vehicle_speed = 30;%Velocidad
tray.n_waypoints = 6;%Numero de waypoints para hacer trayectorias aleatorias si se selecciona random mas abajo
tray.plot = 0;%Si se quiere plotear solo la trayectoria y waypoints
tray.case='case3';

%IMU
IMU_params.plot = 0;%Ploteo de IMU
IMU_params.IMUmodel='litton200';% litton200, HG4930, MPU6500, old, default
IMU_params.biasmodel='complete';
%PARAMETROS SATELITE Y SPOOFING
SATSIM.n_satelites = 4;%Numero de satelites totales
SATSIM.n_satelites_spoof = 4;%Numero de satelites spoofeados
SATSIM.spoof_technique = 'progressive';%none, delta, progressive, stealth_progressive, hard_progressive, complete
SATSIM.spoof_type = 'decalibrate'; % decalibrate, trajectory, dangerzone

%% ---------------- INICIO DEL MODELO ------------------------
tic
%INICIAMOS EL MODELO REAL
[WORLD] = tray_gen(tray,time_params);%Disponibles: 'worstcase','random'

%GENERACION DE MEDIDAS IMU
[IMUSIM] = IMU_sim(WORLD,IMU_params,time_params.dt);
%[IMUSIM] = IMU_3merge(WORLD,IMU_params,time_params.dt);

%GENERACION DE MEDIDAS GNSS
[SATSIM] = satelite_gen(WORLD,SATSIM,time_params);

%% --------- SIN SPOOFING  ------------
%CONDICIONES INICIALES E INTEGRACION
SATSIM.spoof_flag=0;
[Auth_Sim] = propagation(WORLD,SATSIM,IMUSIM,time_params);
%% -----------     SPOOFING     -----------------------
% INTEGRACION
SATSIM.spoof_flag=1;
[Spoof_Sim] = propagation(WORLD,SATSIM,IMUSIM,time_params);
%% -----------  DETECCION DE SPOOFING  -----------------------
Spoof_Sim.plot_postprocess=1;
Auth_Sim.plot_postprocess=1;
[guess_container1] = postprocess(Auth_Sim,time_params);
[guess_container2] = postprocess(Spoof_Sim,time_params);

%% -------   PLOTEO BIAS AUTÉNTICO Y SPOOFEADO PARA VER SI HAY DIFERENCIAS
%Kalman_analysis(Auth_Sim,time_params);
plot_sim(WORLD,time_params,Auth_Sim,Spoof_Sim)
%plot_motion(time_params,WORLD,SATSIM,Auth_Sim,Spoof_Sim)
%statistical_analysis(Spoof_Sim,time_params)

time_params.t_fin=toc;fprintf('Tiempo de calculo : %.3f \n\n', time_params.t_fin);
