clc
clear all
close all force
addpath('core') 

DataFolder='results/'+string(datetime('now','TimeZone','local','Format','yyyy-MM-dd/HH-mm-ss'));
multiparams.DataFolder=DataFolder;
multiparams.n_simulations=10;
multiparams.dt=1/100;
multiparams.n_satelites=4;
multiparams.n_satelites_spoof=4;
multiparams.spoof_type='decalibrate';% decalibrate, trajectory, dangerzone, sigmaincrement
multiparams.spoof_technique='progressive';%none, delta, progressive, complete
multiparams.case='case5';%
multiparams.IMU='litton200';%

multiparams.dtGNSS=1.0;
[dataoutput10] = IMUcompare(multiparams);
 close all force
 
 multiparams.dtGNSS=0.1;
%[dataoutput1] = noincfun(multiparams);
[dataoutput1] = dtGNSSincfun(multiparams);
 close all force
 multiparams.spoof_type='trajectory'; % decalibrate, trajectory, dangerzone
 [dataoutput2] = dtGNSSincfun(multiparams);
 close all force
 multiparams.spoof_type='dangerzone'; % decalibrate, trajectory, dangerzone
 [dataoutput3] = dtGNSSincfun(multiparams);
 close all force
% % 


multiparams.n_simulations = 5;
multiparams.spoof_type='decalibrate'; % decalibrate, trajectory, dangerzone
multiparams.dtGNSS = 0.5;%Frecuencia de acualizacion de medidas GPS, entre dt y 50-60
multiparams.n_satelites = 3;%Numero de satelites totales
multiparams.n_satelites_spoof = 3;%Numero de satelites spoofeados
[dataoutput4] = nsatelitesincfun(multiparams);
close all force
 multiparams.spoof_type='trajectory'; % decalibrate, trajectory, dangerzone
[dataoutput5] = nsatelitesincfun(multiparams);
close all force
 multiparams.spoof_type='dangerzone'; % decalibrate, trajectory, dangerzone
[dataoutput6] = nsatelitesincfun(multiparams);
close all force
% 


multiparams.n_simulations = 7;
multiparams.dtGNSS = 1.0;%Frecuencia de acualizacion de medidas GPS, entre dt y 50-60
multiparams.n_satelites = 8;%Numero de satelites totales
multiparams.n_satelites_spoof = 1;%Numero de satelites spoofeados
multiparams.spoof_type='decalibrate'; % decalibrate, trajectory, dangerzone
[dataoutput7] = nsatelitesspoofincfun(multiparams);
close all force
multiparams.spoof_type='trajectory'; % decalibrate, trajectory, dangerzone
[dataoutput8] = nsatelitesspoofincfun(multiparams);
close all force
multiparams.spoof_type='dangerzone'; % decalibrate, trajectory, dangerzone
[dataoutput9] = nsatelitesspoofincfun(multiparams);
close all force



multiparams.n_simulations = 15;
multiparams.dtGNSS = 0.5;%Frecuencia de acualizacion de medidas GPS, entre dt y 50-60
multiparams.n_satelites = 6;%Numero de satelites totales
multiparams.dt=1/1000;
multiparams.n_satelites_spoof = 3;%Numero de satelites spoofeados
multiparams.spoof_type='decalibrate'; % decalibrate, trajectory, dangerzone
[dataoutput11] = dtincfun(multiparams);
close all force
multiparams.spoof_type='trajectory'; % decalibrate, trajectory, dangerzone
[dataoutput12] = dtincfun(multiparams);
close all force
multiparams.spoof_type='dangerzone'; % decalibrate, trajectory, dangerzone
[dataoutput13] = dtincfun(multiparams);
close all force

