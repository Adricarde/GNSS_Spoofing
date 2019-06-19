function [dataoutput] = nsatelitesspoofincfun(multiparams)

rootfolder=multiparams.DataFolder;
n_simulations=multiparams.n_simulations;
time_params.dt=multiparams.dt;
time_params.dtGNSS=multiparams.dtGNSS;
SAT.n_satelites = multiparams.n_satelites;
SAT.n_satelites_spoof =multiparams.n_satelites_spoof;
SAT.spoof_type=multiparams.spoof_type;
SAT.spoof_technique=multiparams.spoof_technique;
tray.case=multiparams.case;
IMU_params.IMUmodel=multiparams.IMU;


%PARAMETROS TIPICOS
time_params.alineamiento = 20000;%Tiempo de alineamiento inicial, entre 100 y 2000
%PARAMETROS TRAYECTORIA
tray.vehicle_speed=30;%Velocidad, el modelo no incorpora aceleraciones por ahora
tray.n_waypoints=7;%Numero de waypoints para hacer trayectorias aleatorias si se selecciona random mas abajo
tray.radio_min=50;%Radio Minimo
tray.height=10000;%Altitud del vehículo
%IMU
IMU_params.plot = 0;%Ploteo de IMU
IMU_params.biasmodel='complete';
%PLOTEOS
tray.plot = 0;%Si se quiere plotear solo la trayectoria y waypoints

%% ---------------- INICIO DEL MODELO ------------------------
tic
%INICIAMOS EL MODELO REAL
[WORLD] = tray_gen(tray,time_params);%Disponibles: 'circulo','worstcase','random'

%% -----------     SPOOFING     -----------------------
SAT.spoof_flag=1;
guess_vect(n_simulations,7)=0;detect_time(n_simulations,7)=0;
indicator(n_simulations,7)=0;spoof_time(n_simulations)=0;
false_alarm_time(n_simulations)=0;
deviated_distance(n_simulations)=0;
trajectories(n_simulations,11,size(WORLD.true_tray,1))=0;
dtGNSSvect(n_simulations)=0;
endtime(n_simulations)=0;

for i=1:n_simulations

    %SIMULATION
    %[IMUSIM] = IMU_sim(WORLD,IMU_params,time_params.dt,'litton200','complete');
    [IMUSIM] = IMU_sim(WORLD,IMU_params,time_params.dt);

    [SAT] = satelite_gen(WORLD,SAT,time_params);
    [Spoof_Sim] = propagation(WORLD,SAT,IMUSIM,time_params);
    Spoof_Sim.plot_postprocess=0;
    [guess_container] = postprocess(Spoof_Sim,time_params);
    pause(0.1);
    %DATA
    guess_vect(i,:)=guess_container.vect;
    detect_time(i,:)=guess_vect(i,:)-SAT.t_begin_spoof*time_params.dt;
    spoof_time(i)=SAT.t_begin_spoof*time_params.dt;
    false_alarm_time(i)=SAT.false_alarm_start*time_params.dt;
    trajectories(i,:,:)=Spoof_Sim.state_vector;
    
    XReal=WORLD.true_tray;
    deviated_distance(i)=abs(sqrt((Spoof_Sim.state_vector(1,end)-XReal(end,1))^2+(Spoof_Sim.state_vector(2,end)-XReal(end,2))^2));
    
    for j=1:7
        if detect_time(i,j) > 0 && guess_vect(i,j) ~= size(XReal,1)*time_params.dt
            indicator(i,j)=guess_vect(i,j)-SAT.t_begin_spoof*time_params.dt;
        else
            indicator(i,j)=0;
        end
    end
    endtime(i)=size(WORLD.true_tray,1);
    dtGNSSvect(i)=SAT.n_satelites_spoof;
    SAT.n_satelites_spoof=SAT.n_satelites_spoof+1;%Frecuencia de acualizacion de medidas GPS, entre dt y 50-60
    t_fin=toc;t_aux=(t_fin/i*n_simulations-t_fin)/60;
    fprintf('Tiempo restante estimado : %.3f \n\n', t_aux);
    fprintf('Numero de simulacion y total : %.0f / %.0f \n\n', i,n_simulations);
end
clc
t_fin=toc;fprintf('Tiempo de calculo : %.3f \n\n', t_fin);
%close(h)

dataoutput.guess_vect=guess_vect;
dataoutput.trajectories=trajectories;
dataoutput.indicator=indicator;
dataoutput.detect_time=detect_time;
dataoutput.spoof_time=spoof_time;
dataoutput.deviated_distance=deviated_distance;
dataoutput.false_alarm_time=false_alarm_time;
dataoutput.dtGNSSvect=dtGNSSvect;
dataoutput.endtime=endtime;

%%  --------------------   PLOTEO  -----------------------------------

savestring=strcat('Nsim_',string(n_simulations),...
    '_Nsat_',string(SAT.n_satelites),'_dtGNSS_',string(time_params.dtGNSS),...
    '_Stype_',string(SAT.spoof_type),'_Stechn_',string(SAT.spoof_technique),'_');

DataFolder=rootfolder+'/satspoofinc/';

if ~exist(char(DataFolder), 'dir')
    mkdir(char(DataFolder));
end

%----------------------------
plot_multisim(WORLD,dataoutput,DataFolder,savestring,time_params)


end


