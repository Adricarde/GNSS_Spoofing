function [IMU_SIMULATION] = IMU_sim(WORLDMODEL,IMU_params,dt)
%Esta funcion simula las medidas inerciales de un acelerometro y un
%giroscopo a lo largo de la trayectoria del vehiculo
%INPUTS
%bias_param: Se usan los parametros para modelizar el bias en el tiempo
%dt: Paso de tiempo
%XReal: Trayectoria, velocidad, aceleracion, orientacion y velocidad
%angular verdaderas a simular
%ABody: Aceleraciones y velocidades angulares en ejes cuerpo
%bias_type: Tipos de bias (constante o lineal)
%OUTPUTS
%AM: Medidas de los acelerometros en x e y, y velocidad angular z
%biasx, biasy, biasw: Los valores de tres bias en el tiempo
XReal=WORLDMODEL.true_tray;
ABody=WORLDMODEL.forces;
T=size(ABody,1);
%Parametros de los bias
%PARAMETROS BIAS
switch IMU_params.IMUmodel
    case 'default'
        IMU_params.sigma_bias_a = 0.5*9.81/1000;
        IMU_params.sigma_bias_w = 5*pi/180/3600;
        tau_lpf=100;
        tau_lpfw=300;
        RWa=0.3;
        RWw=0.1;
    case 'litton200'
        IMU_params.sigma_bias_a = 0.3*9.81/1000;
        IMU_params.sigma_bias_w = 0.5*pi/180/3600;
        tau_lpf=100;
        tau_lpfw=300;
        RWa=0.02;
        RWw=0.05;
    case 'MPU6500'
        IMU_params.sigma_bias_a = 60*9.81/1000;
        IMU_params.sigma_bias_w = 18000*pi/180/3600;
        tau_lpf=100;
        tau_lpfw=300;
        RWa=0.15;
        RWw=0.6;
    case 'HG4930'
        IMU_params.sigma_bias_a = 2*9.81/1000;
        %IMU_params.sigma_bias_w = 1*9.81/1000;
        IMU_params.sigma_bias_w = 10*pi/180/3600;
        tau_lpf=100;
        tau_lpfw=2000;
        %tau_lpfw=100;
        RWa=0.04;
        RWw=0.05;
        %RWw=0.6;
    case 'old'
        IMU_params.sigma_bias_a=0.0001;
        IMU_params.sigma_bias_w=0.0001;
        tau_lpf=100;
        tau_lpfw=300;
        RWa=0.3;
        RWw=0.1;
end

sigma_bias_a=IMU_params.sigma_bias_a; %Error del bias, ruido blanco
sigma_bias_w=IMU_params.sigma_bias_w; %Error del bias, ruido blanco

%Offset del Bias
bias0x=sigma_bias_a*5*(2*rand-1);
bias0y=sigma_bias_a*5*(2*rand-1);
bias0w=sigma_bias_w*5*(2*rand-1);

%Diferentes modelos de bias
switch IMU_params.biasmodel
    case 'constante'
        bias_driftx = 0;
        bias_drifty = 0;
        bias_driftw = 0;
    case { 'lineal','driftrand'}
        bias_driftx = bias0x/100;
        bias_drifty = bias0y/100;
        bias_driftw = bias0w/100;
end
fprintf('********************************************************\n');
fprintf('IMU Measures Generation...\n');
D=0;D2=0;D1=0;
AM(T,9)=0;drift_counter=1;
biasx(T)=0;biasy(T)=0;biasw(T)=0;
for i=1:T
    switch IMU_params.biasmodel
        case {'constante' , 'lineal'}
            %Formula de los diferentes bias, b=b0+b1*t+ruido
            biasx(i)=bias0x+bias_driftx*(drift_counter-1)*dt+sigma_bias_a*(2*rand-1);
            biasy(i)=bias0y+bias_drifty*(drift_counter-1)*dt+sigma_bias_a*(2*rand-1);
            biasw(i)=(bias0w+bias_driftw*(drift_counter-1)*dt+sigma_bias_w*(2*rand-1));
            drift_counter=drift_counter+1;
            if drift_counter > T/20
                drift_counter=1;
            end
        case 'complete'
            RW=RWa/60*(2*rand-1)*sqrt(2/dt/tau_lpf);
            RW1=RWa/60*(2*rand-1)*sqrt(2/dt/tau_lpf);
            D=dt*(-D/tau_lpf+sigma_bias_a*(2*rand-1)*sqrt(1/dt))+D;
            D1=dt*(-D1/tau_lpf+sigma_bias_a*(2*rand-1)*sqrt(1/dt))+D1;
            biasx(i)=bias0x+RW+D;
            biasy(i)=bias0y+RW1+D1;
            
            RW2=RWw/60*(2*rand-1)*sqrt(2/dt/tau_lpfw);
            D2=dt*(-D2/tau_lpfw+sigma_bias_w*(2*rand-1)*sqrt(1/dt))+D2;
            biasw(i)=bias0w+RW2+D2;
    end
    %Aceleraciones x e y, velocidad angular w, EJES CUERPO
    AM(i,1)=ABody(i,1)+biasx(i);
    AM(i,2)=ABody(i,2)+biasy(i);
    AM(i,3)=ABody(i,3)+biasw(i);
    AM(i,4)=XReal(i,3);
    AM(i,5)=XReal(i,4);
    AM(i,6)=XReal(i,5);
    AM(i,7)=XReal(i,6);
    AM(i,8)=XReal(i,8);
    AM(i,9)=XReal(i,7);
end
fprintf('Done\n');

IMU_SIMULATION.measures=AM;
IMU_SIMULATION.biasx=biasx;
IMU_SIMULATION.biasy=biasy;
IMU_SIMULATION.biasw=biasw;
IMU_SIMULATION.params=IMU_params;
IMU_SIMULATION.discarted_measures=0;
IMU_SIMULATION.checked_measures=zeros(T,1);

if IMU_params.plot
    figure(300)
    subplot(3,2,1),
    hold on
    plot(linspace(1,T,T)*dt,biasy),grid on, axis tight;
    xlabel('Time')
    ylabel('m/s^2')
    title('b_y')
    set(gca, 'FontName', 'Cambria','Fontsize',20)
    subplot(3,2,2),
    hold on
    plot(linspace(1,T,T)*dt,AM(:,2)),grid on, axis tight;
    xlabel('Time')
    ylabel('m/s^2')
    title('Yb Acceleration')
    set(gca, 'FontName', 'Cambria','Fontsize',20)
    subplot(3,2,3),
    hold on
    plot(linspace(1,T,T)*dt,biasx),grid on, axis tight;
    xlabel('Time')
    ylabel('m/s^2')
    title('b_x')
    set(gca, 'FontName', 'Cambria','Fontsize',20)
    subplot(3,2,4),
    hold on
    plot(linspace(1,T,T)*dt,AM(:,1)),grid on, axis tight;
    xlabel('Time')
    ylabel('m/s^2')
    title('Xb Acceleration')
    set(gca, 'FontName', 'Cambria','Fontsize',20)
    subplot(3,2,5),
    hold on
    plot(linspace(1,T,T)*dt,biasw),grid on, axis tight;
    xlabel('Time')
    ylabel('rad/s')
    title('b_w')
    set(gca, 'FontName', 'Cambria','Fontsize',20)
    subplot(3,2,6),
    hold on
    plot(linspace(1,T,T)*dt,AM(:,3)),grid on, axis tight;
    xlabel('Time')
    ylabel('rad/s')
    title('Zb Angular Speed')
    set(gca, 'FontName', 'Cambria','Fontsize',20)
    pause(0.01)
end
end


