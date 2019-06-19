function [SIMULATION] = propagation(WORLDMODEL,SAT,IMU_SIMULATION,time_parameters)
%Esta funcion inicia y propaga el integrador de las ecuaciones del
%movimiento y de los diferentes bias, con actualizaciones de medidas GNSS.
%INPUT
%Xini: Condiciones iniciales
%AM: Medidas de las aceleraciones
%XSAT: Trayectorias de los satelites (Efemerides disponibles o calculadas )
%dSAT: Distancia de los satelites al vehiculo
%dt: Paso de tiempo
%dtGNSS: Paso de tiempo de las actualizaciones GNSS
%SAT: Objeto satelite, usamos unos pocos parametros para las matrices Kalman
%bias_param: Objeto bias, usamos unos pocos parametros para las matrices Kalman
%OUTPUT
%XK: Vector de estado en todos los instantes temporales
%residuos: Residuos calculados en la propagacion
%P: Matriz P de Kalman
Xini=WORLDMODEL.CI;
dt=time_parameters.dt;
%IMU
AM=IMU_SIMULATION.measures;
bias_param=IMU_SIMULATION.params;

sigma_GNSS=SAT.sigma_GNSS;
sigma_doppler=SAT.sigma_doppler; %Error del satélite

sigma_bias=bias_param.sigma_bias_a;
sigma_biasw=bias_param.sigma_bias_w;
n_satelites=size(SAT.position,2)/2;

%% ---------  MATRICES KALMAN  ---------------
R1=eye(n_satelites)*sigma_GNSS^2;
R2=eye(n_satelites)*sigma_doppler^2;
R=[R1 zeros(n_satelites);zeros(n_satelites) R2];

P(1,1)=10;P(2,2)=10;
P(5,5)=3*sigma_bias^2;P(6,6)=3*sigma_bias^2;
P(10,10)=3*sigma_biasw^2;P(11,11)=0;

sigma_v=sigma_bias*dt;
sigma_pos=sigma_bias*dt^2/2+sigma_v*dt;
sigma_ang=sigma_biasw*dt;

Q=zeros(11);
Q(1,1)=0.003;
Q(2,2)=0.003;
Q(3,3)=0.003;
Q(4,4)=0.003;
Q(5,5)=sigma_bias^2;
Q(6,6)=sigma_bias^2;
Q(9,9)=sigma_ang^2/100;
Q(10,10)=10*sigma_biasw^2;


Kalman_Matrix.P=P;
Kalman_Matrix.R=R;
Kalman_Matrix.Q=Q;

%% --------- INICIO DEL INTEGRADOR  --------------
%Se procede a propagar las ecuaciones diferenciales del movimiento con
%diferentes esquemas temporales
integrador='Euler';%Se puede elegir Kutta o Euler
kn=0;j=0;s=dt;T=length(AM);residuos(T)=0;P_buff(11,11,T)=0;
K_buff(11,n_satelites*2,T)=0;
Xini=[Xini(1) Xini(2) Xini(3) Xini(4) 0 0 0 0 Xini(7) 0 0]';
U=zeros(11,T);U(:,1)=Xini;m=0;imu_indicator(T)=0;
U2=zeros(11,T);U2(:,1)=Xini;dU=zeros(11,T);
U3=zeros(11,T);U3(:,1)=Xini;
fprintf('********************************************************\n');
if SAT.spoof_flag==0
    fprintf('Authentic Propagation Progress: \n');
elseif SAT.spoof_flag==1
    fprintf('Spoofed Propagation Progress: \n');
end
[threshold] = probability_interval(n_satelites,0.1);
%[threshold] = probability_interval(n_satelites,0.01);
IMU_SIMULATION.threshold=threshold;
if strcmp(integrador,'Euler')
    for k=2:T
        %Euler
        [K,Pkal,Kkal,j,residuos(k),IMU_SIMULATION] = mechanization(k-1,U(:,k-1),j,Kalman_Matrix,IMU_SIMULATION,SAT,time_parameters,WORLDMODEL);
        U(:,k)=U(:,k-1)+s*K;
        
        dU(:,k)=K;
        P_buff(:,:,k)=Pkal;
        if size(Kkal,2) == n_satelites*2
            K_buff(:,:,k)=Kkal;
        end
        Kalman_Matrix.P=Pkal;
        
        %Segunda propagacion con otra funcion
        [K] = mechanization_IMUonly(k-1,U2(:,k-1),IMU_SIMULATION);
        U2(:,k)=U2(:,k-1)+s*K;
        
        m=m+1;
        if m > 80
            m = 0;
            imu_indicator(k)=sqrt((U(1,k)-U2(1,k))^2+(U(2,k)-U2(2,k))^2);
            U2(:,k)=U(:,k);
        end
        
        %Tercera propagacion con otra funcion
        [K] = mechanization_GNSSonly(k-1,U(:,k-1),j,Kalman_Matrix,IMU_SIMULATION,SAT,time_parameters,WORLDMODEL);
        U3(:,k)=U3(:,k-1)+s*K;
        
        kn=kn+1;
        if kn>T/10
            fprintf('%.0f-',floor(100*k/T));
            kn=0;
        end
    end
    
elseif strcmp(integrador,'Kutta')
    for k=2:T
        %Runge Kutta 4 hasta el final
        [K1] = mechanization(k-1,U(:,k-1),j,Kalman_Matrix,IMU_SIMULATION,SAT,time_parameters,WORLDMODEL);
        [K2] = mechanization(k-1,U(:,k-1)+s*0.5*K1,j,Kalman_Matrix,IMU_SIMULATION,SAT,time_parameters,WORLDMODEL);
        [K3] = mechanization(k-1,U(:,k-1)+s*0.5*K2,j,Kalman_Matrix,IMU_SIMULATION,SAT,time_parameters,WORLDMODEL);
        [K4,Pkal,Kkal,j,residuos(k),IMU_SIMULATION] = mechanization(k-1,U(:,k-1)+s*K3,j,Kalman_Matrix,IMU_SIMULATION,SAT,time_parameters,WORLDMODEL);
        
        U(:,k)=U(:,k-1)+s/6*(K1+2*K2+2*K3+K4);
        
        dU(:,k)=K4;
        P_buff(:,:,k)=Pkal;
        P_buff(:,:,k)=Pkal;
        if size(Kkal,2) == n_satelites*2
            K_buff(:,:,k)=Kkal;
        end
        Kalman_Matrix.P=Pkal;
        
        %Segunda propagacion con otra funcion
        [K1] = mechanization_IMUonly(k-1,U2(:,k-1),IMU_SIMULATION);
        [K2] = mechanization_IMUonly(k-1,U2(:,k-1)+s*0.5*K1,IMU_SIMULATION);
        [K3] = mechanization_IMUonly(k-1,U2(:,k-1)+s*0.5*K2,IMU_SIMULATION);
        [K4] = mechanization_IMUonly(k-1,U2(:,k-1)+s*K3,IMU_SIMULATION);
        
        U2(:,k)=U2(:,k-1)+s/6*(K1+2*K2+2*K3+K4);
        
        m=m+1;
        if m > 80
            m = 0;
            imu_indicator(k)=sqrt((U(1,k)-U2(1,k))^2+(U(2,k)-U2(2,k))^2);
            U2(:,k)=U(:,k);
        end
        
        %Tercera propagacion con otra funcion
        [K1] = mechanization_GNSSonly(k-1,U3(:,k-1),j,Kalman_Matrix,IMU_SIMULATION,SAT,time_parameters,WORLDMODEL);
        [K2] = mechanization_GNSSonly(k-1,U3(:,k-1)+s*0.5*K1,j,Kalman_Matrix,IMU_SIMULATION,SAT,time_parameters,WORLDMODEL);
        [K3] = mechanization_GNSSonly(k-1,U3(:,k-1)+s*0.5*K2,j,Kalman_Matrix,IMU_SIMULATION,SAT,time_parameters,WORLDMODEL);
        [K4] = mechanization_GNSSonly(k-1,U3(:,k-1)+s*K3,j,Kalman_Matrix,IMU_SIMULATION,SAT,time_parameters,WORLDMODEL);
        
        U3(:,k)=U3(:,k-1)+s/6*(K1+2*K2+2*K3+K4);
        
        kn=kn+1;
        if kn>T/10
            fprintf('%.0f-',floor(100*k/T));
            kn=0;
        end
    end
end

XK=U;dXK=dU;
XGNSS=U3;
fprintf(' Done \n');

SIMULATION.state_vector=XK;
SIMULATION.gnss_state_vector=XGNSS;
SIMULATION.residuos=residuos;
SIMULATION.PMatrix=P_buff;
SIMULATION.KMatrix=K_buff;
SIMULATION.imu_only=imu_indicator;
SIMULATION.state_vector_derivate=dXK;
SIMULATION.IMU=IMU_SIMULATION;
SIMULATION.SAT=SAT;
end


function [threshold]=probability_interval(n_satelites,detection_probability)
%
N_div=1000;
x=linspace(0,15,N_div);
[y]=chifun(x,n_satelites);
for i=1:N_div
    if y(N_div-i)>detection_probability
        break
    end
end
threshold=x(N_div-i);
end

function [y]=chifun(x,nu)
%
y=zeros(length(x),1);
for i=1:length(x)
    y(i) = x(i)^((nu-2)/2)*exp(-x(i)/2)/gamma(nu/2)/2^(nu/2);
end
end