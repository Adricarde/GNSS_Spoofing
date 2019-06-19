function [WORLDMODEL] = tray_gen(trajectory_params,time_parameters)
%Esta funcion se encarga de general la trayectoria real del vehiculo, y
%calcular sus velocidades y aceleraciones en ejes inerciales. Tambien
%calcula las aceleraciones en ejes cuerpo.
%INPUTS
%trajectory_params: parametros de la trayectoria, como velocidad o numero
%de waypoints.
%dt: Paso de tiempo
%case_sel: String que contiene el nombre del tipo de trayectoria. Hasta
%ahora se dispone de circulo, worstcase y random.
%OUPUTS:
%XRI: Contiene la trayectoria, velocidad, aceleracion, orientacion y
%velocidad angular real del vehiculo.
%AB: Contiene las aceleraciones en ejes cuerpo Xb e Yb y velocidad angular en Zb
dt=time_parameters.dt;
T_alineamiento=time_parameters.alineamiento;
trajectory_params.radio_min = 50;%Radio Minimo
trajectory_params.height = 10000;%Altitud del vehículo
vehicle_speed=trajectory_params.vehicle_speed;%Velocidad del vehiculo
WORLDMODEL.vehicle_height=trajectory_params.height;%Velocidad del vehiculo
case_sel=trajectory_params.case;
%Casos de una trayectoria con los peores casos recogidos hasta ahora para
%generarla y otro caso en el que se genera aleatoriamente la trayectoria
%con un numero determinado de puntos
magnitude_order=2500;
if strcmp(case_sel, 'worstcase')
    XRW=[0 0; 1 1;0 2;1 3;2 3;1 4;5 5;3 2;5 3.5;2 3;3 5;5 5;7 5;7 7;8 6;9 5.1]*magnitude_order;
    XRWS=[0 0; 1 1;0 2;1 3;2 3;1 4;5 5;3 2;5 3.5;2 3;3 5;5 7;7 9;8 9;8 10;9 15]*magnitude_order;
    n_waypoints=size(XRW,1);
elseif strcmp(case_sel, 'random')
    n_waypoints=trajectory_params.n_waypoints;%Numero de puntos para el generador de trayectorias
    XRW(n_waypoints,2)=0;
    alpha(n_waypoints)=0;
    alpha(1)=0.1*pi*(3+rand*7);
    for i=2:n_waypoints
        radio=2+rand*3;
        alpha(i)=0.1*pi*(3+rand*7);
        alpha_sign=sign(2*rand-1);
        XRW(i,1)=XRW(i-1,1)+radio*cos(-alpha(i-1)+alpha_sign*alpha(i));
        XRW(i,2)=XRW(i-1,2)+radio*sin(-alpha(i-1)+alpha_sign*alpha(i));
        alpha(i)=-alpha(i-1)+alpha_sign*alpha(i);
    end
    XRWS=XRW;
    XRWS(n_waypoints,1)=XRWS(n_waypoints,1)+3;
    XRWS(n_waypoints,2)=XRWS(n_waypoints,2)+10;
    XRW=XRW*magnitude_order;
    XRWS=XRWS*magnitude_order;
elseif strcmp(case_sel, 'case1')
    XRW=[-3 0.1; -2.5 -1;2.5 -1; 4 0.1;2.5 2; -2.5 2;-3 1]*magnitude_order;
    XRWS=[-3 0.1; -2.5 -1;2.5 -1; 4 0.1;5 5; 15 15;25 25]*magnitude_order;
    n_waypoints=size(XRW,1);
elseif strcmp(case_sel, 'case2')
    XRW=[-3 0.1; -2.5 -1;2.5 -1; 4 0.1;2.5 2; -2.5 2;-3 0.1; -2.5 -1;2.5 -1;]*magnitude_order;
    XRWS=[-3 0.1; -2.5 -1;2.5 -1; 4 0.1;5 5; 15 15;25 25;35 35;45 45]*magnitude_order;
    n_waypoints=size(XRW,1);
elseif strcmp(case_sel, 'case3')
    XRW=[1 1;2 3; 1 4; 0 2; 1 1]*magnitude_order*2;
    XRWS=[1 1;2 3; 5 5; 7 7; 10 10]*magnitude_order*2;
    n_waypoints=size(XRW,1);
elseif strcmp(case_sel, 'case4')
    %     XRW=[1 1;2 1; 3 2; 4 1; 5 0]*10000;
    %     XRWS=[1 1;2 1; 3 3; 4 4; 5 5]*10000;
    n_vueltas=3;
    XRW((n_vueltas-1)*6+1,2)=0;
    XRWS((n_vueltas-1)*6+1,2)=0;
    for i=1:n_vueltas
        XRW((i-1)*6+1,1:2)=[0 0];
        XRW((i-1)*6+2,1:2)=[1 1];
        XRW((i-1)*6+3,1:2)=[1 2];
        XRW((i-1)*6+4,1:2)=[0 3];
        XRW((i-1)*6+5,1:2)=[-1 2];
        XRW((i-1)*6+6,1:2)=[-1 1];
        
        XRWS((i-1)*6+1,1:2)=[(i-1)*6 (i-1)*6];
        XRWS((i-1)*6+2,1:2)=[(i-1)*6+1 (i-1)*6+1];
        XRWS((i-1)*6+3,1:2)=[(i-1)*6+2 (i-1)*6+2];
        XRWS((i-1)*6+4,1:2)=[(i-1)*6+3 (i-1)*6+3];
        XRWS((i-1)*6+5,1:2)=[(i-1)*6+4 (i-1)*6+4];
        XRWS((i-1)*6+6,1:2)=[(i-1)*6+5 (i-1)*6+5];
    end
    XRW=XRW*5000;
    XRWS=XRWS*5000;
    n_waypoints=size(XRW,1);
elseif strcmp(case_sel, 'case5')
    XRW=[1 1;2 1; 3 2; 4 2; 5 2]*10000;
    XRWS=[1 1;2 1; 3 3; 4 4; 5 5]*10000;
    n_waypoints=size(XRW,1);
    
end
%Comienzo del generador de trayectoria que une tres puntos mediante
%un arco de circuferencia. Se realizan los calculos necesarios y se
%calcula la circunferencia tangente a las dos rectas que unen los
%tres puntos. Con la velocidad del vehiculo se discretiza la
%trayectoria y al ser rectas y circunferencias se sabe en todo
%momento las aceleraciones en ejes cuerpo.
vehicle_accel(n_waypoints-1,2)=0;
vehicle_accel(n_waypoints-1,1)=0.;
vehicle_accel(n_waypoints-1,2)=0;
% vehicle_accel(1,1)=0.50;
% vehicle_accel(2,1)=-0.10;
% vehicle_accel(3,1)=+0.50;
fprintf('********************************************************\n');
fprintf('Trajectory Generation...\n');
[XRI,AB] = state_generation(XRW,dt,T_alineamiento,vehicle_speed,vehicle_accel,trajectory_params.radio_min);
fprintf('Spoofed Trajectory Generation...\n');
[XRIS,ABS] = state_generation(XRWS,dt,T_alineamiento,vehicle_speed,vehicle_accel,trajectory_params.radio_min);
fprintf('Done\n');

if size(XRI,1) > size(XRIS,1)
    XRI=XRI(1:size(XRIS,1),:);
else
    XRIS=XRIS(1:size(XRI,1),:);
end
WORLDMODEL.true_tray=XRI;
WORLDMODEL.spoof_tray=XRIS;
WORLDMODEL.forces=AB;
WORLDMODEL.CI=XRI(1,:)';

if trajectory_params.plot
    figure(200)
    hold on
    plot(XRW(:,1)./1000,XRW(:,2)./1000);
    scatter(XRW(:,1)./1000,XRW(:,2)./1000);
    plot(XRI(:,1)./1000,XRI(:,2)./1000), grid on;
    %     plot(XRIS(:,1)./1000,XRIS(:,2)./1000), axis tight, grid on;
    xlabel('X (km)')
    ylabel('Y (km)')
    title('Trajectory')
    set(gca, 'FontName', 'Cambria','Fontsize',20)
    
    %axis tight
    axis equal
    pause(0.5)
end
end



function [XRI,AB] = state_generation(XRW,dt,T_align,vehicle_speed,vehicle_accel,Rmin)
XRI(1e8,8)=0;
AB(1e8,3)=0;
n_waypoints=size(XRW,1);
factor_distancia=0.75;k=1;

for i=1:n_waypoints-2
    d1=sqrt((XRW(i+1,1)-XRW(i,1))^2+(XRW(i+1,2)-XRW(i,2))^2);
    d2=sqrt((XRW(i+2,1)-XRW(i+1,1))^2+(XRW(i+2,2)-XRW(i+1,2))^2);
    theta1=atan2(XRW(i+1,2)-XRW(i,2),XRW(i+1,1)-XRW(i,1));
    theta2=atan2(XRW(i+2,2)-XRW(i+1,2),XRW(i+2,1)-XRW(i+1,1));
    
    xa=(XRW(i+1,1)-XRW(i,1))*factor_distancia+XRW(i,1);
    ya=(XRW(i+1,2)-XRW(i,2))*factor_distancia+XRW(i,2);
    xb=XRW(i+1,1)+d1*(1-factor_distancia)*cos(theta2);
    yb=XRW(i+1,2)+d1*(1-factor_distancia)*sin(theta2);
    
    a1=-(XRW(i+1,1)-XRW(i,1))/(XRW(i+1,2)-XRW(i,2));
    a2=-(XRW(i+2,1)-XRW(i+1,1))/(XRW(i+2,2)-XRW(i+1,2));
    b1=ya-a1*xa;
    b2=yb-a2*xb;
    
    if 1/a1 == 0
        xR=xa;
        yR=a2*xR+b2;
    elseif 1/a2 == 0
        xR=xb;
        yR=a1*xR+b1;
    else
        xR=(b1-b2)/(a2-a1);
        yR=a1*xR+b1;
    end
    R=sqrt((xR-xa)^2+(yR-ya)^2);
    if R < Rmin
        fprintf('***--------------ALARM---------------***\n');
        fprintf('       Trajectoy radius too small       \n');
        fprintf('***----------------------------------***\n');
    end
end

for k=1:T_align
    theta1=atan2(XRW(2,2)-XRW(1,2),XRW(2,1)-XRW(1,1));
    XRI(k,1)=XRW(1,1);
    XRI(k,2)=XRW(1,2);
    XRI(k,3)=0;
    XRI(k,4)=0;
    XRI(k,5)=0;XRI(k,6)=0;
    XRI(k,7)=theta1;XRI(k,8)=0;
    AB(k,1)=0;AB(k,2)=0;AB(k,3)=0;
end
accel_ini=2.1;
t_accel_ini=vehicle_speed/accel_ini;
N_i=floor(t_accel_ini/dt);
for j=1:N_i
    theta1=atan2(XRW(2,2)-XRW(1,2),XRW(2,1)-XRW(1,1));
    k=k+1;
    XRI(k,1)=XRW(1,1)+cos(theta1)*0.5*accel_ini*(j*dt)^2;
    XRI(k,2)=XRW(1,2)+sin(theta1)*0.5*accel_ini*(j*dt)^2;
    XRI(k,3)=cos(theta1)*accel_ini*(j*dt);
    XRI(k,4)=sin(theta1)*accel_ini*(j*dt);
    XRI(k,5)=cos(theta1)*accel_ini;XRI(k,6)=sin(theta1)*accel_ini;
    XRI(k,7)=theta1;XRI(k,8)=0;
    AB(k,1)=accel_ini;AB(k,2)=0;AB(k,3)=0;
end
T_init=k;
d_alin=sqrt((XRI(k,1)-XRI(1,1))^2+(XRI(k,2)-XRI(1,2))^2);
fprintf('If minus, Init Accel too small:  %.2f   \n',d1-d_alin);

k=k+1;
for i=1:n_waypoints-2
    d1=sqrt((XRW(i+1,1)-XRW(i,1))^2+(XRW(i+1,2)-XRW(i,2))^2);
    d2=sqrt((XRW(i+2,1)-XRW(i+1,1))^2+(XRW(i+2,2)-XRW(i+1,2))^2);
    theta1=atan2(XRW(i+1,2)-XRW(i,2),XRW(i+1,1)-XRW(i,1));
    theta2=atan2(XRW(i+2,2)-XRW(i+1,2),XRW(i+2,1)-XRW(i+1,1));
    
    xa=(XRW(i+1,1)-XRW(i,1))*factor_distancia+XRW(i,1);
    ya=(XRW(i+1,2)-XRW(i,2))*factor_distancia+XRW(i,2);
    xb=XRW(i+1,1)+d1*(1-factor_distancia)*cos(theta2);
    yb=XRW(i+1,2)+d1*(1-factor_distancia)*sin(theta2);
    
    a1=-(XRW(i+1,1)-XRW(i,1))/(XRW(i+1,2)-XRW(i,2));
    a2=-(XRW(i+2,1)-XRW(i+1,1))/(XRW(i+2,2)-XRW(i+1,2));
    b1=ya-a1*xa;
    b2=yb-a2*xb;
    
    if 1/a1 == 0
        xR=xa;
        yR=a2*xR+b2;
    elseif 1/a2 == 0
        xR=xb;
        yR=a1*xR+b1;
    else
        xR=(b1-b2)/(a2-a1);
        yR=a1*xR+b1;
    end
    if theta1 < 0
        theta1_aux=theta1+2*pi;
    else
        theta1_aux=theta1;
    end
    if theta2 < 0
        theta2_aux=theta2+2*pi;
    else
        theta2_aux=theta2;
    end
    R=sqrt((xR-xa)^2+(yR-ya)^2);
    
    thetaR=(theta2-theta1);
    thetaR_aux=(theta2_aux-theta1_aux);
    
    if thetaR < -pi
        thetaR = 2*pi+thetaR;
    end
    if thetaR > pi
        thetaR = -2*pi+thetaR;
    end
    vehicle_speed=sqrt(XRI(k-1,3)^2+XRI(k-1,4)^2);
    if thetaR ~= 0
        if vehicle_accel(i,2) == 0
            t_total=(R*abs(thetaR))/vehicle_speed;
        else
            t_total=(sqrt(vehicle_speed^2+2*vehicle_accel(i,2)*(R*abs(thetaR)))-vehicle_speed)/(vehicle_accel(i,2));
        end
        NR=floor(t_total/dt);
        
        %NR=floor(R*abs(thetaR)/(ship_speed*dt+0.5*ship_accel(i,2)*dt^2));
        delta_thetaR=0;delta_thetaR(NR)=0;
        for j=1:NR
            delta_thetaR(j)=abs(vehicle_speed*(j*dt)+0.5*vehicle_accel(i,2)*(dt*j)^2)/R;
        end
        omega=0;omega(NR)=0;
        if thetaR_aux  > pi || (thetaR_aux  < 0 && thetaR_aux  > -pi)
            for j = 1:NR
                omega(j)=-abs((vehicle_speed+(j-1)*vehicle_accel(i,2)*dt)/R);
            end
        else
            for j = 1:NR
                omega(j)=abs((vehicle_speed+(j-1)*vehicle_accel(i,2)*dt)/R);
            end
        end
    end
    
    if i==1
        if thetaR ==0
            if vehicle_accel(i,1) == 0
                t_total=(d1*(2-factor_distancia)-d_alin)/vehicle_speed;
            else
                t_total=(sqrt(vehicle_speed^2+2*vehicle_accel(i,1)*(d1*(2-factor_distancia)-d_alin))-vehicle_speed)/(vehicle_accel(i,1));
            end
            N_i=floor(t_total/dt);
            %N_i=floor((d1*(2-factor_distancia)-d_alin)/(ship_speed*dt+0.5*ship_accel(i,1)*dt^2));
            for j=1:N_i
                XRI(k,1)=XRI(T_init,1)+j*vehicle_speed*dt*cos(theta1)+0.5*vehicle_accel(i,1)*(j*dt)^2*cos(theta1);
                XRI(k,2)=XRI(T_init,2)+j*vehicle_speed*dt*sin(theta1)+0.5*vehicle_accel(i,1)*(j*dt)^2*sin(theta1);
                XRI(k,3)=vehicle_speed*cos(theta1)+j*vehicle_accel(i,1)*dt*cos(theta1);
                XRI(k,4)=vehicle_speed*sin(theta1)+j*vehicle_accel(i,1)*dt*sin(theta1);
                XRI(k,5)=vehicle_accel(i,1)*cos(theta1);XRI(k,6)=vehicle_accel(i,1)*sin(theta1);
                XRI(k,7)=theta1;XRI(k,8)=0;
                AB(k,1)=vehicle_accel(i,1);
                AB(k,2)=0;
                AB(k,3)=0;
                k=k+1;
            end
        else
            if vehicle_accel(i,1) == 0
                t_total=(d1*factor_distancia-d_alin)/vehicle_speed;
            else
                t_total=(sqrt(vehicle_speed^2+2*vehicle_accel(i,1)*(d1*factor_distancia-d_alin))-vehicle_speed)/(vehicle_accel(i,1));
            end
            N_i=floor(t_total/dt);
            %N_i=floor((d1*factor_distancia-d_alin)/(ship_speed*dt+0.5*ship_accel(i,1)*dt^2));
            for j=1:N_i
                XRI(k,1)=XRI(T_init,1)+j*vehicle_speed*dt*cos(theta1)+0.5*vehicle_accel(i,1)*(j*dt)^2*cos(theta1);
                XRI(k,2)=XRI(T_init,2)+j*vehicle_speed*dt*sin(theta1)+0.5*vehicle_accel(i,1)*(j*dt)^2*sin(theta1);
                XRI(k,3)=vehicle_speed*cos(theta1)+j*vehicle_accel(i,1)*dt*cos(theta1);
                XRI(k,4)=vehicle_speed*sin(theta1)+j*vehicle_accel(i,1)*dt*sin(theta1);
                XRI(k,5)=vehicle_accel(i,1)*cos(theta1);XRI(k,6)=vehicle_accel(i,1)*sin(theta1);
                XRI(k,7)=theta1;XRI(k,8)=0;
                AB(k,1)=vehicle_accel(i,1);
                AB(k,2)=0;
                AB(k,3)=0;
                k=k+1;
            end
            vehicle_speed=sqrt(XRI(k-1,3)^2+XRI(k-1,4)^2);
            if vehicle_accel(i,2) == 0
                t_total=(R*abs(thetaR))/vehicle_speed;
            else
                t_total=(sqrt(vehicle_speed^2+2*vehicle_accel(i,2)*(R*abs(thetaR)))-vehicle_speed)/(vehicle_accel(i,2));
            end
            NR=floor(t_total/dt);
            
            delta_thetaR=0;delta_thetaR(NR)=0;
            for j=1:NR
                delta_thetaR(j)=abs(vehicle_speed*(j*dt)+0.5*vehicle_accel(i,2)*(dt*j)^2)/R;
            end
            omega=0;omega(NR)=0;
            if thetaR_aux  > pi || (thetaR_aux  < 0 && thetaR_aux  > -pi)
                for j = 1:NR
                    omega(j)=-abs((vehicle_speed+(j-1)*vehicle_accel(i,2)*dt)/R);
                end
            else
                for j = 1:NR
                    omega(j)=abs((vehicle_speed+(j-1)*vehicle_accel(i,2)*dt)/R);
                end
            end
            for j=1:NR
                XRI(k,1)=xR+sign(omega(j))*R*sin(theta1+sign(omega(j))*delta_thetaR(j));
                XRI(k,2)=yR-sign(omega(j))*R*cos(theta1+sign(omega(j))*delta_thetaR(j));
                XRI(k,3)=sign(omega(j))*R*omega(j)*cos(theta1+sign(omega(j))*delta_thetaR(j));
                XRI(k,4)=sign(omega(j))*R*omega(j)*sin(theta1+sign(omega(j))*delta_thetaR(j));
                XRI(k,5)=-sign(omega(j))*R*omega(j)*omega(j)*sin(theta1+sign(omega(j))*delta_thetaR(j));
                XRI(k,6)=sign(omega(j))*R*omega(j)*omega(j)*cos(theta1+sign(omega(j))*delta_thetaR(j));
                XRI(k,7)=theta1+sign(omega(j))*delta_thetaR(j);
                XRI(k,8)=omega(j);
                AB(k,1)=vehicle_accel(i,2);
                AB(k,2)=sign(omega(j))*omega(j)^2*R;
                AB(k,3)=omega(j);
                k=k+1;
            end
        end
        vehicle_speed=sqrt(XRI(k-1,3)^2+XRI(k-1,4)^2);
        if i==n_waypoints-2
            if vehicle_accel(i+1,1) == 0
                t_total=(d2-d1*(1-factor_distancia))/vehicle_speed;
            else
                t_total=(sqrt(vehicle_speed^2+2*vehicle_accel(i+1,1)*(d2-d1*(1-factor_distancia)))-vehicle_speed)/(vehicle_accel(i+1,1));
            end
            N_i=floor(t_total/dt);
            %N_i=floor((d2-d1*(1-factor_distancia))/(ship_speed*dt+0.5*ship_accel(i+1,1)*dt^2));
            for j=1:N_i
                XRI(k,1)=xb+j*vehicle_speed*dt*cos(theta2)+0.5*vehicle_accel(i+1,1)*(j*dt)^2*cos(theta2);
                XRI(k,2)=yb+j*vehicle_speed*dt*sin(theta2)+0.5*vehicle_accel(i+1,1)*(j*dt)^2*sin(theta2);
                XRI(k,3)=vehicle_speed*cos(theta2)+j*vehicle_accel(i+1,1)*dt*cos(theta2);
                XRI(k,4)=vehicle_speed*sin(theta2)+j*vehicle_accel(i+1,1)*dt*sin(theta2);
                XRI(k,5)=vehicle_accel(i+1,1)*cos(theta2);XRI(k,6)=vehicle_accel(i+1,1)*sin(theta2);
                XRI(k,7)=theta2;XRI(k,8)=0;
                AB(k,1)=vehicle_accel(i+1,1);
                AB(k,2)=0;
                AB(k,3)=0;
                k=k+1;
            end
        else
            if vehicle_accel(i+1,1) == 0
                t_total=(d2-d1*(1-factor_distancia)-d2*(1-factor_distancia))/vehicle_speed;
            else
                t_total=(sqrt(vehicle_speed^2+2*vehicle_accel(i+1,1)*(d2-d1*(1-factor_distancia)-d2*(1-factor_distancia)))-vehicle_speed)/(vehicle_accel(i+1,1));
            end
            N_i=floor(t_total/dt);
            %N_i=floor((d2-d1*(1-factor_distancia)-d2*(1-factor_distancia))/(ship_speed*dt+0.5*ship_accel(i+1,1)*dt^2));
            for j=1:N_i
                XRI(k,1)=xb+j*vehicle_speed*dt*cos(theta2)+0.5*vehicle_accel(i+1,1)*(j*dt)^2*cos(theta2);
                XRI(k,2)=yb+j*vehicle_speed*dt*sin(theta2)+0.5*vehicle_accel(i+1,1)*(j*dt)^2*sin(theta2);
                XRI(k,3)=vehicle_speed*cos(theta2)+j*vehicle_accel(i+1,1)*dt*cos(theta2);
                XRI(k,4)=vehicle_speed*sin(theta2)+j*vehicle_accel(i+1,1)*dt*sin(theta2);
                XRI(k,5)=vehicle_accel(i+1,1)*cos(theta2);XRI(k,6)=vehicle_accel(i+1,1)*sin(theta2);
                XRI(k,7)=theta2;XRI(k,8)=0;
                AB(k,1)=vehicle_accel(i+1,1);
                AB(k,2)=0;
                AB(k,3)=0;
                k=k+1;
            end
        end
        
    else
        %i > 1
        if thetaR ==0
            if vehicle_accel(i+1,1) == 0
                t_total=(2*d1*(1-factor_distancia))/vehicle_speed;
            else
                t_total=(sqrt(vehicle_speed^2+2*vehicle_accel(i+1,1)*(2*d1*(1-factor_distancia)))-vehicle_speed)/(vehicle_accel(i+1,1));
            end
            N_i=floor(t_total/dt);
            %N_i=floor((2*d1*(1-factor_distancia))/(ship_speed*dt+0.5*ship_accel(i,1)*dt^2));
            for j=1:N_i
                XRI(k,1)=xa+j*vehicle_speed*dt*cos(theta1)+0.5*vehicle_accel(i,1)*(j*dt)^2*cos(theta1);
                XRI(k,2)=ya+j*vehicle_speed*dt*sin(theta1)+0.5*vehicle_accel(i,1)*(j*dt)^2*sin(theta1);
                XRI(k,3)=vehicle_speed*cos(theta1)+j*vehicle_accel(i,1)*dt*cos(theta1);
                XRI(k,4)=vehicle_speed*sin(theta1)+j*vehicle_accel(i,1)*dt*sin(theta1);
                XRI(k,5)=vehicle_accel(i,1)*cos(theta1);XRI(k,6)=vehicle_accel(i,1)*sin(theta1);
                XRI(k,7)=theta1;XRI(k,8)=0;
                AB(k,1)=vehicle_accel(i,1);
                AB(k,2)=0;
                AB(k,3)=0;
                k=k+1;
            end
        else
            %                                      figure(1)
            %                                      hold on
            %                 plot(XRI(1:k-1,1),XRI(1:k-1,2)), axis equal, grid on;
            
            for j=1:NR
                XRI(k,1)=xR+sign(omega(j))*R*sin(theta1+sign(omega(j))*delta_thetaR(j));
                XRI(k,2)=yR-sign(omega(j))*R*cos(theta1+sign(omega(j))*delta_thetaR(j));
                XRI(k,3)=sign(omega(j))*R*omega(j)*cos(theta1+sign(omega(j))*delta_thetaR(j));
                XRI(k,4)=sign(omega(j))*R*omega(j)*sin(theta1+sign(omega(j))*delta_thetaR(j));
                XRI(k,5)=-sign(omega(j))*R*omega(j)*omega(j)*sin(theta1+sign(omega(j))*delta_thetaR(j));
                XRI(k,6)=sign(omega(j))*R*omega(j)*omega(j)*cos(theta1+sign(omega(j))*delta_thetaR(j));
                XRI(k,7)=theta1+sign(omega(j))*delta_thetaR(j);
                XRI(k,8)=omega(j);
                AB(k,1)=vehicle_accel(i,2);
                AB(k,2)=sign(omega(j))*omega(j)^2*R;
                AB(k,3)=omega(j);
                %                         scatter(XRI(k,1),XRI(k,2)), axis tight, grid on;
                %                         pause(0.01)
                k=k+1;
            end
        end
        
        vehicle_speed=sqrt(XRI(k-1,3)^2+XRI(k-1,4)^2);
        if i==n_waypoints-2
            if vehicle_accel(i+1,1) == 0
                t_total=(d2-d1*(1-factor_distancia))/vehicle_speed;
            else
                t_total=(sqrt(vehicle_speed^2+2*vehicle_accel(i+1,1)*(d2-d1*(1-factor_distancia)))-vehicle_speed)/(vehicle_accel(i+1,1));
            end
            N_i=floor(t_total/dt);
            for j=1:N_i
                XRI(k,1)=xb+j*vehicle_speed*dt*cos(theta2)+0.5*vehicle_accel(i+1,1)*(j*dt)^2*cos(theta2);
                XRI(k,2)=yb+j*vehicle_speed*dt*sin(theta2)+0.5*vehicle_accel(i+1,1)*(j*dt)^2*sin(theta2);
                XRI(k,3)=vehicle_speed*cos(theta2)+j*vehicle_accel(i+1,1)*dt*cos(theta2);
                XRI(k,4)=vehicle_speed*sin(theta2)+j*vehicle_accel(i+1,1)*dt*sin(theta2);
                XRI(k,5)=vehicle_accel(i+1,1)*cos(theta2);XRI(k,6)=vehicle_accel(i+1,1)*sin(theta2);
                XRI(k,7)=theta2;XRI(k,8)=0;
                AB(k,1)=vehicle_accel(i+1,1);
                AB(k,2)=0;
                AB(k,3)=0;
                k=k+1;
            end
        else
            if vehicle_accel(i+1,1) == 0
                t_total=(d2-d1*(1-factor_distancia)-d2*(1-factor_distancia))/vehicle_speed;
            else
                t_total=(sqrt(vehicle_speed^2+2*vehicle_accel(i+1,1)*(d2-d1*(1-factor_distancia)-d2*(1-factor_distancia)))-vehicle_speed)/(vehicle_accel(i+1,1));
            end
            N_i=floor(t_total/dt);
            % N_i=floor((d2-d1*(1-factor_distancia)-d2*(1-factor_distancia))/(ship_speed*dt+0.5*ship_accel(i+1,1)*dt^2));
            alpha=0.5;
            beta=0.05;
            XP(N_i)=0;tP(N_i)=0;
            %             for j=1:N_i
            %                 XP(j)=200*exp(-(beta*(j-N_i*alpha)*dt)^2);
            %             end
            AP=diff(diff(XP))/dt^2;
            AP(N_i)=0;
            
            for j=1:N_i
                tP(j)=dt*j;
                
                XRI(k,1)=xb+j*vehicle_speed*dt*cos(theta2)+0.5*vehicle_accel(i+1,1)*(j*dt)^2*cos(theta2);
                XRI(k,2)=yb+j*vehicle_speed*dt*sin(theta2)+0.5*vehicle_accel(i+1,1)*(j*dt)^2*sin(theta2);
                XRI(k,3)=vehicle_speed*cos(theta2)+j*vehicle_accel(i+1,1)*dt*cos(theta2);
                XRI(k,4)=vehicle_speed*sin(theta2)+j*vehicle_accel(i+1,1)*dt*sin(theta2);
                XRI(k,5)=vehicle_accel(i+1,1)*cos(theta2);
                XRI(k,6)=vehicle_accel(i+1,1)*sin(theta2);
                %                 XRI(k,3)=(XRI(k,1)-XRI(k-1,1))/dt;
                %                 XRI(k,4)=(XRI(k,2)-XRI(k-1,2))/dt;
                %                 XRI(k,5)=(XRI(k,3)-XRI(k-1,3))/dt;
                %                 XRI(k,6)=(XRI(k,4)-XRI(k-1,4))/dt;
                
                %                 XRI(k,7)=atan2(XRI(k,4),XRI(k,3));
                XRI(k,7)=theta2;
                %                 XRI(k,8)=(XRI(k,7)-XRI(k-1,7))/dt;
                XRI(k,8)=0;
                AB(k,1)=vehicle_accel(i+1,1);
                AB(k,2)=0;
                AB(k,3)=XRI(k,8);
                k=k+1;
            end
            %plot(XRI(:,1),XRI(:,2)), grid on, axis tight;
            %             figure(10)
            %             hold on
            %             plot(tP,XP), grid on, axis tight;
            %             plot(tP,AP), grid on, axis tight;
            %             pause(0.01)
        end
    end
end
XRI=XRI(1:k-1,:);AB(k-1,3)=0;AB=AB(1:k-1,:);
% for i=2:k-1
% XRI(i,3)=(XRI(i,1)-XRI(i-1,1))/dt;
% XRI(i,4)=(XRI(i,2)-XRI(i-1,2))/dt;
% XRI(i,5)=(XRI(i,3)-XRI(i-1,3))/dt;
% XRI(i,6)=(XRI(i,4)-XRI(i-1,4))/dt;
% XRI(i,7)=atan2(XRI(i,4),XRI(i,3));
% XRI(i,8)=(XRI(i,7)-XRI(i-1,7))/dt;
% AB(i,1)=XRI(i,5)*cos(XRI(i,7))+XRI(i,6)*sin(XRI(i,7));
% AB(i,1)=-XRI(i,5)*sin(XRI(i,7))+XRI(i,6)*cos(XRI(i,7));
% AB(i,3)=XRI(i,8);
% end
end
