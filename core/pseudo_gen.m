function [SATPARAMS] = pseudo_gen(WORLDMODEL,SATPARAMS,tiempo)
%Esta funcion genera trayectorias rectilineas de un numero determinado
%de satelites uniendo dos puntos aleatorios generados entre las dimensiones
%de la trayectoria del vehiculo
%Tambien calcula la distancia de los satelites al vehiculo, e introduce
%errores en tiempos determinados para simular spoofing
%INPUTS
%XReal: Trayectoria, velocidad, aceleracion, orientacion y velocidad
%angular verdaderas a simular
%SATPARAMS: Objeto parametros del satelite
%OUTPUTS
%XSAT: Trayectoria de los satelites
%dSAT: Distancia de los satelites al vehiculo
XReal=WORLDMODEL.true_tray;
vehicle_height=WORLDMODEL.vehicle_height;
dt=tiempo.dt;dtGNSS=tiempo.dtGNSS;T_alineamiento=tiempo.alineamiento;T=size(XReal,1);
switch SATPARAMS.spoof_technique
    case 'none'
        SATPARAMS.coeff_spoof_0=0;%Offset del Spoofing
        SATPARAMS.coeff_spoof_1=0;%Termino lineal
    case 'delta'
        SATPARAMS.coeff_spoof_0=1e-4;%Offset del Spoofing
        SATPARAMS.coeff_spoof_1=0;%Termino lineal
    case 'progressive'
        SATPARAMS.coeff_spoof_0=0;%Offset del Spoofing
        SATPARAMS.coeff_spoof_1=1e-9;%Termino lineal
    case 'stealth_progressive'
        SATPARAMS.coeff_spoof_0=0;%Offset del Spoofing
        SATPARAMS.coeff_spoof_1=1e-11;%Termino lineal
    case 'hard_progressive'
        SATPARAMS.coeff_spoof_0=0;%Offset del Spoofing
        SATPARAMS.coeff_spoof_1=1e-3;%Termino lineal
    case 'complete'
        SATPARAMS.coeff_spoof_0=1e-3;%Offset del Spoofing
        SATPARAMS.coeff_spoof_1=1e-9;%Termino lineal
end
SATPARAMS.doppler_spoof=0;%Spoofing en Doppler
SATPARAMS.sigma_doppler=0.05; %Error del efecto Doppler
SATPARAMS.sigma_GNSS=4;%Error tipico en las actualizaciones GNSS
SATPARAMS.altitude=20e6;%
XSAT=SATPARAMS.position;
VSAT=SATPARAMS.velocity;

false_alarm_intensity=1.0;

n_satelites=SATPARAMS.n_satelites; % Numero de satelites
n_satelites_spoof=SATPARAMS.n_satelites_spoof; % Numero de satelites falsos
coeff_spoof0=SATPARAMS.coeff_spoof_0; % Offset del bias
coeff_spoof1=SATPARAMS.coeff_spoof_1; % Termino lineal del bias
sigma_GNSS=SATPARAMS.sigma_GNSS; % Error del satélite
sigma_doppler=SATPARAMS.sigma_doppler; % Error del satélite
SATPARAMS.lever_arm = 10;

fprintf('********************************************************\n');
fprintf('Satelites Pseudorange Generation...\n');

%Calculo de las distancias de los satelites al vehiculo
spoofing_type=SATPARAMS.spoof_type;
dSAT(T,n_satelites)=0;dSATs(T,n_satelites)=0;dopplerSAT(T,n_satelites)=0;dopplerSATs(T,n_satelites)=0;
switch spoofing_type
    case 'decalibrate'
        SATPARAMS.t_begin_spoof=floor(T*(0.25+rand*0.4));SATPARAMS.t_end_spoof=floor(3*T/3);
        t_begin_spoof=SATPARAMS.t_begin_spoof; %Inicio del Spoofing
        t_end_spoof=SATPARAMS.t_end_spoof; %Fin del Spoofing
        SATPARAMS.false_alarm_start = floor(SATPARAMS.t_begin_spoof*(0.4+rand*0.3));
        SATPARAMS.false_alarm_end = SATPARAMS.t_begin_spoof*(1+rand*0.3);
        for i=1:T
            if i > SATPARAMS.false_alarm_start && i < SATPARAMS.false_alarm_end
                sigma_GNSS=SATPARAMS.sigma_GNSS*false_alarm_intensity; % Error del satélite
            else
                sigma_GNSS=SATPARAMS.sigma_GNSS; % Error del satélite
            end
            for j=1:n_satelites
                %                 dSAT(i,j)=sqrt((XSAT(i,1+2*(j-1))-XReal(i,1))^2+(XSAT(i,2+2*(j-1))-XReal(i,2))^2 ...
                %                     +(SATPARAMS.altitude-vehicle_height)^2)+sigma_GNSS*(2*rand-1);
                dSAT(i,j)=sqrt((XSAT(i,1+2*(j-1))-(XReal(i,1)-SATPARAMS.lever_arm*cos(XReal(i,7))))^2 ...
                    +(XSAT(i,2+2*(j-1))-(XReal(i,2)-SATPARAMS.lever_arm*sin(XReal(i,7))))^2 ...
                    +(SATPARAMS.altitude-vehicle_height)^2)+sigma_GNSS*(2*rand-1);
                %n_vect=([XSAT(i,1+2*(j-1))-XReal(i,1) XSAT(i,2+2*(j-1))-XReal(i,2)]);
                n_vect=([XSAT(i,1+2*(j-1))-(XReal(i,1)-SATPARAMS.lever_arm*cos(XReal(i,7))) XSAT(i,2+2*(j-1))-(XReal(i,2)-SATPARAMS.lever_arm*sin(XReal(i,7)))]);
                n_vect=n_vect/norm(n_vect,2);
                %                 dopplerSAT(i,j)=dot_mod(XReal(i,3:4)-([VSAT(i,1+2*(j-1)) VSAT(i,2+2*(j-1))]),n_vect)+sigma_doppler*(2*rand-1);
                dopplerSAT(i,j)=dot_mod(XReal(i,3:4)-([VSAT(i,1+2*(j-1))-...
                    SATPARAMS.lever_arm*sin(XReal(i,7))*XReal(i,8) VSAT(i,2+2*(j-1))+...
                    SATPARAMS.lever_arm*cos(XReal(i,7))*XReal(i,8)]),n_vect)+sigma_doppler*(2*rand-1);
                
                if i>t_begin_spoof && i<t_end_spoof && j<=n_satelites_spoof
                    dSATs(i,j)=(1+coeff_spoof0)*dSAT(i,j);
                    if SATPARAMS.doppler_spoof == 1
                        dopplerSATs(i,j)=(1+coeff_spoof0)*dopplerSAT(i,j);
                    else
                        dopplerSATs(i,j)=dopplerSAT(i,j);
                    end
                    if coeff_spoof0 < 15
                        coeff_spoof0=coeff_spoof0+coeff_spoof1;
                    end
                else
                    dSATs(i,j)=dSAT(i,j);
                    dopplerSATs(i,j)=dopplerSAT(i,j);
                end
                
            end
        end
    case 'trajectory'
        XSpoof=WORLDMODEL.spoof_tray;
        for i=1:T
            if XReal(i,1) ~= XSpoof(i,1)
                SATPARAMS.t_begin_spoof=i;SATPARAMS.t_end_spoof=floor(3*T/3);
                break
            end
        end
        SATPARAMS.false_alarm_start = floor(SATPARAMS.t_begin_spoof*(0.5+rand*0.3));
        SATPARAMS.false_alarm_end = SATPARAMS.t_begin_spoof*(1+rand*0.3);
        for i=1:T
            if i > SATPARAMS.false_alarm_start && i < SATPARAMS.false_alarm_end
                sigma_GNSS=SATPARAMS.sigma_GNSS*false_alarm_intensity; % Error del satélite
            else
                sigma_GNSS=SATPARAMS.sigma_GNSS; % Error del satélite
            end
            for j=1:n_satelites
                %                 dSAT(i,j)=sqrt((XSAT(i,1+2*(j-1))-XReal(i,1))^2+(XSAT(i,2+2*(j-1))-XReal(i,2))^2 ...
                %                     +(SATPARAMS.altitude-vehicle_height)^2)+sigma_GNSS*(2*rand-1);
                dSAT(i,j)=sqrt((XSAT(i,1+2*(j-1))-(XReal(i,1)-SATPARAMS.lever_arm*cos(XReal(i,7))))^2 ...
                    +(XSAT(i,2+2*(j-1))-(XReal(i,2)-SATPARAMS.lever_arm*sin(XReal(i,7))))^2 ...
                    +(SATPARAMS.altitude-vehicle_height)^2)+sigma_GNSS*(2*rand-1);
                
                %                 n_vect=([XSAT(i,1+2*(j-1))-XReal(i,1) XSAT(i,2+2*(j-1))-XReal(i,2)]);
                %                 n_vect=n_vect/norm(n_vect,2);
                %                 %dopplerSAT(i,j)=dot_mod(XReal(i,3:4),n_vect)-dot_mod(([VSAT(i,1+2*(j-1)) VSAT(i,2+2*(j-1))]),n_vect)+sigma_doppler*(2*rand-1);
                %                 dopplerSAT(i,j)=dot_mod(XReal(i,3:4)-([VSAT(i,1+2*(j-1)) VSAT(i,2+2*(j-1))]),n_vect)+sigma_doppler*(2*rand-1);
                n_vect=([XSAT(i,1+2*(j-1))-(XReal(i,1)-SATPARAMS.lever_arm*cos(XReal(i,7))) XSAT(i,2+2*(j-1))-(XReal(i,2)-SATPARAMS.lever_arm*sin(XReal(i,7)))]);
                n_vect=n_vect/norm(n_vect,2);
                %                 dopplerSAT(i,j)=dot_mod(XReal(i,3:4)-([VSAT(i,1+2*(j-1)) VSAT(i,2+2*(j-1))]),n_vect)+sigma_doppler*(2*rand-1);
                dopplerSAT(i,j)=dot_mod(XReal(i,3:4)-([VSAT(i,1+2*(j-1))-...
                    SATPARAMS.lever_arm*sin(XReal(i,7))*XReal(i,8) VSAT(i,2+2*(j-1))+...
                    SATPARAMS.lever_arm*cos(XReal(i,7))*XReal(i,8)]),n_vect)+sigma_doppler*(2*rand-1);
                
                if  i>SATPARAMS.t_begin_spoof && i<SATPARAMS.t_end_spoof && j<=n_satelites_spoof
                    %                     dSATs(i,j)=sqrt((XSAT(i,1+2*(j-1))-XSpoof(i,1))^2+(XSAT(i,2+2*(j-1))-XSpoof(i,2))^2 ...
                    %                         +(SATPARAMS.altitude-vehicle_height)^2)+sigma_GNSS*(2*rand-1);
                    %                     n_vect=([XSAT(i,1+2*(j-1))-XSpoof(i,1) XSAT(i,2+2*(j-1))-XSpoof(i,2)]);
                    %                     n_vect=n_vect/norm(n_vect,2);
                    %                     dopplerSATs(i,j)=dot_mod(XSpoof(i,3:4),n_vect)-dot_mod(([VSAT(i,1+2*(j-1)) VSAT(i,2+2*(j-1))]),n_vect)+sigma_doppler*(2*rand-1);
                    %
                    dSATs(i,j)=sqrt((XSAT(i,1+2*(j-1))-(XSpoof(i,1)-SATPARAMS.lever_arm*cos(XSpoof(i,7))))^2 ...
                        +(XSAT(i,2+2*(j-1))-(XSpoof(i,2)-SATPARAMS.lever_arm*sin(XSpoof(i,7))))^2 ...
                        +(SATPARAMS.altitude-vehicle_height)^2)+sigma_GNSS*(2*rand-1);
                    n_vect=([XSAT(i,1+2*(j-1))-(XSpoof(i,1)-SATPARAMS.lever_arm*cos(XSpoof(i,7))) XSAT(i,2+2*(j-1))-(XSpoof(i,2)-SATPARAMS.lever_arm*sin(XReal(i,7)))]);
                    n_vect=n_vect/norm(n_vect,2);
                    %                 dopplerSAT(i,j)=dot_mod(XReal(i,3:4)-([VSAT(i,1+2*(j-1)) VSAT(i,2+2*(j-1))]),n_vect)+sigma_doppler*(2*rand-1);
                    dopplerSATs(i,j)=dot_mod(XSpoof(i,3:4)-([VSAT(i,1+2*(j-1))-...
                        SATPARAMS.lever_arm*sin(XSpoof(i,7))*XSpoof(i,8) VSAT(i,2+2*(j-1))+...
                        SATPARAMS.lever_arm*cos(XSpoof(i,7))*XSpoof(i,8)]),n_vect)+sigma_doppler*(2*rand-1);
                    
                    
                    %dopplerSATs(i,j)=dopplerSAT(i,j);
                    dSATs(i,j)=(1+coeff_spoof0)*dSAT(i,j);
                    %                     if SATPARAMS.doppler_spoof == 1
                    %                         dopplerSATs(i,j)=(1+coeff_spoof0)*dopplerSAT(i,j);
                    %                     else
                    %                         dopplerSATs(i,j)=dopplerSAT(i,j);
                    %                     end
                    if coeff_spoof0 < 15
                        coeff_spoof0=coeff_spoof0+coeff_spoof1;
                    end
                    
                else
                    dSATs(i,j)=dSAT(i,j);
                    dopplerSATs(i,j)=dopplerSAT(i,j);
                end
            end
        end
    case 'dangerzone'
        SATPARAMS.t_begin_spoof=floor(T*(0.25+rand*0.4));SATPARAMS.t_end_spoof=T;
        SATPARAMS.false_alarm_start = floor(SATPARAMS.t_begin_spoof*(0.4+rand*0.3));
        SATPARAMS.false_alarm_end = SATPARAMS.t_begin_spoof*(1+rand*0.3);
        t_begin_spoof=SATPARAMS.t_begin_spoof; %Inicio del Spoofing
        t_end_spoof=SATPARAMS.t_end_spoof; %Fin del Spoofing
        %SATPARAMS.pos_spoof=[min(XReal(:,1))-1e3 min(XReal(:,2))-1e3];
        SATPARAMS.pos_spoof=[max(XReal(:,1))+1e5 max(XReal(:,2))+1e5];
        for i=1:T
            if i > SATPARAMS.false_alarm_start && i < SATPARAMS.false_alarm_end
                sigma_GNSS=SATPARAMS.sigma_GNSS*false_alarm_intensity; % Error del satélite
            else
                sigma_GNSS=SATPARAMS.sigma_GNSS; % Error del satélite
            end
            for j=1:n_satelites
                %                 dSAT(i,j)=sqrt((XSAT(i,1+2*(j-1))-XReal(i,1))^2+(XSAT(i,2+2*(j-1))-XReal(i,2))^2 ...
                %                     +(SATPARAMS.altitude-vehicle_height)^2)+sigma_GNSS*(2*rand-1);
                %                 n_vect=([XSAT(i,1+2*(j-1))-XReal(i,1) XSAT(i,2+2*(j-1))-XReal(i,2)]);
                %                 n_vect=n_vect/norm(n_vect,2);
                %                 dopplerSAT(i,j)=dot_mod(XReal(i,3:4),n_vect)-dot_mod(([VSAT(i,1+2*(j-1)) VSAT(i,2+2*(j-1))]),n_vect)+sigma_doppler*(2*rand-1);
                %
                dSAT(i,j)=sqrt((XSAT(i,1+2*(j-1))-(XReal(i,1)-SATPARAMS.lever_arm*cos(XReal(i,7))))^2 ...
                    +(XSAT(i,2+2*(j-1))-(XReal(i,2)-SATPARAMS.lever_arm*sin(XReal(i,7))))^2 ...
                    +(SATPARAMS.altitude-vehicle_height)^2)+sigma_GNSS*(2*rand-1);
                n_vect=([XSAT(i,1+2*(j-1))-(XReal(i,1)-SATPARAMS.lever_arm*cos(XReal(i,7))) XSAT(i,2+2*(j-1))-(XReal(i,2)-SATPARAMS.lever_arm*sin(XReal(i,7)))]);
                n_vect=n_vect/norm(n_vect,2);
                dopplerSAT(i,j)=dot_mod(XReal(i,3:4)-([VSAT(i,1+2*(j-1))-...
                    SATPARAMS.lever_arm*sin(XReal(i,7))*XReal(i,8) VSAT(i,2+2*(j-1))+...
                    SATPARAMS.lever_arm*cos(XReal(i,7))*XReal(i,8)]),n_vect)+sigma_doppler*(2*rand-1);
                
                
                if  i>t_begin_spoof && i<t_end_spoof && j<=n_satelites_spoof
                    dSATs(i,j)=sqrt((XSAT(i,1+2*(j-1))-SATPARAMS.pos_spoof(1))^2+...
                        (XSAT(i,2+2*(j-1))-SATPARAMS.pos_spoof(2))^2+(SATPARAMS.altitude-vehicle_height)^2)+sigma_GNSS*(2*rand-1);
                    
                    dopplerSATs(i,j)=dopplerSAT(i,j);
                    
                    coeff_sign=0;
                    xspoof=SATPARAMS.pos_spoof(1);
                    xreal=XReal(i,1);
                    xsat=XSAT(i,1+2*(j-1));
                    yspoof=SATPARAMS.pos_spoof(2);
                    yreal=XReal(i,2);
                    ysat=XSAT(i,2+2*(j-1));
                    if xspoof > xreal
                        if yspoof > yreal
                            %ZONE 1, 2.1, 3
                            if ysat > yreal && xsat > xreal
                                coeff_sign=-1;
                            elseif ysat < yreal && xsat < xreal
                                coeff_sign=1;
                            end
                        else
                            if ysat < yreal && xsat > xreal
                                coeff_sign=-1;
                            elseif ysat > yreal && xsat < xreal
                                coeff_sign=1;
                            end
                        end
                    else
                        if yspoof > yreal
                            %ZONE 1, 2.1, 3
                            if ysat > yreal && xsat < xreal
                                coeff_sign=-1;
                            elseif ysat < yreal && xsat > xreal
                                coeff_sign=1;
                            end
                        else
                            if ysat < yreal && xsat < xreal
                                coeff_sign=-1;
                            elseif ysat > yreal && xsat > xreal
                                coeff_sign=1;
                            end
                        end
                    end
                    
                    dSATs(i,j)=(1+coeff_sign*coeff_spoof0)*dSATs(i,j);
                    if coeff_spoof0 < 15
                        coeff_spoof0=coeff_spoof0+coeff_spoof1;
                    end
                else
                    dSATs(i,j)=dSAT(i,j);
                    dopplerSATs(i,j)=dopplerSAT(i,j);
                end
            end
        end
    case 'sigmaincrement'
        SATPARAMS.t_begin_spoof=floor(T*(0.25+rand*0.4));SATPARAMS.t_end_spoof=floor(3*T/3);
        t_begin_spoof=SATPARAMS.t_begin_spoof; %Inicio del Spoofing
        t_end_spoof=SATPARAMS.t_end_spoof; %Fin del Spoofing
        SATPARAMS.false_alarm_start = floor(SATPARAMS.t_begin_spoof*(0.4+rand*0.3));
        SATPARAMS.false_alarm_end = SATPARAMS.t_begin_spoof*(1+rand*0.3);
        for i=1:T
            if i > SATPARAMS.false_alarm_start && i < SATPARAMS.false_alarm_end
                sigma_GNSS=SATPARAMS.sigma_GNSS*false_alarm_intensity; % Error del satélite
            else
                sigma_GNSS=SATPARAMS.sigma_GNSS; % Error del satélite
            end
            for j=1:n_satelites
                
                d=sqrt((XSAT(i,1+2*(j-1))-(XReal(i,1)-SATPARAMS.lever_arm*cos(XReal(i,7))))^2 ...
                    +(XSAT(i,2+2*(j-1))-(XReal(i,2)-SATPARAMS.lever_arm*sin(XReal(i,7))))^2 ...
                    +(SATPARAMS.altitude-vehicle_height)^2);
                dSAT(i,j)=d+sigma_GNSS*(2*rand-1);
                
                n_vect=([XSAT(i,1+2*(j-1))-(XReal(i,1)-SATPARAMS.lever_arm*cos(XReal(i,7))) XSAT(i,2+2*(j-1))-(XReal(i,2)-SATPARAMS.lever_arm*sin(XReal(i,7)))]);
                n_vect=n_vect/norm(n_vect,2);
                dopplerSAT(i,j)=dot_mod(XReal(i,3:4)-([VSAT(i,1+2*(j-1))-...
                    SATPARAMS.lever_arm*sin(XReal(i,7))*XReal(i,8) VSAT(i,2+2*(j-1))+...
                    SATPARAMS.lever_arm*cos(XReal(i,7))*XReal(i,8)]),n_vect)+sigma_doppler*(2*rand-1);
                
                if i>t_begin_spoof && i<t_end_spoof && j<=n_satelites_spoof
                    dSATs(i,j)=d+(1+coeff_spoof0)*sigma_GNSS*(2*rand-1);
                    if SATPARAMS.doppler_spoof == 1
                        dopplerSATs(i,j)=(1+coeff_spoof0)*dopplerSAT(i,j);
                    else
                        dopplerSATs(i,j)=dopplerSAT(i,j);
                    end
                    if coeff_spoof0 < 1000
                        coeff_spoof0=coeff_spoof0+coeff_spoof1;
                    end
                else
                    dSATs(i,j)=dSAT(i,j);
                    dopplerSATs(i,j)=dopplerSAT(i,j);
                end
                
            end
        end
end

fprintf('Final Spoofing Coefficient: %.4f  ...\n',coeff_spoof0);

fprintf('LSM GNSS Positions...\n');
warning('off','all');
[x,y,xs,ys,tGNSS] = position_determination_lsm(dSAT,dSATs,WORLDMODEL,SATPARAMS,tiempo,XSAT);
XGNSS=[x' y'];XGNSS_spoofed=[xs' ys'];
fprintf('Done\n');

SATPARAMS.plot=0;
if SATPARAMS.plot
    figure(12)
    close 12
    figure(12)
    hold on
    scatter(XGNSS(:,1),XGNSS(:,2)),grid on, axis tight;
    scatter(XGNSS_spoofed(:,1),XGNSS_spoofed(:,2)),grid on, axis tight;
    pause(0.001)
end

SATPARAMS.distance=dSAT;
SATPARAMS.distance_spoofed=dSATs;
SATPARAMS.doppler=dopplerSAT;
SATPARAMS.doppler_spoofed=dopplerSATs;
SATPARAMS.vehicle_position.x=XGNSS(:,1);
SATPARAMS.vehicle_position.y=XGNSS(:,2);
SATPARAMS.vehicle_position.x_spoofed=XGNSS_spoofed(:,1);
SATPARAMS.vehicle_position.y_spoofed=XGNSS_spoofed(:,2);
SATPARAMS.tGNSS=tGNSS;
end

function [result] = dot_mod(A,B)
result=A(1)*B(1)+A(2)*B(2);
end

function [x,y,xs,ys,t] = position_determination_lsm(dSAT,dSATs,WORLDMODEL,SATPARAMS,tiempo,XSAT)
XReal=WORLDMODEL.true_tray;
XSpoof=WORLDMODEL.spoof_tray;
vehicle_height=WORLDMODEL.vehicle_height;

dt=tiempo.dt;
dtGNSS=tiempo.dtGNSS;
T=size(XReal,1);
n_satelites=SATPARAMS.n_satelites; % Numero de satelites
sigma_GNSS=SATPARAMS.sigma_GNSS; % Error del satélite
b(n_satelites,1)=0;A(n_satelites,2)=0;delta(T,2)=0;x(T)=0;y(T)=0;
count=1;count2=1;t(T)=0;watchdog=0;
for k=1:T
    if count2 >= dtGNSS/dt
        count2=1;
        if count ==1
            x(count)=XReal(k,1);
            y(count)=XReal(k,2);
        else
            x(count)=x(count-1);
            y(count)=y(count-1);
        end
        delta(count,:)=[1000 1000]';
        while norm(delta(count,:),2) > 0.1
            for i=1:n_satelites
                rsat = norm([XSAT(k,1+2*(i-1))-x(count),XSAT(k,2+2*(i-1))-y(count),SATPARAMS.altitude-vehicle_height],2);
                b(i,1)= dSAT(k,i)-rsat;
                A(i,1) = (x(count) - XSAT(k,1+2*(i-1)))/rsat;
                A(i,2) = (y(count) - XSAT(k,2+2*(i-1)))/rsat;
            end
            aux=(inv(A'*A)*A'*b)';
            delta(count,:)=aux;
            x(count)= x(count) + delta(count,1);
            y(count)= y(count) + delta(count,2);
            watchdog=watchdog+1;
            if watchdog > 100
                break
            end
        end
        t(count)=k*dt;
        count=count+1;
        watchdog=0;
    end
    count2=count2+1;
end
x=x(1:count-1);y=y(1:count-1);t=t(1:count-1);
count=1;count2=1;

xs(T)=0;ys(T)=0;
for k=1:T
    if count2 >= dtGNSS/dt
        count2=1;
        if count ==1
            xs(count)=XReal(k,1);
            ys(count)=XReal(k,2);
        else
            xs(count)=xs(count-1);
            ys(count)=ys(count-1);
        end
        delta(count,:)=[1000 1000]';
        while norm(delta(count,:),2) > 0.1
            for i=1:n_satelites
                rsat = norm([XSAT(k,1+2*(i-1))-xs(count),XSAT(k,2+2*(i-1))-ys(count),SATPARAMS.altitude-vehicle_height],2);
                b(i,1)= dSATs(k,i)-rsat;
                A(i,1) = (xs(count) - XSAT(k,1+2*(i-1)))/rsat;
                A(i,2) = (ys(count) - XSAT(k,2+2*(i-1)))/rsat;
            end
            aux=(inv(A'*A)*A'*b)';
            delta(count,:)=aux;
            xs(count)= xs(count) + delta(count,1);
            ys(count)= ys(count) + delta(count,2);
            watchdog=watchdog+1;
            if watchdog > 100
                break
            end
        end
        count=count+1;
        watchdog=0;
        
    end
    count2=count2+1;
end
xs=xs(1:count-1);ys=ys(1:count-1);

end

