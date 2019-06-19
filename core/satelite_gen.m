function [SATPARAMS] = satelite_gen(WORLDMODEL,SATPARAMS,tiempo)
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
dt=tiempo.dt;dtGNSS=tiempo.dtGNSS;
T_alineamiento=tiempo.alineamiento;T=size(XReal,1);

try
    n_satelites_start=size(SATPARAMS.distance,2)+1;
    XSAT=SATPARAMS.position;
    VSAT=SATPARAMS.velocity;
catch
    n_satelites_start=1;
end
n_satelites=SATPARAMS.n_satelites; % Numero de satelites

fprintf('********************************************************\n');
fprintf('Satelites Positions Generation...\n');
sat_gen='onepoint';
V_order=3000;%m/s
if strcmp(sat_gen,'twopoint')
    XSAT(T,2*n_satelites)=0;sat_angle(n_satelites)=0;VSAT(T,2*n_satelites)=0;
    %xmagnitude_order=max([abs(max(XReal(:,1))) abs(min(XReal(:,1)))]);
    xmagnitude_order=V_order*dt*T/2;
    ymagnitude_order=max([abs(max(XReal(:,2))) abs(min(XReal(:,2))) xmagnitude_order]);
    %Generacion de trayectorias aleatorias de satelites
    for j=1:n_satelites
        XSAT(:,1+2*(j-1))=linspace(xmagnitude_order*(2*rand-1),xmagnitude_order*(2*rand-1),T);
        XSAT(:,2+2*(j-1))=linspace(ymagnitude_order*(2*rand-1),ymagnitude_order*(2*rand-1),T);
        sat_angle(j)=atan2(XSAT(end,2+2*(j-1))-XSAT(1,2+2*(j-1)),XSAT(end,1+2*(j-1))-XSAT(1,1+2*(j-1)));
        VSAT(:,1+2*(j-1))=sqrt((XSAT(end,2+2*(j-1))-XSAT(1,2+2*(j-1)))^2+...
            (XSAT(end,1+2*(j-1))-XSAT(1,1+2*(j-1)))^2)/(T*dt)*cos(sat_angle(j));
        VSAT(:,2+2*(j-1))=sqrt((XSAT(end,2+2*(j-1))-XSAT(1,2+2*(j-1)))^2+...
            (XSAT(end,1+2*(j-1))-XSAT(1,1+2*(j-1)))^2)/(T*dt)*sin(sat_angle(j));
    end
elseif strcmp(sat_gen,'onepoint')
    XSAT(T,2*n_satelites)=0;sat_angle(n_satelites)=0;VSAT(T,2*n_satelites)=0;
    %xmagnitude_order=max([abs(max(XReal(:,1))) abs(min(XReal(:,1)))]);
    xmagnitude_order=V_order*dt*T/2;
    ymagnitude_order=max([abs(max(XReal(:,2))) abs(min(XReal(:,2))) xmagnitude_order]);
    %Generacion de trayectorias aleatorias de satelites
    for j=n_satelites_start:n_satelites
        x1=sign(rand-0.5)*xmagnitude_order*(2*rand-1);
        y1=sign(rand-0.5)*ymagnitude_order*(2*rand-1);
        V_sat=V_order+(2*rand-1);theta1=1.2*rand;r=V_sat*dt*T;
        x2=x1+sign(rand-0.5)*r*cos(theta1);
        y2=y1+sign(rand-0.5)*r*sin(theta1);
        XSAT(:,1+2*(j-1))=linspace(x1,x2,T);
        XSAT(:,2+2*(j-1))=linspace(y1,y2,T);
        sat_angle(j)=atan2(XSAT(end,2+2*(j-1))-XSAT(1,2+2*(j-1)),XSAT(end,1+2*(j-1))-XSAT(1,1+2*(j-1)));
        VSAT(:,1+2*(j-1))=V_sat*cos(theta1);
        VSAT(:,2+2*(j-1))=V_sat*sin(theta1);
    end
elseif strcmp(sat_gen,'ellipses')
    V_order=3000;%m/s
    a=50*max(XReal(:,1));b=50*max(XReal(:,2));eps=sqrt(1-b*b/a/a);
    XSAT(T,2*n_satelites)=0;VSAT(T,2*n_satelites)=0;
    %Generacion de trayectorias aleatorias de satelites
    for j=1:n_satelites
        theta0=2*pi*rand;
        V_sat=V_order+(2*rand-1);
        for i=1:T
            r=b/sqrt(1-(eps*cos(theta0))^2);
            XSAT(i,1+2*(j-1))=a*cos(theta0);
            XSAT(i,2+2*(j-1))=b*sin(theta0);
            VSAT(i,1+2*(j-1))=-V_sat*a*sin(theta0)/r;
            VSAT(i,2+2*(j-1))=V_sat*b*cos(theta0)/r;
            theta0=theta0+V_sat/r*dt;
        end
    end
end

SATPARAMS.position=XSAT;
SATPARAMS.velocity=VSAT;

[SATPARAMS] = pseudo_gen(WORLDMODEL,SATPARAMS,tiempo);

end
