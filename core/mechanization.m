function [dy,P,K,j,residues,IMU_SIMULATION] = mechanization(k,y,j,Kalman_Matrix,IMU_SIMULATION,SAT,timeparameters,WORLDMODEL)
%Esta funcion genera la derivada temporal del vector de estado para
%inyectar en el esquema temporal
%INPUT
%k: Instante de tiempo (indice)
%y: Vector de estado
%j: Contador hasta la proxima medida GNSS
%P, Q, R: Matrices de Kalman
%AM: Medidas de los acelerometros en x e y, y velocidad angular z
%XSAT: Trayectorias de los satelites (Efemerides disponibles o calculadas )
%dSAT: Distancia de los satelites al vehiculo
%dt: Paso de tiempo
%dtGNSS: Paso de tiempo de las actualizaciones GNSS
%OUTPUT
%dy: Derivada temporal del vector de estado
%P: Matriz P de Kalman
%j: Contador hasta la proxima medida GNSS
%residuos: Residuos calculados en el instante k
dt=timeparameters.dt;
dtGNSS=timeparameters.dtGNSS;
T_align=timeparameters.alineamiento;
vehicle_height=WORLDMODEL.vehicle_height;
AM=IMU_SIMULATION.measures;

XSAT=SAT.position;
if SAT.spoof_flag
    dSAT=SAT.distance_spoofed;
    doppler=SAT.doppler_spoofed;
else
    dSAT=SAT.distance;
    doppler=SAT.doppler;
end
VSAT=SAT.velocity;

P=Kalman_Matrix.P;R=Kalman_Matrix.R;Q=Kalman_Matrix.Q;
n_satelites=size(XSAT,2)/2;
H(2*n_satelites,11)=0;K=0;
j=j+1;

%Transformacion de ejes cuerpo a ejes fijos
% AMX=(AM(k,1)-y(5)-y(7))*cos(y(9))-(AM(k,2)-y(6)-y(8))*sin(y(9));
% AMY=(AM(k,1)-y(5)-y(7))*sin(y(9))+(AM(k,2)-y(6)-y(8))*cos(y(9));
% AMW=AM(k,3)-y(10)-y(11);
AMX=(AM(k,1)-y(5))*cos(y(9))-(AM(k,2)-y(6))*sin(y(9));
AMY=(AM(k,1)-y(5))*sin(y(9))+(AM(k,2)-y(6))*cos(y(9));
AMW=AM(k,3)-y(10);

%Construcción de matriz A
% A=eye(11);
% A(1,3)=dt;A(1,5)=-dt^2/2*cos(y(9));A(1,6)=dt^2/2*sin(y(9));A(1,7)=-dt^2/2*cos(y(9));A(1,8)=dt^2/2*sin(y(9));
% A(1,9)=-dt^2/2*((AM(k,1)-y(5)-y(7))*sin(y(9))+(AM(k,2)-y(6)-y(8))*cos(y(9)));
% A(2,4)=dt;A(2,5)=-dt^2/2*sin(y(9));A(2,6)=-dt^2/2*cos(y(9));A(2,7)=-dt^2/2*sin(y(9));A(2,8)=-dt^2/2*cos(y(9));
% A(2,9)=dt^2/2*((AM(k,1)-y(5)-y(7))*cos(y(9))-(AM(k,2)-y(6)-y(8))*sin(y(9)));
% A(3,5)=-dt*cos(y(9));A(3,6)=dt*sin(y(9));A(3,7)=-dt*cos(y(9));A(3,8)=dt*sin(y(9));
% A(3,9)=-dt*((AM(k,1)-y(5)-y(7))*sin(y(9))+(AM(k,2)-y(6)-y(8))*cos(y(9)));
% A(4,5)=-dt*sin(y(9));A(4,6)=-dt*cos(y(9));A(4,7)=-dt*sin(y(9));A(4,8)=-dt*cos(y(9));
% A(4,9)=dt*((AM(k,1)-y(5)-y(7))*cos(y(9))-(AM(k,2)-y(6)-y(8))*sin(y(9)));
% A(9,10)=-dt;A(9,11)=-dt;
A=eye(11);
A(1,3)=dt;A(1,5)=-dt^2/2*cos(y(9));A(1,6)=dt^2/2*sin(y(9));
A(1,9)=-dt^2/2*((AM(k,1)-y(5))*sin(y(9))+(AM(k,2)-y(6))*cos(y(9)));
A(2,4)=dt;A(2,5)=-dt^2/2*sin(y(9));A(2,6)=-dt^2/2*cos(y(9));
A(2,9)=dt^2/2*((AM(k,1)-y(5))*cos(y(9))-(AM(k,2)-y(6))*sin(y(9)));
A(3,5)=-dt*cos(y(9));A(3,6)=dt*sin(y(9));
A(3,9)=-dt*((AM(k,1)-y(5))*sin(y(9))+(AM(k,2)-y(6))*cos(y(9)));
A(4,5)=-dt*sin(y(9));A(4,6)=-dt*cos(y(9));
A(4,9)=dt*((AM(k,1)-y(5))*cos(y(9))-(AM(k,2)-y(6))*sin(y(9)));
A(9,10)=-dt;



%Estimacion de la derivada
dy(11,1)=0;
% dy(1,1)=y(3);dy(2,1)=y(4);
% dy(3,1)=(AMX);dy(4,1)=(AMY);
% dy(9,1)=(AMW);
dy(1,1)=y(3);dy(2,1)=y(4);
dy(3,1)=(AMX);dy(4,1)=(AMY);
%dy(5,1)=y(7);
%dy(6,1)=y(8);
dy(9,1)=(AMW);
%dy(10,1)=y(11);


%Condicion para actualizacion GNSS
if k < T_align
%Calculo de Pmenos
P_minus=A*P*A'+Q;

%     Hstatic = [0 0 1 0 0 0 0 0 0 0 0;
%                0 0 0 1 0 0 0 0 0 0 0;
%                0 0 0 0 -1 0 -1 0 0 0 0;
%                0 0 0 0 0 -1 0 -1 0 0 0;
%                0 0 0 0 0 0 0 0 0 -1 -1];
    Hstatic = [0 0 1 0 0 0 0 0 0 0 0;
               0 0 0 1 0 0 0 0 0 0 0;
               0 0 0 0 -1 0 0 0 0 0 0;
               0 0 0 0 0 -1 0 0 0 0 0;
               0 0 0 0 0 0 0 0 0 -1 0];
    Rparado = [0 0 0 0 0;
               0 0 0 0 0;
               0 0 IMU_SIMULATION.params.sigma_bias_a 0 0;
               0 0 0 IMU_SIMULATION.params.sigma_bias_a 0;
               0 0 0 0 IMU_SIMULATION.params.sigma_bias_w]*0.1;
    S=Hstatic*P_minus*Hstatic'+Rparado;
    Sinv=inv(S);
    K=P_minus*Hstatic'*Sinv;
    P=(eye(length(dy))-K*Hstatic)*P_minus;
    %Y_measured = [0 0 AMX AMY AMW]';
    Y_measured = [0 0 0 0 0]';
    %Y_Xminus = [y(3) y(4) AM(k,1)-y(5)-y(7) AM(k,2)-y(6)-y(8) AM(k,3)-y(10)-y(11)]';
    %Y_Xminus = [y(3) y(4) AM(k,1)-y(5)-y(7) AM(k,2)-y(6)-y(8) AM(k,3)-y(10)-y(11)]';
    %Y_Xminus = [y(3) y(4) dy(3) dy(4) dy(9)]';
    Y_Xminus = [y(3) y(4) AM(k,1)-y(5) AM(k,2)-y(6) AM(k,3)-y(10)]';
    nu=Y_measured-Y_Xminus;
    %nu=Y_measured-Hstatic*y;
    inovation=K*(nu);
    residues=(nu)'*(Sinv)*(nu);
    dy=dy+inovation;
else
%Calculo de Pmenos
P_minus=A*P*A'+Q;
    if (j>=dtGNSS/dt)
        j=0;Y_measured(2*n_satelites,1)=0;Y_Xminus(2*n_satelites,1)=0;
        %Construccion de la matriz C y vectores de medidas y estimaciones
        for i=1:size(XSAT,2)
            if i <=size(XSAT,2)/2
%                 d1=sqrt((XSAT(k,1+2*(i-1))-y(1))^2+(XSAT(k,2+2*(i-1))-y(2))^2+(SAT.altitude-vehicle_height)^2);
%                 H(i,1)=-(XSAT(k,1+2*(i-1))-y(1))/d1;
%                 H(i,2)=-(XSAT(k,2+2*(i-1))-y(2))/d1;
                d1=sqrt((XSAT(k,1+2*(i-1))-(y(1)-SAT.lever_arm*cos(y(9))))^2+(XSAT(k,2+2*(i-1))-(y(2)-SAT.lever_arm*sin(y(9))))^2+(SAT.altitude-vehicle_height)^2);
                H(i,1)=-(XSAT(k,1+2*(i-1))-(y(1)-SAT.lever_arm*cos(y(9))))/d1;
                H(i,2)=-(XSAT(k,2+2*(i-1))-(y(2)-SAT.lever_arm*sin(y(9))))/d1;
                H(i,9)=H(i,1)*SAT.lever_arm*sin(y(9))-H(i,2)*SAT.lever_arm*cos(y(9));
                Y_measured(i,1)=dSAT(k,i);
                Y_Xminus(i,1)=d1;
            else
%                 d1=sqrt((XSAT(k,1+2*(i-1-size(XSAT,2)/2))-y(1))^2+...
%                     (XSAT(k,2+2*(i-1-size(XSAT,2)/2))-y(2))^2+(SAT.altitude-vehicle_height)^2);
                d1=sqrt((XSAT(k,1+2*(i-1-size(XSAT,2)/2))-(y(1)-SAT.lever_arm*cos(y(9))))^2+...
                    (XSAT(k,2+2*(i-1-size(XSAT,2)/2))-(y(2)-SAT.lever_arm*sin(y(9))))^2+(SAT.altitude-vehicle_height)^2);
                
%                 H(i,1)=-(y(3)-VSAT(k,1+2*(i-1-size(XSAT,2)/2)))/d1;
%                 H(i,2)=-(y(4)-VSAT(k,2+2*(i-1-size(XSAT,2)/2)))/d1;
%                 H(i,3)=(XSAT(k,1+2*(i-1-size(XSAT,2)/2))-y(1))/d1;
%                 H(i,4)=(XSAT(k,2+2*(i-1-size(XSAT,2)/2))-y(2))/d1;
                H(i,1)=-((y(3)+SAT.lever_arm*sin(y(9))*dy(9,1))-VSAT(k,1+2*(i-1-size(XSAT,2)/2)))/d1;
                H(i,2)=-((y(4)-SAT.lever_arm*cos(y(9))*dy(9,1))-VSAT(k,2+2*(i-1-size(XSAT,2)/2)))/d1;
                H(i,3)=(XSAT(k,1+2*(i-1-size(XSAT,2)/2))-(y(1)-SAT.lever_arm*cos(y(9))))/d1;
                H(i,4)=(XSAT(k,2+2*(i-1-size(XSAT,2)/2))-(y(2)-SAT.lever_arm*sin(y(9))))/d1;
                H(i,9)=H(i,1)*SAT.lever_arm*sin(y(9))-H(i,2)*SAT.lever_arm*cos(y(9));
                H(i,9)=H(i,9)+H(i,3)*SAT.lever_arm*cos(y(9))*dy(9,1)+H(i,4)*SAT.lever_arm*sin(y(9))*dy(9,1);
                
                Y_measured(i,1)=doppler(k,i-size(XSAT,2)/2);
                %n_vect=([XSAT(k,1+2*(i-1-size(XSAT,2)/2))-y(1) XSAT(k,2+2*(i-1-size(XSAT,2)/2))-y(2)]);
                n_vect=([XSAT(k,1+2*(i-1-size(XSAT,2)/2))-(y(1)-SAT.lever_arm*cos(y(9))) XSAT(k,2+2*(i-1-size(XSAT,2)/2))-(y(2)-SAT.lever_arm*sin(y(9)))]);
                n_vect=n_vect/norm(n_vect,2);
                %Y_Xminus(i,1)=dot_mod(y(3:4),n_vect)-dot_mod(([VSAT(k,1+2*(i-1-size(XSAT,2)/2)) VSAT(k,2+2*(i-1-size(XSAT,2)/2))]),n_vect);
                Y_Xminus(i,1)=dot_mod(y(3:4)'-([VSAT(k,1+2*(i-1-size(XSAT,2)/2))-SAT.lever_arm*sin(y(9))*dy(9,1)...
                    VSAT(k,2+2*(i-1-size(XSAT,2)/2))+SAT.lever_arm*cos(y(9))*dy(9,1)]),n_vect);
            end
        end
        %Calculos de inovacion de Kalman
        S=H*P_minus*H'+R;
        Sinv=inv(S);
        K=P_minus*H'*Sinv;
        P=(eye(length(dy))-K*H)*P_minus;
        %(eye(length(dy))-K*H)
        %P_minus(7,7)
        %P(7,7)
        
        nu=Y_measured-Y_Xminus;
        residues=(nu)'*(Sinv)*(nu);
        inovation=K*(nu);
        if residues > IMU_SIMULATION.threshold
            %[check,P,residues,inovation] = satelite_check(P_minus,H,R,SAT,IMU_SIMULATION,nu);
            [check,P,residues,inovation] = double_satelite_check(P_minus,H,R,SAT,IMU_SIMULATION,nu);
            if check < 0
                inovation=K*(nu);
                residues=(nu)'*(Sinv)*(nu);
            end
            IMU_SIMULATION.discarted_measures=IMU_SIMULATION.discarted_measures+1;
            IMU_SIMULATION.checked_measures(k)=check;
        end
        dy=dy+inovation;
    else
        residues=0;
        P=P_minus;
    end
end

end

function [result] = dot_mod(A,B)
result=A(1)*B(1)+A(2)*B(2);
end


function [check,P,residues,inovation] = satelite_check(P_minus,H,R,SAT,IMU_SIMULATION,nu)
XSAT=SAT.position;
n_satelites=size(XSAT,2)/2;
check=-1;
residues(n_satelites)=0;
for i=1:n_satelites
    
    Haux=H;
    Haux(i+n_satelites,:)=[];
    Haux(i,:)=[];
    
    Raux=R;
    Raux(i+n_satelites,:)=[];
    Raux(:,i+n_satelites)=[];
    Raux(i,:)=[];
    Raux(:,i)=[];
    
    nuaux=nu;
    nuaux(i+n_satelites,:)=[];
    nuaux(i,:)=[];
    
    S=Haux*P_minus*Haux'+Raux;
    Sinv=inv(S);
    K=P_minus*Haux'*Sinv;
    P=(eye(size(P_minus,1))-K*Haux)*P_minus;
    
    residues(i)=(nuaux)'*(Sinv)*(nuaux);
    inovation=K*(nuaux);
    if residues(i) < IMU_SIMULATION.threshold
        check=i;
        break
    end
end
residues=residues(i);
end

function [check,P,residues,inovation] = double_satelite_check(P_minus,H,R,SAT,IMU_SIMULATION,nu)
XSAT=SAT.position;
n_satelites=size(XSAT,2)/2;
check=-1;
residues(n_satelites)=0;
for i=1:n_satelites
    
    Haux=H;
    Haux(i+n_satelites,:)=[];
    Haux(i,:)=[];
    
    Raux=R;
    Raux(i+n_satelites,:)=[];
    Raux(:,i+n_satelites)=[];
    Raux(i,:)=[];
    Raux(:,i)=[];
    
    nuaux=nu;
    nuaux(i+n_satelites,:)=[];
    nuaux(i,:)=[];
    
    S=Haux*P_minus*Haux'+Raux;
    Sinv=inv(S);
    K=P_minus*Haux'*Sinv;
    P=(eye(size(P_minus,1))-K*Haux)*P_minus;
    
    residues(i)=(nuaux)'*(Sinv)*(nuaux);
    inovation=K*(nuaux);
    if residues(i) < IMU_SIMULATION.threshold
        check=i;
        break
    end
end
if check == -1
    for i=1:n_satelites
        for j=i+1:n_satelites
            Haux=H;
            Haux(j+n_satelites,:)=[];
            Haux(i+n_satelites,:)=[];
            Haux(j,:)=[];
            Haux(i,:)=[];
            
            
            Raux=R;
            Raux(j+n_satelites,:)=[];
            Raux(i+n_satelites,:)=[];
            Raux(:,j+n_satelites)=[];
            Raux(:,i+n_satelites)=[];
            Raux(j,:)=[];
            Raux(i,:)=[];
            Raux(:,j)=[];
            Raux(:,i)=[];
            
            nuaux=nu;
            nuaux(j+n_satelites,:)=[];
            nuaux(i+n_satelites,:)=[];
            nuaux(j,:)=[];
            nuaux(i,:)=[];
            
            S=Haux*P_minus*Haux'+Raux;
            Sinv=inv(S);
            K=P_minus*Haux'*Sinv;
            P=(eye(size(P_minus,1))-K*Haux)*P_minus;
            
            residues(i)=(nuaux)'*(Sinv)*(nuaux);
            inovation=K*(nuaux);
            if residues(i) < IMU_SIMULATION.threshold
                check=i*10+j;
                break
            end
        end
    end
end
residues=residues(i);
end