function [dy] = mechanization_IMUonly(k,y,IMU_SIMULATION)
%Esta funcion genera la derivada temporal del vector de estado para
%inyectar en el esquema temporal
%INPUT
%k: Instante de tiempo (indice)
%y: Vector de estado
%AM: Medidas de los acelerometros en x e y, y velocidad angular z
%OUTPUT
%dy: Derivada temporal del vector de estado
AM=IMU_SIMULATION.measures;

%Transformacion de ejes cuerpo a ejes fijos
AMX=(AM(k,1)-y(5)-y(7))*cos(y(9))-(AM(k,2)-y(6)-y(8))*sin(y(9));
AMY=(AM(k,1)-y(5)-y(7))*sin(y(9))+(AM(k,2)-y(6)-y(8))*cos(y(9));
AMW=AM(k,3)-y(10)-y(11);
%Estimacion de la derivada
dy(11,1)=0;
dy(1,1)=y(3);dy(2,1)=y(4);
dy(3,1)=(AMX);dy(4,1)=(AMY);
dy(5,1)=0;dy(6,1)=0;
dy(7,1)=0;dy(8,1)=0;
dy(9,1)=(AMW);dy(10,1)=0;dy(11,1)=0;
end

