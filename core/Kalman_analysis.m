function [guess_container] = Kalman_analysis(SIMULATION,tiempo)
%Esta funcion realiza calculos simulados en tiempo real, es decir, los
%calculos que deberia realizar el vehiculo a bordo con los datos
%suministrados por el filtro de kalman.
%Por lo tanto los datos deben ser analizados en bloques de tiempo desde el
%inicio hasta el final pese a tener todos los datos en todos los instantes
%de tiempo.
%INPUT
%XKalman: Vector de estado en todos los tiempos de Kalman
%residuos: Los residuos generados en la propagacion
%imu_indicator: La diferencia de distancias entre la propagacion con Kalman
%y una puramente inercial
%tiempo: objeto tiempo que contiene parametros temporales
%SATPARAMS: Parametros principales del satelite, como el tiempo de spoofing
%OUTPUT
%Los diferentes tiempo predichos por los metodos desarrollados a
%continuacion. Deberian de ser lo mas parecidos al tiempo de inicio de
%Spoofing

dt=tiempo.dt;
dtGNSS=tiempo.dtGNSS;
T_alineamiento=tiempo.alineamiento;

XKalman=SIMULATION.state_vector;
residuos=SIMULATION.residuos;
P_buff=SIMULATION.PMatrix;
K_buff=SIMULATION.KMatrix;
imu_indicator=SIMULATION.imu_only;
dXK=SIMULATION.state_vector_derivate;
%AM=SIMULATION.measures;

T=size(XKalman,2);t=0:dt:T*dt-dt;
XK=XKalman(1,:);YK=XKalman(2,:);
biaskx0=XKalman(5,:);biasky0=XKalman(6,:);
biaskx1=XKalman(7,:);biasky1=XKalman(8,:);
biaskw0=XKalman(10,:);biaskw1=XKalman(11,:);

SATPARAMS=SIMULATION.SAT;
n_satelites=SATPARAMS.n_satelites;
n_satelites_spoof=SATPARAMS.n_satelites_spoof;
coeff_spoof0=SATPARAMS.coeff_spoof_0;
coeff_spoof1=SATPARAMS.coeff_spoof_1;
sigma_GNSS=SATPARAMS.sigma_GNSS;
t_begin_spoof=SATPARAMS.t_begin_spoof;
t_end_spoof=SATPARAMS.t_end_spoof;
t_begin_false_alarm=SATPARAMS.false_alarm_start;
t_end_false_alarm=SATPARAMS.false_alarm_end;
bias_param=SIMULATION.IMU.params;
AM=SIMULATION.IMU.measures;



fprintf('********************************************************\n');
fprintf('Kalman Analysis...\n');

%% ---------------  PLOTEO -------------

figure(450)
close 450
figure(450)

subplot(3,2,1),
hold on
aux(T,1,1)=0;
aux(:,1,1)=P_buff(5,5,:);
taux=t;
plot(taux,aux),axis tight, grid minor;
plot(linspace(min(taux),max(taux),length(taux)),zeros(1,length(taux)));
plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
title('P 5')

subplot(3,2,2),
hold on
aux(T,1,1)=0;
aux(:,1,1)=P_buff(6,6,:);
taux=t;
plot(taux,aux),axis tight, grid minor;
plot(linspace(min(taux),max(taux),length(taux)),zeros(1,length(taux)));
plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
title('P 6')

subplot(3,2,3),
hold on
aux(T,1,1)=0;
aux(:,1,1)=P_buff(7,7,:);
taux=t;
plot(taux,aux),axis tight, grid minor;
plot(linspace(min(taux),max(taux),length(taux)),zeros(1,length(taux)));
plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
title('P 7')

subplot(3,2,4),
hold on
aux(T,1,1)=0;
aux(:,1,1)=P_buff(8,8,:);
taux=t;
plot(taux,aux),axis tight, grid minor;
plot(linspace(min(taux),max(taux),length(taux)),zeros(1,length(taux)));
plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
title('P 8')


subplot(3,2,6),
hold on
aux(T,1,1)=0;
aux(:,1,1)=P_buff(10,10,:);
taux=t;
plot(taux,aux),axis tight, grid minor;
plot(linspace(min(taux),max(taux),length(taux)),zeros(1,length(taux)));
plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
title('P 10')

subplot(3,2,5),
hold on
aux(T,1,1)=0;
aux(:,1,1)=P_buff(9,9,:);
taux=t;
plot(taux,aux),axis tight, grid minor;
plot(linspace(min(taux),max(taux),length(taux)),zeros(1,length(taux)));
plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
title('P 9')

figure(451)
close 451
figure(451)

subplot(2,2,1),
hold on
aux(T,1,1)=0;
aux(:,1,1)=P_buff(1,1,:);
taux=t;
plot(taux,aux),axis tight, grid minor;
plot(linspace(min(taux),max(taux),length(taux)),zeros(1,length(taux)));
plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
title('P 1')

subplot(2,2,2),
hold on
aux(T,1,1)=0;
aux(:,1,1)=P_buff(2,2,:);
taux=t;
plot(taux,aux),axis tight, grid minor;
plot(linspace(min(taux),max(taux),length(taux)),zeros(1,length(taux)));
plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
title('P 2')

subplot(2,2,3),
hold on
aux(T,1,1)=0;
aux(:,1,1)=P_buff(3,3,:);
taux=t;
plot(taux,aux),axis tight, grid minor;
plot(linspace(min(taux),max(taux),length(taux)),zeros(1,length(taux)));
plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
title('P 3')

subplot(2,2,4),
hold on
aux(T,1,1)=0;
aux(:,1,1)=P_buff(4,4,:);
taux=t;
plot(taux,aux),axis tight, grid minor;
plot(linspace(min(taux),max(taux),length(taux)),zeros(1,length(taux)));
plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(min(aux(10:end)),max(aux(10:end)),100)), grid on, axis tight;
title('P 4')



% figure(452)
% close 452
% figure(452)
% % P_buff(1,1,:)=P_buff(1,1,:)/100;
% % P_buff(2,2,:)=P_buff(2,2,:)/100;
% for i=1:floor(T/200):T
%    surf(P_buff(5:11,5:11,i)), grid on; 
%    %axis([0 11 0 11 -1 1])
%    pause(0.01);
% end



end

function [guess,i_guess,ema_long,ema_short,macd]=macd_calculate(data_intput,t_input,dtGNSS,period_long,period_short,minhash)
[ema_long] = ema_fun(data_intput,period_long);
[ema_short] = ema_fun(data_intput,period_short);
macd=0;macd(length(ema_long))=0;
for i=1:length(ema_long)
    macd(i)=(ema_short(i)-ema_long(i))/ema_long(i);
end
%hash=floor(length(ema_long)/10)
hash=floor(max([minhash dtGNSS*10])/(t_input(2)-t_input(1)));
guess=0;i_guess=1;
for i=1+hash:hash:length(ema_long)-hash
    min(macd(i:i+hash));
    if min(macd(i:i+hash)) > 0 && i/length(ema_long) > 0.1%max(macd(1:floor(length(ema_long)/100)))+0.3
        guess=t_input(i+hash);
        i_guess=i;
        break
    end
end

end

function [guess,i_guess]=macd_jump(macd,t_input,dtGNSS,sign_gap,GNSS_periods)
hash=floor(max([20 dtGNSS*GNSS_periods])/(t_input(2)-t_input(1)));
guess=0;i_guess=1;
for i=1+hash:hash:length(macd)-hash
    min(macd(i:i+hash));
    if sign_gap > 0
        if min(macd(i:i+hash)) > 0
            guess=t_input(i+hash);
            i_guess=i;
            break
        end
    else
        if max(macd(i:i+hash)) < 0
            guess=t_input(i+hash);
            i_guess=i;
            break
        end
    end
end
end

function [guess,i_guess]=tranquility_catch(input,t_input,dtGNSS,factor)
hash=floor((dtGNSS*10)/(t_input(2)-t_input(1)));
hash=4;
guess=0;i_guess=1;
for i=5*hash:length(input)-5*hash
    if max(input(10:i))*factor > max(input(i:i+hash))
        if min(input(10:i))*factor < min(input(i:i+hash))
            guess=t_input(i+hash);
            i_guess=i;
            break
        end
    end
end
end


function [data_input,t_input]=get_var(data_input,variance_factor,dt)
T=length(data_input);
[data_input,t_input] = variance(data_input,variance_factor,dt);

aux(T)=0;t_aux(T)=0;count=1;
for i = 1:length(data_input)
    if data_input(i) ~= 0
        aux(count)=data_input(i);
        t_aux(count)=t_input(i);
        count=count+1;
    end
end
data_input=aux(1:count-1);
if count ==1
    data_input=0;
end
t_input=t_aux(1:count-1);
end

function [guess]=get_jump(data_input,jump_factor)
guess=0;
threshold=mean(data_input(1:floor(length(data_input)/10)));
for i = floor(length(data_input)/10):length(data_input)
    if data_input(i) > jump_factor*threshold
        guess=i;
        break
    end
end
end

function [guess]=bias_jump(data_input,jump_gap)
guess=1;
for i = floor(length(data_input)/20):length(data_input)
    if abs(data_input(i)) > jump_gap
        guess=i;
        break
    end
end
end


function [variance,t_var] = variance(data,period,dt)

T=length(data);
kn=1;aux=0;variance=0;
count_var=1;count_i=1;
for i=floor(T/10):T
    aux(count_i)=data(i);
    count_i=count_i+1;
    kn=kn+1;
    if kn>period
        variance(count_var)=var(aux);
        aux=0;
        t_var(count_var)=i*dt;
        count_var=count_var+1;
        count_i=1;
        kn=1;
    end
end
end

function [ema] = ema_fun(data,period)

T=length(data);
cff=2/(period+1);
ema(T)=0;
for i=1:T
    if i>period
        ema(i)=(data(i)-ema(i-1))*cff+ema(i-1);
    elseif i==period
        ema(i)=mean(data(1:i));
    end
end
end

function [yf,nuf,paramf]=chiopt(t,residuos)

minim=1e6;nu(100)=0;diff(100,100)=0;yf=0;nuf=0;paramf=0;
for i=1:100
    nu(i)=0.5+i/10;
    clc
    param(100)=0;
    for j=1:100
        
        param(j)=2+j;
        count=1;y(length(t))=0;
        for k=1:length(t)
            if residuos(count)==0
                y(count)=0;
            else
                y(count) = chifun(t(k)/param(j),nu(i));
            end
            count=count+1;
        end
        y=y(1:count-1);
        %         figure(10)
        %         plot(t,residuos,t,y/max(y)),grid on, axis tight;
        %         pause(0.01)
        diff(i,j)=norm(residuos-y/max(y)*max(residuos),2);
        if diff(i,j)<minim
            minim=diff(i,j);
            nuf=nu(i);
            paramf=param(j);
            yf=y;
        end
    end
end
end


function [y]=chifun(x,nu)

y=zeros(length(x),1);
for i=1:length(x)
    y(i) = x(i)^((nu(i)-2)/2)*exp(-x(i)/2)/gamma(nu(i)/2)/2^(nu(i)/2);
end
end

function [C,t_corr]=correlation_coeficient(x,y,hash,dt)

hash=floor(hash);count=1;C(length(x))=0;t_corr(length(x))=0;
for i=1:hash:length(x)-hash
    x_samples=x(i:i+hash);
    y_samples=y(i:i+hash);
    x_tau=mean(x_samples);
    y_tau=mean(y_samples);
    aux_top=mean((x_samples-x_tau).*(y_samples-y_tau));
    aux_bottom=sqrt(mean((x_samples-x_tau).^2)*mean((y_samples-y_tau).^2));
    C(count)=aux_top/aux_bottom;
    t_corr(count)=i*dt;
    count=count+1;
end
C=abs(C(1:count-1));t_corr=t_corr(1:count-1);
end


function [prob_res,eje_x] = probability_distribution(data_input,time_begin,N_div)
%ORDENA LOS RESIDUOS EN DIFERENTES NIVELES DE CANTIDADES
aux=data_input/max(data_input);
aux=aux(ceil(time_begin):end);
%N_div=10;
prob_res(N_div)=0;div=linspace(1,0,N_div);
maxaux=max(aux);

for i=1:length(aux)
    
    for j=1:length(div)
        if aux(i) ~=0
            if aux(i)>(maxaux*div(j))
                prob_res(j)=prob_res(j)+1;
                break
            end
        end
    end
end
prob_res=(prob_res/length(aux));
eje_x=max(aux);
end