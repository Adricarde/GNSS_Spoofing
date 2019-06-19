function [guess_container] = postprocess(SIMULATION,tiempo)
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
residues=SIMULATION.residuos;
P_buff=SIMULATION.PMatrix;
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
fprintf('Postprocess...\n');
%% --------------      METODO DE VARIANZA DE BIAS  -------------
fprintf('Bias... ');

biaskx=biaskx0+biaskx1;
biasky=biasky0+biasky1;
biaskw=biaskw0+biaskw1;
[i_bias_x]=bias_jump(biaskx,bias_param.sigma_bias_a*20);
[i_bias_y]=bias_jump(biasky,bias_param.sigma_bias_a*20);
[i_bias_w]=bias_jump(biaskw,bias_param.sigma_bias_w*20);
if  biaskx(i_bias_x) > biasky(i_bias_y)
    bias_guess=(i_bias_x)*dt;
else
    bias_guess=(i_bias_y)*dt;
end

%% --------------      METODO DE VARIANZA DE RESIDUOS  -------------
fprintf('Kalman Residues... ');
[residues_var,t_varr]=get_var(residues,20,dt);

[res_guess,i_res_guess,ema_long_res,ema_short_res,macd_res]=macd_calculate(residues_var,t_varr,5.5*(t_varr(2)-t_varr(1)),20,15,20);
[res_guess,i_res_guess,upres,downres]=tranquility_catch(diff(macd_res),t_varr(1:end-1),dtGNSS,0.15);

%% --------------      METODO DE CHI CUADRADO  -------------
fprintf('Chi Squared...\n');
%ORDENA LOS RESIDUOS EN DIFERENTES NIVELES DE CANTIDADES
count=1;residuos2(T)=0;t2(T)=0;
for i=2*T_alineamiento:T
    if residues(i) > 0
        residuos2(count)=residues(i);
        t2(count)=dt*i;
        count=count+1;
    end
end
residuos2=residuos2(1:count-1);
t2=t2(1:count-1);
[ema_residuos] = ema_fun(residuos2,3);

aux=ema_residuos;
N_div=100;
prob_res(N_div)=0;div=linspace(1,0,N_div);
maxaux=max(aux);
count_res=1;
for k=1:length(aux)
    for j=1:length(div)
        if aux(k) ~=0
            if aux(k)>(maxaux*div(j))
                prob_res(j)=prob_res(j)+1;
                count_res=count_res+1;
                break
            end
        end
    end
end
prob_res=prob_res/count_res;

[spoofing_threshold] = probability_interval(n_satelites,0.0001);
spoofing_threshold=spoofing_threshold*10;
chi_quad_guess=1;
for i=1:length(t2)
    if ema_residuos(i) > spoofing_threshold
        chi_quad_guess=t2(i);
        break
    end
end

%% --------------      METODO DE DOBLE PROPAGACIÓN  -------------
fprintf('Double Propagation... ');

%[imu_var,t_varimu]=get_var(imu_indicator,10,dt);
%[imu_guess,i_imu_guess,ema_long_imu,ema_short_imu,macd_imu]=macd_calculate(imu_var,t_varimu,3*(t_varimu(2)-t_varimu(1)),10,5,10*(t_varimu(2)-t_varimu(1)));
%[imu_guess,i_imu_guess]=tranquility_catch(diff(macd_imu),t_varimu(1:end-1),dtGNSS,0.1);

[i_imu_guess]=bias_jump(imu_indicator,100*max(imu_indicator(1:floor(20/dt))));
imu_guess=dt*i_imu_guess;
%% --------------      METODO DE ELIPSES  -------------
fprintf('Ellipses... ');

if SATPARAMS.spoof_flag==1
    xGNSS=SATPARAMS.vehicle_position.x_spoofed;
    yGNSS=SATPARAMS.vehicle_position.y_spoofed;
else
    xGNSS=SATPARAMS.vehicle_position.x;
    yGNSS=SATPARAMS.vehicle_position.y;
end
XKelip=interp1(t,XK,SATPARAMS.tGNSS,'PCHIP')';
YKelip=interp1(t,YK,SATPARAMS.tGNSS,'PCHIP')';
d=[xGNSS yGNSS]-[XKelip YKelip];
thetaKelip=XKalman(9,:);
a=100;b=80;eps=sqrt(1-b*b/a/a);
thetab(length(XKelip))=0;d_vect(length(XKelip))=0;
count2=1;elip_vect(length(XKelip))=0;relip(length(XKelip))=0;
elip_count(length(XKelip))=0;

for i = 1:length(XKelip)
    thetab(i)=atan2(d(i,2),d(i,1))-thetaKelip(i);
    relip(i)=b/sqrt(1-(eps*cos(thetab(i)))^2);
    d_vect(i)=norm(d(i,:),2);
    if d_vect(i) > relip(i)
        count2=count2+1;
    end
    elip_vect(i)=d_vect(i);
    elip_count(i)=count2;
end
% figure(700)
% subplot(2,1,1),
% hold on
% plot(SATPARAMS.tGNSS,elip_count),axis tight, grid on;
% subplot(2,1,2),
% hold on
% plot(SATPARAMS.tGNSS(1:end-1),diff(elip_count)),axis tight, grid on;
%[elip_var,t_varelip]=get_var(elip_vect,2,SATPARAMS.tGNSS(2)-SATPARAMS.tGNSS(1));
%[distance_guess,i_distance_guess,ema_long_elip,ema_short_elip,macd_elip]=macd_calculate(elip_var,t_varelip,0.8*(t_varelip(2)-t_varelip(1)),6,4,2*(t_varelip(2)-t_varelip(1)));
%[distance_guess,i_distance_guess]=tranquility_catch(diff(macd_elip),t_varelip(1:end-1),dtGNSS,0.15);
[i_distance_guess,elip_thresh]=elip_jump(elip_vect);
distance_guess=SATPARAMS.tGNSS(i_distance_guess);
[elip_guess,i_distance_guess]=nonzero_catch(diff(elip_count),SATPARAMS.tGNSS(1:end-1));
%% --------------      METODO VARIANZA DE ALLAN  -------------

% N_allan=floor(T/8);
% aux_allan=biaskx0(N_allan:2*N_allan)+biaskx1(N_allan:2*N_allan);
%
% count=1;sigma_allan(N_allan)=0;tau_allan(N_allan)=0;
% for i=5000:1000:floor(40000)
%     n_allan=floor(i);
%     sum_allan=0;
%     for k=1:N_allan-n_allan*2
%         omega_allan0 = mean(aux_allan(k:k+n_allan));
%         omega_allan1 = mean(aux_allan(k+n_allan:k+2*n_allan));
%         sum_allan = sum_allan + (omega_allan1-omega_allan0)^2;
%     end
%     sigma_allan(count)=sum_allan/2/(N_allan-n_allan*2);
%     tau_allan(count)=dt*n_allan;
%     count=count+1;
%     i
% end
% sigma_allan=sigma_allan(1:count-1);
% tau_allan=tau_allan(1:count-1);
% figure(599)
% loglog(tau_allan,sigma_allan),axis tight, grid on;
%plot(log10(tau_allan),log10(sigma_allan)),axis tight, grid on;

%% --------------      METODO DE NÚMERO DE SATÉLITES DESCARTADOS  -------------
fprintf('Satelites Discard... ');

aux_matrix(T,n_satelites)=0;
t_aux(T,1)=0;count=1;
detect_time=floor(10/dtGNSS);
discard_guess=0;

for i=1:T
    if SIMULATION.IMU.checked_measures(i) == -1
        count=count+1;
    elseif SIMULATION.IMU.checked_measures(i) > 0
        count=1;
    end
    if count > detect_time
        discard_guess=i*dt;
        break
    end
end


%% ------------------------      RESULTADOS -----------------------------
fprintf('Done\n');

guess_container.bias=bias_guess;guess_container.res=res_guess;
guess_container.imu=imu_guess;guess_container.chi_quad=chi_quad_guess;
guess_container.gnss=elip_guess;guess_container.elip=distance_guess;
guess_container.discard=discard_guess;


guess_vect=[bias_guess res_guess imu_guess chi_quad_guess distance_guess elip_guess discard_guess];
fail_vect(length(guess_vect))=0;
for i=1:length(guess_vect)
    if guess_vect(i)-t_begin_spoof*tiempo.dt<0
        fail_vect(i)=1;
    end
    if guess_vect(i)<2
        guess_vect(i)=0;
    end
end
guess_container.vect=guess_vect;
guess_container.fail=fail_vect;

%% ------------------------  PLOTEO ---------------------------
if SIMULATION.plot_postprocess
    
    if SIMULATION.SAT.spoof_flag==1
        figure(399)
        %         close 400
        %         figure(400)
        subplot(4,1,1:2),
        
        hold on
        bar(guess_vect,'FaceColor',[.6 .2 .2],'EdgeColor',[.8 .2 .2],'LineWidth',1.5);
        bar(t_begin_spoof*tiempo.dt*ones(7,1),'FaceColor',[.2 .6 .2],'EdgeColor',[.2 .8 .2],'LineWidth',1.5),grid minor, axis tight;
        for i=1:length(guess_vect)
            if fail_vect(i)==1
                bar(i,t_begin_spoof*tiempo.dt,'FaceColor',[.2 .2 .6],'EdgeColor',[.2 .2 .8],'LineWidth',1.5);
            end
        end
        xticks(linspace(1,7,7))
        xticklabels({'Bias','Residues','Double Prop','Chi','Distance','Ellipses','Discard'})
        ylabel('Simulation Time (s)');
    else
        figure(399)
        close 399
        figure(399)
        
    end
    %     subplot(3,1,1:2),
    
    
    if SIMULATION.SAT.spoof_flag==1
        ha=subplot(4,1,3);
    else
        ha=subplot(4,1,4);
    end
    
    pos = get(ha,'Position');
    un = get(ha,'Units');
    delete(ha)
    data=[{'Bias','Residues','Double Prop','Chi','Distance','Ellipses','Discard'}' num2cell(guess_vect)']';
    uitable('Data',data,'Units',un,'Position',pos);
    
    
    figure(411)
    close 411
    figure(411)
    
    %subplot(2,2,1),
    
    plot(t,biaskx,t,biasky,t,biaskw),grid minor, axis tight;
    hold on
    axis manual
    plot(linspace(min(t),max(t),length(t)),zeros(1,length(t)));
    plot(linspace(min(t),max(t),length(t)),bias_param.sigma_bias_a*50*ones(1,length(t)));
    plot(linspace(min(t),max(t),length(t)),-bias_param.sigma_bias_a*50*ones(1,length(t)));
    
    plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(-10,10,100));
    plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(-10,10,100));
    xlabel('Time (s)');
    ylabel('Bias Value');
    title('BIASX, BIASY and BIASW')
    set(gca, 'FontName', 'Cambria','Fontsize',20)
    
    figure(412)
    close 412
    figure(412)
    %subplot(2,2,2),
    plot(t_varr(1:end-1),diff(macd_res)),grid minor, axis tight;
    hold on
    axis manual
    plot(linspace(min(t_varr),max(t_varr),length(t_varr)),zeros(1,length(t_varr)));
    plot(linspace(min(t),max(t),length(t)),upres*ones(1,length(t)));
    plot(linspace(min(t),max(t),length(t)),downres*ones(1,length(t)));
    plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(-10,10,100));
    plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(-10,10,100));
    plot(t_end_false_alarm*dt*ones(1,100),1.1*linspace(-10,10,100));
    legend('Indicator','Zero','Upper Threshold','Lower Threshold','Spoofing Begin Time','FA Begin Time','FA End Time');
    xlabel('Time (s)');
    ylabel('MACD Derivate')
    title('Residues Derivate Variance')
    set(gca, 'FontName', 'Cambria','Fontsize',20)
    
    figure(415)
    close 415
    figure(415)
    subplot(3,1,1:2),
    plot(t_varr,(residues_var)),grid minor, axis tight;
    hold on
    axis manual
    plot(t_varr,(ema_long_res)),grid minor, axis tight;
    plot(t_varr,(ema_short_res)),grid minor, axis tight;
    plot(linspace(min(t_varr),max(t_varr),length(t_varr)),zeros(1,length(t_varr)));
    legend('Residues Variance','Long EMA','Short EMA');
    title('Residues Variance and EMAs')
    set(gca, 'FontName', 'Cambria','Fontsize',20)
    
    subplot(3,1,3),
    title('Residues MACD')
    plot(t_varr,(macd_res)),grid minor, axis tight;
    legend('Residues PPO');
    xlabel('Time (s)');
    set(gca, 'FontName', 'Cambria','Fontsize',20)
    
    
    figure(413)
    close 413
    figure(413)
    taux=t;
    aux=imu_indicator;
    plot(taux,aux),grid minor, axis tight;
    hold on
    axis manual
    plot(linspace(min(t),max(t),length(t)),100*max(imu_indicator(1:floor(20/dt)))*ones(1,length(t)));
    plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(-10,10,100));
    plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(-10,10,100));
    plot(t_end_false_alarm*dt*ones(1,100),1.1*linspace(-10,10,100));
    legend('Indicator','Threshold','Spoofing Begin Time','FA Begin Time','FA End Time');
    xlabel('Time (s)');
    ylabel('DP Indicator');
    title('Double Propagation Indicator')
    set(gca, 'FontName', 'Cambria','Fontsize',20)
    
    if length(elip_vect) > 15
        
        figure(414)
        close 414
        figure(414)
        
        plot(SATPARAMS.tGNSS,elip_vect),grid minor, axis tight;
        hold on
        axis manual
        plot(linspace(min(t),max(t),length(t)),elip_thresh*ones(1,length(t)));
        plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(-max(elip_vect),max(elip_vect),100));
        plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(-max(elip_vect),max(elip_vect),100));
        plot(t_end_false_alarm*dt*ones(1,100),1.1*linspace(-max(elip_vect),max(elip_vect),100));
        legend('Indicator','Threshold','Spoofing Begin Time','FA Begin Time','FA End Time');
        xlabel('Time (s)');
        ylabel('Distance to Vehicle');
        title('Ellipses Distance Indicator')
        set(gca, 'FontName', 'Cambria','Fontsize',20)
        
        
        figure(416)
        close 416
        figure(416)
        
        subplot(2,1,1),
        plot(SATPARAMS.tGNSS,elip_count),grid minor, axis tight;
        hold on
        axis manual
        plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(-max(elip_count),max(elip_count),100));
        plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(-max(elip_count),max(elip_count),100));
        plot(t_end_false_alarm*dt*ones(1,100),1.1*linspace(-max(elip_count),max(elip_count),100));
        xlabel('Time (s)');
        title('Ellipses Count');
        legend('Count','Spoofing Begin Time','FA Begin Time','FA End Time');
        set(gca, 'FontName', 'Cambria','Fontsize',20)
        
        subplot(2,1,2),
        plot(SATPARAMS.tGNSS(1:end-1),diff(elip_count)),grid minor, axis tight;
        hold on
        axis manual
        plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(-max(elip_count),max(elip_count),100));
        plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(-max(elip_count),max(elip_count),100));
        plot(t_end_false_alarm*dt*ones(1,100),1.1*linspace(-10,10,100));
        axis([min(t) max(t) -1 2])
        xlabel('Time (s)');
        legend('Indicator','Spoofing Begin Time','FA Begin Time','FA End Time');
        title('Ellipses Count Derivate Indicator');
        set(gca, 'FontName', 'Cambria','Fontsize',20)
        
    end
    
    
    if SIMULATION.SAT.spoof_flag==0
        figure(402)
        close 402
        figure(402)
        %         subplot(2,1,1),
        x=linspace(0,15,N_div);
        y=chifun(x,n_satelites);
        hold on
        area(maxaux*div,prob_res/max(prob_res)*max(y),'FaceAlpha',0.2);
        area(maxaux*div,prob_res,'FaceAlpha',0.2);
        area(x,y,'FaceAlpha',0.2), axis tight, grid minor;
        legend('Scaled','True','N_{satelites}');
        title('Chi Squared Method')
        xlabel('Residue Magnitude')
        ylabel('Probability')
        
        figure(403)
        close 403
        figure(403)
        x=linspace(0,15,N_div);
        y2=chifun(x,2);
        y3=chifun(x,3);
        y4=chifun(x,4);
        y5=chifun(x,5);
        y6=chifun(x,6);
        hold on
        plot(x,y2,x,y3,x,y4,x,y5,x,y6), axis tight, grid minor;
        legend('2 Satellites','3 Satellites','4 Satellites','5 Satellites','6 Satellites');
        xlabel('Residue Magnitude')
        ylabel('Probability')
        set(gca, 'FontName', 'Cambria','Fontsize',20)
        
    end
    
end
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

function [guess,i_guess,upbound,downbound]=tranquility_catch(input,t_input,dtGNSS,factor)
hash=floor((dtGNSS*30)/(t_input(2)-t_input(1)));
hash=floor(15/(t_input(2)-t_input(1)));
%hash=7;
guess=0;i_guess=1;
for i=2*hash:length(input)-2*hash
    if max((input(50:i)))*factor > max(input(i:i+hash))
        if min(input(50:i))*factor < min(input(i:i+hash))
            guess=t_input(i+hash);
            i_guess=i;
            break
        end
    end
end
upbound=max((input(50:i)))*factor;
downbound=min(input(50:i))*factor;
end

function [guess,i_guess]=nonzero_catch(input,t_input)
hash=floor(20/(t_input(2)-t_input(1)));
%hash=7;
guess=0;i_guess=1;
for i=20*hash:length(input)-5*hash
    if  mean(input(i:i+hash)) == 1
        guess=t_input(i+hash);
        i_guess=i;
        break
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

function [guess,threshold]=elip_jump(data_input)
guess=1;
threshold=15*mean(data_input(1:floor(length(data_input)/20)));
for i = floor(length(data_input)/20):length(data_input)
    if abs(data_input(i)) > threshold
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


function [y]=chisum(x,nu)
%
% y=zeros(length(x),1);
% for i=1:length(x)
%      y(i) = gammainc(x(i)/2,nu/2,'lower')/gamma(nu/2);
%     %y(i) = 1-gammainc(x(i)/2,nu/2,'upper')/gamma(nu/2);
% end
y=chifun(x,nu);
y=cumsum(y);
y=y/max(y);
end

function [threshold]=probability_interval(n_satelites,detection_probability)
%
N_div=1000;
x=linspace(0,20,N_div);
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