function plot_sim(WORLDMODEL,tiempo,SIMULATION,SIMULATION_s)
%% -------------- CALCULOS PREVIOS  -----------------
XReal=WORLDMODEL.true_tray;
XR=XReal(:,1);YR=XReal(:,2);THR=XReal(:,7);
AXR=XReal(:,5);AYR=XReal(:,6);WTHR=XReal(:,8);

dt=tiempo.dt;
dtGNSS=tiempo.dtGNSS;
T_alineamiento=tiempo.alineamiento;
T=size(XReal,1);t=0:dt:T*dt-dt;

%IMU
IMU_SIMULATION=SIMULATION.IMU;

AM=IMU_SIMULATION.measures;
bias_param=IMU_SIMULATION.params;
biasx=IMU_SIMULATION.biasx;
biasy=IMU_SIMULATION.biasy;
biasw=IMU_SIMULATION.biasw;

XF=SIMULATION.state_vector;
XGNSS=SIMULATION.gnss_state_vector;
dXK=SIMULATION.state_vector_derivate;
residuos=SIMULATION.residuos;

XFS=SIMULATION_s.state_vector;
XGNSS_S=SIMULATION_s.gnss_state_vector;
dXKs=SIMULATION_s.state_vector_derivate;
residuos_spoof=SIMULATION_s.residuos;

radio=max(XR);
XK=XF(1,:);YK=XF(2,:);
VXK=XF(3,:);VYK=XF(4,:);
biaskx0=XF(5,:);biasky0=XF(6,:);
biaskx1=XF(7,:);biasky1=XF(8,:);
biaskw0=XF(10,:);biaskw1=XF(11,:);
THK=XF(9,:);

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
XSAT=SATPARAMS.position;

XKS=XFS(1,:);YKS=XFS(2,:);
biaskx0_s=XFS(5,:);biasky0_s=XFS(6,:);
biaskx1_s=XFS(7,:);biasky1_s=XFS(8,:);
biaskw0_s=XFS(10,:);biaskw1_s=XFS(11,:);
THKS=XFS(9,:);WKS=XFS(10,:);

N_div=10;div=linspace(1,0,N_div);

%% -------  PLOTEO -------------
n_fig=1;
% figure(n_fig)
% close(n_fig)
% figure(n_fig)
figure('Name','Trajectories','NumberTitle','off','units','normalized','outerposition',[0 0 1 1],'Color',[234 253 234]/255)


% subplot(6,1,1:4),
hold on
plot(XR/1000,YR/1000,'k-');
plot(XK/1000,YK/1000,'b-',XK(t_begin_spoof:end)/1000,YK(t_begin_spoof:end)/1000,'g-',XKS(t_begin_spoof:end)/1000,YKS(t_begin_spoof:end)/1000,'r-.'), grid on;%, axis([1.2*min(XKS) 1.2*max(XKS) 1.2*min(YKS) 1.2*max(YKS)]);
plot(XGNSS(1,:)/1000,XGNSS(2,:)/1000,'c--',XGNSS_S(1,:)/1000,XGNSS_S(2,:)/1000,'m--'), grid on;%, axis([1.2*min(XKS) 1.2*max(XKS) 1.2*min(YKS) 1.2*max(YKS)]);
n_sat=size(XSAT,2)/2;
for i=1:n_sat
    if i <= n_satelites_spoof
        colorp=[153, 0, 0]/255;
    else
        colorp=[0, 153, 76]/255;
    end
    plot(XSAT(:,1+2*(i-1))/1000,XSAT(:,2+2*(i-1))/1000,'color',colorp,'LineStyle','-','Marker','.'), grid on;%, axis([1.2*min(XKS) 1.2*max(XKS) 1.2*min(YKS) 1.2*max(YKS)]);
end
%scatter(SATPARAMS.pos_spoof(1),SATPARAMS.pos_spoof(2)), grid on;%, axis([1.2*min(XKS) 1.2*max(XKS) 1.2*min(YKS) 1.2*max(YKS)]);

% xmax=max(XSAT(:,1+2*(i-1))
axis equal
auxmax=max([max(XKS/1000) max(YKS/1000)]);
auxmin=min([min(XKS/1000) min(YKS/1000)]);
%axis([1.2*auxmin 1.2*auxmax 1.2*auxmin 1.2*auxmax])
%legend('Model Trajectory','Non-Spoofed Trajectory','Authentic Trajectory','Spoofed Trajectory','Spoofed Satelite','Non-Spoofed Satelite')
legend('Model Trajectory','Non-Spoofed Trajectory','Authentic Trajectory','Spoofed Trajectory','GNSS only','GNSS only Spoofed')
title('X-Y (km)')
set(gca, 'FontName', 'Cambria','Fontsize',20)


figure('Name','Residues','NumberTitle','off','units','normalized','outerposition',[0 0 1 1],'Color',[234 253 234]/255)

subplot(2,1,1),
hold on
aux=0;aux(T)=0;taux=0;taux(T)=0;count=1;
for i=1:T
    if residuos(i) > 0
        aux(count)=(residuos(i));
        taux(count)=i*dt;
        count=count+1;
    end
end
aux=aux(1:count-1);taux=taux(1:count-1);
plot(taux,aux), grid minor, axis tight;
plot(t_begin_spoof*dt*ones(1,100),max(aux)*1.1*linspace(0,1,100),...
    t_end_spoof*dt*ones(1,100),max(aux)*1.1*linspace(0,1,100)), grid on, axis tight;
plot(SATPARAMS.false_alarm_start*dt*ones(1,100),max(aux)*1.1*linspace(0,1,100),...
    SATPARAMS.false_alarm_end*dt*ones(1,100),max(aux)*1.1*linspace(0,1,100)), grid on, axis tight;
legend('Authentic Residues','Begin Spoof','End Spoof','Begin FA','End FA');
xlabel('Time')
ylabel('Residues')
title('Non Spoofed Residues')
set(gca, 'FontName', 'Cambria','Fontsize',20)


subplot(2,1,2),
hold on
aux=0;aux(T)=0;taux=0;taux(T)=0;count=1;
for i=1:T
    if residuos_spoof(i) > 0
        aux(count)=log10(residuos_spoof(i));
        taux(count)=i*dt;
        count=count+1;
    end
end
aux=aux(1:count-1);taux=taux(1:count-1);
plot(taux,(aux)), grid minor, axis tight;
plot(t_begin_spoof*dt*ones(1,100),max(aux)*1.1*linspace( -1,1,100),...
    t_end_spoof*dt*ones(1,100),max(aux)*1.1*linspace(-1,1,100)), grid on, axis tight;
plot(SATPARAMS.false_alarm_start*dt*ones(1,100),max(aux)*1.1*linspace(-1,1,100),...
    SATPARAMS.false_alarm_end*dt*ones(1,100),max(aux)*1.1*linspace(-1,1,100)), grid on, axis tight;
legend('Spoofed Residues','Begin Spoof','End Spoof','Begin FA','End FA');
xlabel('Time')
ylabel('Residues (log10)')
title('Spoofed Residues')
set(gca, 'FontName', 'Cambria','Fontsize',20)


% figure(n_fig)
% close(n_fig)
% figure(n_fig)
figure('Name','Bias Comparison','NumberTitle','off','units','normalized','outerposition',[0 0 1 1],'Color',[234 253 234]/255)

subplot(3,2,1),
hold on
plot(t,biasy,t,biasky0+biasky1), grid minor, axis tight;
legend('Modelled Yb','Propagated Yb')
title('b_y:  Modelled vs Non-Spoofed')
xlabel('Time');ylabel('m/s^2')
set(gca, 'FontName', 'Cambria','Fontsize',20)

subplot(3,2,2),
hold on
aux=biasky0_s(floor(T_alineamiento*1.5):end)+biasky1_s(floor(T_alineamiento*1.5):end);
plot(t,biasky0+biasky1,t,biasky0_s+biasky1_s), grid on, axis([min(t) max(t) min(aux) max(aux)]);
%plot(t(1:end-1),diff(biasy),t,biasky1), grid minor, axis tight;
ax = gca;
aux1=[ax.XLim ax.YLim];
plot(t_begin_spoof*dt*ones(1,100),linspace(-1,1,100),t_end_spoof*dt*ones(1,100),linspace(-1,1,100)),axis(aux1);
legend('Authentic Yb','Spoofed Yb')
title('b_y')
xlabel('Time');ylabel('m/s^2')
set(gca, 'FontName', 'Cambria','Fontsize',20)

subplot(3,2,3),
hold on
plot(t,biasx,t,biaskx0+biaskx1), grid minor, axis tight;
legend('Modelled Xb','Propagated Xb')
title('b_x:  Modelled vs Non-Spoofed')
xlabel('Time');ylabel('m/s^2')
set(gca, 'FontName', 'Cambria','Fontsize',20)

subplot(3,2,4),
hold on
aux=biaskx0_s(floor(T_alineamiento*1.5):end)+biaskx1_s(floor(T_alineamiento*1.5):end);
plot(t,biaskx0+biaskx1,t,biaskx0_s+biaskx1_s), grid on, axis([min(t) max(t) min(aux) max(aux)]);
ax = gca;
aux1=[ax.XLim ax.YLim];
plot(t_begin_spoof*dt*ones(1,100),linspace(-1,1,100),t_end_spoof*dt*ones(1,100),linspace(-1,1,100)),axis(aux1);
legend('Authentic Xb','Spoofed Xb')
title('b_x')
xlabel('Time');ylabel('m/s^2')
set(gca, 'FontName', 'Cambria','Fontsize',20)

subplot(3,2,5),
hold on
plot(t,biasw,t,biaskw0+biaskw1), grid minor, axis tight;
legend('Modelled Wb','Propagated Wb')
title('b_w:  Modelled vs Non-Spoofed')
xlabel('Time');ylabel('m/s^2')
set(gca, 'FontName', 'Cambria','Fontsize',20)

subplot(3,2,6),
hold on
aux=biaskw0_s(floor(T_alineamiento*1.5):end)+biaskw1_s(floor(T_alineamiento*1.5):end);
plot(t,biaskw0+biaskw1,t,biaskw0_s+biaskw1_s), grid on;%, axis([min(t) max(t) min(aux) max(aux)]);
ax = gca;
aux1=[ax.XLim ax.YLim];
plot(t_begin_spoof*dt*ones(1,100),linspace(-1,1,100),t_end_spoof*dt*ones(1,100),linspace(-1,1,100)),axis(aux1);
legend('Authentic Wb','Spoofed Wb')
title('b_w')
xlabel('Time');ylabel('m/s^2')
set(gca, 'FontName', 'Cambria','Fontsize',20)

n_fig=n_fig+1;

figure('Name','Accel','NumberTitle','off','units','normalized','outerposition',[0 0 1 1],'Color',[234 253 234]/255)
subplot(2,1,1),
plot(t,dXK(3,:),t,AXR),grid minor, axis tight;
subplot(2,1,2),
plot(t,dXK(4,:),t,AYR),grid minor, axis tight;
set(gca, 'FontName', 'Cambria','Fontsize',20)


figure('Name','Satelites Integrity3','NumberTitle','off','units','normalized','outerposition',[0 0 1 1],'Color',[234 253 234]/255)

aux_matrix(T,n_satelites)=0;aux_matrix_s(T,n_satelites)=0;
t_aux(T,1)=0;t_aux_s(T,1)=0;
count=1;
for i=1:T
    if SIMULATION.IMU.checked_measures(i) > 0
        j=floor(SIMULATION.IMU.checked_measures(i)/10);
        if j > 0
            aux_matrix(count,j)=j;
            t_aux(count)=i*dt;
        end
        j=SIMULATION.IMU.checked_measures(i)-floor(SIMULATION.IMU.checked_measures(i)/10)*10;
        if j > 0
            aux_matrix(count,j)=j;
            t_aux(count)=i*dt;
        end
        count=count+1;
    elseif SIMULATION.IMU.checked_measures(i) == -1
        aux_matrix(count,:)=-1;
        t_aux(count)=i*dt;
        count=count+1;
    elseif SIMULATION.IMU.checked_measures(i) == 0
        aux_matrix(count,:)=0;
        t_aux(count)=i*dt;
        count=count+1;
    end
end
aux_matrix=aux_matrix(1:count-1,:);
t_aux=t_aux(1:count-1);

count=1;
for i=1:T
    if SIMULATION_s.IMU.checked_measures(i) > 0
        j=floor(SIMULATION_s.IMU.checked_measures(i)/10);
        if j > 0
            aux_matrix_s(count,j)=j;
            t_aux_s(count)=i*dt;
        end
        j=SIMULATION_s.IMU.checked_measures(i)-floor(SIMULATION_s.IMU.checked_measures(i)/10)*10;
        if j > 0
            aux_matrix_s(count,j)=j;
            t_aux_s(count)=i*dt;
        end
        count=count+1;
    elseif SIMULATION_s.IMU.checked_measures(i) == -1
        aux_matrix_s(count,:)=-1;
        t_aux_s(count)=i*dt;
        count=count+1;
        
    elseif SIMULATION_s.IMU.checked_measures(i) == 0
        aux_matrix_s(count,:)=0;
        t_aux_s(count)=i*dt;
        count=count+1;
    end
end
aux_matrix_s=aux_matrix_s(1:count-1,:);
t_aux_s=t_aux_s(1:count-1);

subplot(2,1,1),
hold on
plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(-10,10,100),'LineWidth',2);
plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(-10,10,100),'LineWidth',2);
plot(t_end_false_alarm*dt*ones(1,100),1.1*linspace(-10,10,100),'LineWidth',2);
for j=1:n_satelites
    scatter(t_aux,aux_matrix(:,j)),grid on,axis tight;
end
axis([min(t) max(t) -1 n_satelites])
legend('Spoofing Begin Time','FA Begin Time','FA End Time');
title('Non Spoofed Simulation Satelite Check')
xlabel('time');ylabel('index')
set(gca, 'FontName', 'Cambria','Fontsize',20)

subplot(2,1,2),
hold on
plot(t_begin_spoof*dt*ones(1,100),1.1*linspace(-10,10,100),'LineWidth',2);
plot(t_begin_false_alarm*dt*ones(1,100),1.1*linspace(-10,10,100),'LineWidth',2);
plot(t_end_false_alarm*dt*ones(1,100),1.1*linspace(-10,10,100),'LineWidth',2);
for j=1:n_satelites
    scatter(t_aux_s,aux_matrix_s(:,j)),grid on;
end
axis([min(t) max(t) -1 n_satelites])
legend('Spoofing Begin Time','FA Begin Time','FA End Time');
title('Spoofed Simulation Satelite Check')
xlabel('time');ylabel('index')
set(gca, 'FontName', 'Cambria','Fontsize',20)


n_fig=n_fig+1;

figure('Name','Distance Error','NumberTitle','off','units','normalized','outerposition',[0 0 1 1],'Color',[234 253 234]/255)

subplot(2,1,1),
hold on
aux=(XK'-XR)./XR;
aux2=(XKS'-XR)./XR;
taux=t;
plot(taux,aux), grid minor, axis tight;
plot(t_begin_spoof*dt*ones(1,100),max(aux)*1.1*linspace(-1,1,100),...
    t_end_spoof*dt*ones(1,100),max(aux)*1.1*linspace(-1,1,100)), grid on, axis tight;
plot(SATPARAMS.false_alarm_start*dt*ones(1,100),max(aux)*1.1*linspace(-1,1,100),...
    SATPARAMS.false_alarm_end*dt*ones(1,100),max(aux)*1.1*linspace(-1,1,100)), grid on, axis tight;
legend('Error','Begin Spoof','End Spoof','Begin FA','End FA');
xlabel('Time')
ylabel('Error')
title('X Error')
set(gca, 'FontName', 'Cambria','Fontsize',20)

subplot(2,1,2),
hold on
aux=(YK'-YR)./YR;
taux=t;
plot(taux,aux), grid minor, axis tight;
plot(t_begin_spoof*dt*ones(1,100),max(aux)*1.1*linspace(-1,1,100),...
    t_end_spoof*dt*ones(1,100),max(aux)*1.1*linspace(-1,1,100)), grid on, axis tight;
plot(SATPARAMS.false_alarm_start*dt*ones(1,100),max(aux)*1.1*linspace(-1,1,100),...
    SATPARAMS.false_alarm_end*dt*ones(1,100),max(aux)*1.1*linspace(-1,1,100)), grid on, axis tight;
legend('Error','Begin Spoof','End Spoof','Begin FA','End FA');
xlabel('Time')
ylabel('Error')
title('Y Error')
set(gca, 'FontName', 'Cambria','Fontsize',20)


end


function [prob_res,eje_x] = probability_distribution(data_input,time_begin,N_div)
%ORDENA LOS RESIDUOS EN DIFERENTES NIVELES DE CANTIDADES
aux=data_input/max(data_input);
aux=aux(floor(time_begin):end);
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