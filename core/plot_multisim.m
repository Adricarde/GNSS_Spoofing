function plot_multisim(WORLDMODEL,dataoutput,DataFolder,savestring,time_params)
%% -------------- CALCULOS PREVIOS  -----------------
XReal=WORLDMODEL.true_tray;
XR=XReal(:,1);YR=XReal(:,2);

guess_vect=dataoutput.guess_vect;
trajectories=dataoutput.trajectories;
indicator=dataoutput.indicator;
detect_time=dataoutput.detect_time;
spoof_time=dataoutput.spoof_time;
deviated_distance=dataoutput.deviated_distance;
false_alarm_time=dataoutput.false_alarm_time;
dtGNSSvect=dataoutput.dtGNSSvect;
simulation_time=size(XReal,1)*time_params.dt;
endtime=dataoutput.endtime;
n_simulations=size(trajectories,1);

%% -------  PLOTEO -------------

figure('Name','Trajectories','units','normalized','outerposition',[0 0 1 1])
hold on
plot(XReal(:,1)/1000,XReal(:,2)/1000,'k-');
for i=1:n_simulations
    XKS(:,1,1)=trajectories(i,1,1:endtime(i));YKS(:,1,1)=trajectories(i,2,1:endtime(i));
    if max(abs(XKS)) < 500*max(abs(XReal(:,1))) && max(abs(YKS)) < 500*max(abs(XReal(:,2)))
        plot(XKS/1000,YKS/1000,'-.','LineWidth',3,'Color',[i/n_simulations 0.2 0.35]), grid on;%
        legstring(i+1)='Simulation: '+string(i);
    end
    clear XKS YKS
end
legstring(1)='Real Trajectory';
axis equal
% legend(legstring)
title('X-Y (km)')
set(gca, 'FontName', 'Cambria','Fontsize',20)
saveas(gcf,strcat(DataFolder,savestring,'trajectories.png'))

%----------------------------

figure('Name','Detect Time','units','normalized','outerposition',[0 0 1 1])
hold on
for i=1:7
    aux=gca;
    plot(guess_vect(:,i),'--o','MarkerFaceColor',aux.ColorOrder(i,:),'Color',aux.ColorOrder(i,:)),axis tight, grid on;
end
area(spoof_time,'FaceAlpha',0.2,'EdgeAlpha',0.2,'EdgeColor','red'),axis tight, grid on;
area(false_alarm_time,'FaceAlpha',0.2,'EdgeAlpha',0.2,'EdgeColor','black'),axis tight, grid on;
plot(size(XReal,1)*time_params.dt*ones(1,n_simulations),':b'),axis tight, grid on;
xlabel('Nº de Simulación');
ylabel('Tiempo de deteccion de Spoofing');
legend('Bias Drift','Residues','Double Prop','Chi','Distance','Ellipses','Discard','Init Spoofing','Init False Alarm')
set(gca, 'FontName', 'Cambria','Fontsize',20)
saveas(gcf,strcat(DataFolder,savestring,'detect_time.png'))

%----------------------------
figure('Name','Indicator','units','normalized','outerposition',[0 0 1 1])
hold on
bar(indicator,'stacked'),axis tight, grid on;
xlabel('Nº de Simulación');
ylabel('Tiempo de deteccion de Spoofing');
legend('Bias Drift','Residues','Double Prop','Chi','Distance','Ellipses','Discard')
set(gca, 'FontName', 'Cambria','Fontsize',20)
saveas(gcf,strcat(DataFolder,savestring,'indicator.png'))

%----------------------------
figure('Name','Total Time','units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1),
hold on
for i=1:7
    scatter(linspace(1,n_simulations,n_simulations),guess_vect(:,i),'filled'),axis tight, grid on;
end
plot(spoof_time,'--r'),axis tight, grid on;
plot(false_alarm_time,'--k'),axis tight, grid on;
plot(size(XReal,1)*time_params.dt*ones(1,n_simulations),':b'),axis tight, grid on;
xlabel('Nº de Simulación');
ylabel('Tiempo');
legend('Bias Drift','Residues','Double Prop','Chi','Distance','Ellipses','Discard','Init Spoofing','Init False Alarm')
set(gca, 'FontName', 'Cambria','Fontsize',20)
subplot(2,1,2),
hold on
bar(deviated_distance),axis tight, grid on;
xlabel('Nº de Simulación');
ylabel('Distancia Desviada');
set(gca, 'FontName', 'Cambria','Fontsize',20)
saveas(gcf,strcat(DataFolder,savestring,'total_time.png'))

%----------------------------

figure('Name','Histogram','units','normalized','outerposition',[0 0 1 1])
hold on
area(indicator),axis tight, grid on;
xlabel('Nº de Simulación');
ylabel('Tiempo de deteccion de Spoofing');
legend('Bias Drift','Residues','Double Prop','Chi','Distance','Ellipses','Discard')
set(gca, 'FontName', 'Cambria','Fontsize',20)
saveas(gcf,strcat(DataFolder,savestring,'histogram.png'))

%----------------------------

for i=1:n_simulations
   total_time=sum(indicator(i,:));
   percent_indicator(i,:)= (indicator(i,:))/total_time*100;
end
figure('Name','Histogram Percentage','units','normalized','outerposition',[0 0 1 1])
hold on
area(percent_indicator),axis tight, grid on;
xlabel('Nº de Simulación');
ylabel('Tiempo de deteccion de Spoofing');
legend('Bias Drift','Residues','Double Prop','Chi','Distance','Ellipses','Discard')
set(gca, 'FontName', 'Cambria','Fontsize',20)
saveas(gcf,strcat(DataFolder,savestring,'histogram_percent.png'))
end