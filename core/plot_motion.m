function plot_motion(tiempo,WORLDMODEL,SATPARAMS,SIMULATION,SIMULATION_s)


%% -------------- CALCULOS PREVIOS  -----------------
XReal=WORLDMODEL.true_tray;

XR=XReal(:,1);YR=XReal(:,2);
dt=tiempo.dt;
dtGNSS=tiempo.dtGNSS;
T_alineamiento=tiempo.alineamiento;
T=size(XReal,1);t=0:dt:T*dt-dt;

XF=SIMULATION.state_vector;
dXK=SIMULATION.state_vector_derivate;
residuos=SIMULATION.residuos;

XFS=SIMULATION_s.state_vector;
dXKs=SIMULATION_s.state_vector_derivate;
residuos_spoof=SIMULATION_s.residuos;

radio=max(XR);
XK=XF(1,:);YK=XF(2,:);
biaskx0=XF(5,:);biasky0=XF(6,:);
biaskx1=XF(7,:);biasky1=XF(8,:);
THK=XF(9,:);

n_satelites=SATPARAMS.n_satelites;
n_satelites_spoof=SATPARAMS.n_satelites_spoof;
coeff_spoof0=SATPARAMS.coeff_spoof_0;
coeff_spoof1=SATPARAMS.coeff_spoof_1;
sigma_GNSS=SATPARAMS.sigma_GNSS;
t_begin_spoof=SATPARAMS.t_begin_spoof;
t_end_spoof=SATPARAMS.t_end_spoof;
XSAT=SATPARAMS.position;

radio=max(XR);
XK=XF(1,:);
YK=XF(2,:);
biaskx0=XF(5,:);
biasky0=XF(6,:);
biaskx1=XF(7,:);
biasky1=XF(8,:);

XKS=XFS(1,:);
YKS=XFS(2,:);
biaskx0_s=XFS(5,:);
biasky0_s=XFS(6,:);
biaskx1_s=XFS(7,:);
biasky1_s=XFS(8,:);


%APAÑO QUE HAGO PARA ORDENAR LOS RESIDUOS POR NIVELES

aux=residuos/max(residuos);
aux=aux(floor(T/4):end);
N_div=100;
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


aux=residuos_spoof/max(residuos_spoof);
aux=aux(floor(T/4):end);
prob_res_s(N_div)=0;div=linspace(1,0,N_div);
maxaux=max(aux);

for i=1:length(aux)
    for j=1:length(div)
        if aux(i) ~=0
            if aux(i)>(maxaux*div(j))
                prob_res_s(j)=prob_res_s(j)+1;
                break
            end
        end
    end
end
prob_res_s=(prob_res_s/length(aux));


%% -------  PLOTEO -------------
n_fig=1;


figure('units','normalized','outerposition',[0 0 1 1])
straux={'Real','Authentic','Spoofed'};
for j=1:n_satelites
    straux(j+3)={strcat('Satelite',string(j))};
end
limit=max(max(XReal))*2;
hold on
plot(XK,YK,'g-'), grid on, axis([2*min(XK)-0.2*max(XK) 2*max(XK) 1.2*min(YK)-0.2*max(YK) 2*max(YK)]);

%for i=floor(t_begin_spoof*0.6):floor(T/100):floor(t_end_spoof*1)
for i=1:floor(T/100):floor(t_end_spoof*1)
    % plot(XR(i),YR(i),XK(i),YK(i),XKS(i),YKS(i),':')
%     scatter(XR(i),YR(i),20,'MarkerEdgeColor',[0 .5 .5],...
%         'MarkerFaceColor',[0 .7 .7],...
%         'LineWidth',1.5);
    scatter(XK(i),YK(i),20,'MarkerEdgeColor',[0.5 .5 .0],...
        'MarkerFaceColor',[0.7 .7 .0],...
        'LineWidth',1.5);
    scatter(XKS(i),YKS(i),20,'MarkerEdgeColor',[0.5 .0 .5],...
        'MarkerFaceColor',[0.7 .0 .7],...
        'LineWidth',1.5);
    grid on, axis([2*min(XK)-0.2*max(XK) 2*max(XK) 1.2*min(YK)-0.2*max(YK) 2*max(YK)]);
    legend(straux)
    axis equal
    title(i/T*100)
%     for j=1:n_satelites
%         if j<= n_satelites_spoof
%             colorp=[153, 0, 0]/255;
%         else
%             colorp=[0, 153, 76]/255;
%         end
%         scatter(XSAT(i,1+2*(j-1)),XSAT(i,2+2*(j-1)),20,'MarkerEdgeColor',[0 .5 .5],...
%             'MarkerFaceColor',colorp,...
%             'LineWidth',0.5);
%     end
    pause(0.001)
end


end

