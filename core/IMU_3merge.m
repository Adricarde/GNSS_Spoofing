function [IMU_SIMULATION] = IMU_3merge(WORLDMODEL,IMU_params,dt)


ABody=WORLDMODEL.forces;

[IMUSIM1] = IMU_sim(WORLDMODEL,IMU_params,dt);
[IMUSIM2] = IMU_sim(WORLDMODEL,IMU_params,dt);
[IMUSIM3] = IMU_sim(WORLDMODEL,IMU_params,dt);
AM1=IMUSIM1.measures;
AM2=IMUSIM2.measures;
AM3=IMUSIM3.measures;
d1=AM1-AM2;
d2=AM1-AM3;
d3=AM2-AM3;
mean1=(AM1+AM2)/2;
mean2=(AM1+AM3)/2;
mean3=(AM2+AM3)/2;
mean_total=(AM1+AM2+AM3)/3;

threshold=IMUSIM1.params.sigma_bias_a*5;
AM_merge=AM1;T=size(AM1,1);
i1=1;i2=1;i3=1;
for i=1:size(AM1,1)
    if abs(d1(i)) > threshold && abs(d2(i)) > threshold
        AM_merge(i,1:3)=mean3(i,1:3);
        i1=i1+1;
    elseif abs(d2(i)) > threshold && abs(d3(i)) > threshold
        AM_merge(i,1:3)=mean1(i,1:3);
        i2=i2+1;
    elseif abs(d1(i)) > threshold && abs(d3(i)) > threshold
        AM_merge(i,1:3)=mean2(i,1:3);
        i3=i3+1;
    else
        AM_merge(i,1:3)=mean_total(i,1:3);
    end
end

biasx=AM_merge(:,1)-ABody(:,1);
biasy=AM_merge(:,2)-ABody(:,2);
biasw=AM_merge(:,3)-ABody(:,3);
IMU_SIMULATION=IMUSIM1;
IMU_SIMULATION.measures=AM_merge;
IMU_SIMULATION.biasx=biasx;
IMU_SIMULATION.biasy=biasy;
IMU_SIMULATION.biasw=biasw;

if IMU_params.plot
    close 300
    figure('Name','3 IMU Merge','NumberTitle','off','units','normalized','outerposition',[0 0 1 1],'Color',[234 253 234]/255)
    subplot(6,1,1),
    hold on
    plot(linspace(1,T,T)*dt,biasy),grid on, axis tight;
    
    subplot(6,1,2),
    hold on
    plot(linspace(1,T,T)*dt,AM1(:,2)),grid on, axis tight;
    plot(linspace(1,T,T)*dt,AM2(:,2)),grid on, axis tight;
    plot(linspace(1,T,T)*dt,AM3(:,2)),grid on, axis tight;
    plot(linspace(1,T,T)*dt,AM_merge(:,2)),grid on, axis tight;
    
    subplot(6,1,3),
    hold on
    plot(linspace(1,T,T)*dt,biasx),grid on, axis tight;
    
    subplot(6,1,4),
    hold on
    plot(linspace(1,T,T)*dt,AM1(:,1)),grid on, axis tight;
    plot(linspace(1,T,T)*dt,AM2(:,1)),grid on, axis tight;
    plot(linspace(1,T,T)*dt,AM3(:,1)),grid on, axis tight;
    plot(linspace(1,T,T)*dt,AM_merge(:,1)),grid on, axis tight;
    
    subplot(6,1,5),
    hold on
    plot(linspace(1,T,T)*dt,biasw),grid on, axis tight;
    
    subplot(6,1,6),
    hold on
    plot(linspace(1,T,T)*dt,AM1(:,3)),grid on, axis tight;
    plot(linspace(1,T,T)*dt,AM2(:,3)),grid on, axis tight;
    plot(linspace(1,T,T)*dt,AM3(:,3)),grid on, axis tight;
    plot(linspace(1,T,T)*dt,AM_merge(:,3)),grid on, axis tight;
    pause(0.01)
end
end


