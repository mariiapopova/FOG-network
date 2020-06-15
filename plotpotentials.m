%plotting utility for potential in each part of the network

% figure;
% subplot(2,2,1); 
% ax=plotyy(t,vth(1,:),t,Istim(1:tmax/dt+1)); %plot for 1st neuron both Istim and Vth
% set(ax(1),'XLim',[0 tmax],'YLim',[-100 20],'Visible','on')
% set(ax(2),'XLim',[0 tmax],'YLim',[-2 30],'Visible','off')
% title('Thalamus')
% ylabel('Vm (mV)'); xlabel('Time (msec)');
% 
% subplot(2,2,2); %for 1st STN neuron
% plot(t,vsn(1,:));axis([0 tmax -100 80 ]) 
% title('STN');
% ylabel('Vm (mV)'); xlabel('Time (msec)');
% 
% subplot(2,2,3); %for 1st GPe neuron
% plot(t,vge(1,:));axis([0 tmax -100 80 ]) 
% title('GPe')
% ylabel('Vm (mV)'); xlabel('Time (msec)');
% 
% subplot(2,2,4); %for 1st GPi neuron
% plot(t,vgi(1,:));axis([0 tmax -100 80 ]) 
% title('GPi')
% ylabel('Vm (mV)'); xlabel('Time (msec)');

figure;
subplot(2,2,1); 
ax=plotyy(t,vppn(1,:),t,Istim(1:tmax/dt+1)); %plot for 1st neuron both Istim and Vppn
set(ax(1),'XLim',[0 tmax],'YLim',[-100 80],'Visible','on')
set(ax(2),'XLim',[0 tmax],'YLim',[-2 30],'Visible','off')
title('PPN')
ylabel('Vm (mV)'); xlabel('Time (msec)');

subplot(2,2,2); %for 1st STN neuron
plot(t,vsn(1,:));axis([0 tmax -100 80 ]) 
title('STN');
ylabel('Vm (mV)'); xlabel('Time (msec)');

subplot(2,2,3); %for 1st GPe neuron
plot(t,vge(1,:));axis([0 tmax -100 80 ]) 
title('GPe')
ylabel('Vm (mV)'); xlabel('Time (msec)');

subplot(2,2,4); %for 1st GPi neuron
plot(t,vgi(1,:));axis([0 tmax -100 80 ]) 
title('GPi')
ylabel('Vm (mV)'); xlabel('Time (msec)');
