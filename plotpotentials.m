%plotting utility for potential in each part of the network
figure;
subplot(2,5,1); 
%ax=plotyy(t,vppn(1,:),t,Istim(1:tmax/dt+1)); %plot for 1st neuron both Istim and Vppn
%set(ax(1),'XLim',[0 tmax],'YLim',[-100 80],'Visible','on')
plot(t,vppn(1,:));axis([0 tmax -100 80 ])
%set(ax(2),'XLim',[0 tmax],'YLim',[-2 30],'Visible','off')
title(['PPN, FR: ' num2str(int64(fr4)) 'Hz'])
ylabel('Vm (mV)'); xlabel('Time (msec)');

subplot(2,5,2); 
ax=plotyy(t,vth(1,:),t,Istim(1:tmax/dt+1)); %plot for 1st neuron both Istim and Vth
set(ax(1),'XLim',[0 tmax],'YLim',[-100 20],'Visible','on')
set(ax(2),'XLim',[0 tmax],'YLim',[-2 30],'Visible','off')
title(['Thalamus, FR: ' num2str(int64(fr5)) 'Hz'])
ylabel('Vm (mV)'); xlabel('Time (msec)');

subplot(2,5,3); %for 1st STN neuron
plot(t,vsn(1,:));axis([0 tmax -100 80 ]) 
title(['STN, FR: ' num2str(int64(fr1)) 'Hz']);
ylabel('Vm (mV)'); xlabel('Time (msec)');

subplot(2,5,4); %for 1st GPe neuron
plot(t,vge(1,:));axis([0 tmax -100 80 ]) 
title(['GPe, FR: ' num2str(int64(fr2)) 'Hz'])
ylabel('Vm (mV)'); xlabel('Time (msec)');

subplot(2,5,5); %for 1st GPi neuron
plot(t,vgi(1,:));axis([0 tmax -100 80 ]) 
title(['GPi, FR: ' num2str(int64(fr3)) 'Hz'])
ylabel('Vm (mV)'); xlabel('Time (msec)');

subplot(2,5,6); %for 1st striatum neuron
plot(t,vstr(1,:));axis([0 tmax -100 80 ]) 
title(['Striatum, FR: ' num2str(int64(fr6)) 'Hz'])
ylabel('Vm (mV)'); xlabel('Time (msec)');

subplot(2,5,7); %for 1st SNr neuron
plot(t,vsnr(1,:));axis([0 tmax -100 80 ]) 
title(['SNr, FR: ' num2str(int64(fr7)) 'Hz'])
ylabel('Vm (mV)'); xlabel('Time (msec)');

subplot(2,5,8); %for 1st PRF neuron
plot(t,vprf(1,:));axis([0 tmax -100 80 ]) 
title(['PRF, FR: ' num2str(int64(fr8)) 'Hz'])
ylabel('Vm (mV)'); xlabel('Time (msec)');

subplot(2,5,9); %for 1st CNF neuron
plot(t,vcnf(1,:));axis([0 tmax -100 80 ]) 
title(['CNF, FR: ' num2str(int64(fr9)) 'Hz'])
ylabel('Vm (mV)'); xlabel('Time (msec)');

subplot(2,5,10); %for 1st LC neuron
plot(t,vlc(1,:));axis([0 tmax -100 80 ]) 
title(['LC, FR: ' num2str(int64(fr10)) 'Hz'])
ylabel('Vm (mV)'); xlabel('Time (msec)');

suptitle({['Firing patterns in freezing of gait network']
           ['Thalamic relay EI: ' num2str(titleval)]});