%LIF model taken from Farokhniaee 2019

% % transmission + synaptic delay: td 
% td=2; %2 ms for trasmission and .5 ms for synaptic delay
% 
% dt=.1; ti=dt; tf=1000+td;%tf=1500+td;%tf=61000+td; %in mili seconds 
% t=ti:dt:tf;
% 
% %noise term
% % wght=0;   %no noise
% wght=.5; %default noise
% % wght=5;   %high noise
% kisi=wght*randn(1,length(t));
% 
% %input spike train
% sp=zeros(1,length(t));
% 
% % Neuron parameters: (for ~20 Hz base firing .56 and for ~8-10 Hz choose .26)
% Cm= 1; Rm=100; Ie=.56; %(for deterministic model)
% % Ie=.16; %subthreshold firing (for noise purpose, stochastic model)
% El=-70; Vth=-54; 
% Vreset=-80;
% 
% % % Neuron parameters: (for 62.5 Hz base firing)
% % Cm= 1; Rm=100; Ie=1.52; %(for deterministic model)
% % % Ie=.18; %subthreshold firing (for noise purpose, stochastic model)
% % El=-70; Vth=-54; 
% % Vreset=-80;
% 
% % Compute neuron firing pattern with and without synaptic input:
% Vin=zeros(1,length(t));
% 
% wk=10; %Poissonian weight
% poiss=wk*rand(1,length(sp)).*sp(1,:);
% for i=1:length(t)-1
% Vin(i+1) = Vin(i) + (dt/Cm)*(((El-Vin(i))/Rm) + Ie + poiss(i) + kisi(i));
% if Vin(i+1)>= Vth+kisi(i)
%     Vin(i)=20+kisi(i);
%     Vin(i+1)=Vreset+kisi(i);
% end
% end
% 
% plot(Vin)
%% Implementation for our code
clear all
close all

Cm= 1; Rm=100; Ie=.56;
El=-62+randn(1,1)*5; Vth=-54; 
Vreset=-62+randn(1,1)*5;

%time variables
tmax=1000; %maximum time (ms)
dt=0.01; %timestep (ms)
t=0:dt:tmax; %time vector
n=10; %number of neurons in each nucleus (TH, STN, GPe, GPi)

v2=-62+randn(n,1)*5; %for STN
vsn=zeros(n,length(t)); %STN membrane voltage
vsn(:,1)=v2; %for STN
sp=zeros(1,length(t)); %check out

%noise term
wght=.5; %default noise
kisi=wght*randn(n,length(t));


for i=1:length(t)-1  

    if i==1
        V2=vsn(:,i);
    else
        V2=vsn(:,i-1);
    end
    

    %STN
    vsn(:,i+1) = vsn(:,i) + (dt/Cm)*(((El-V2)/Rm) + Ie);
    if vsn(:,i+1)>= Vth
        vsn(:,i)=40;
        vsn(:,i+1)=Vreset;
    end
    
    
end

figure
plot(vsn(1,:));