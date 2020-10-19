clear all
close all
%file for lc check 
Csnr=1;% muF

gna=[120]; Ena=[55];
gk=[20]; Ek=[-72];
gl=[0.3]; El=[-17];
ga=[47.7];
B=0.21.*(ga/gk);


%time variables
tmax=1000; %maximum time (ms)
dt=0.01; %timestep (ms)
t=0:dt:tmax; %time vector
n=10; %number of neurons in each nucleus (TH, STN, GPe, GPi)

v10=-62+randn(n,1)*5; %for LC
vlc=zeros(n,length(t)); %LC membrane voltage
vlc(:,1)=v10; %for LC

q=lc_qinf(vlc(:,1));%q for lc 

for i=2:length(t)  
    V10=vlc(:,i-1);
    
    qinf=lc_qinf(V10); m10=lc_minf(V10); b10=lc_binf(V10);
    tq=lc_tauq(V10); 
    
    %LC cell currents
    Ina10=gna(1)*(m10.^3).*(-3.*(q-B.*b10)+0.85).*(V10-Ena(1));
    Ik10=gk(1)*q.*(V10-Ek(1));
    Il10=gl(1).*(V10-El(1));
    Ib=5;
    
    %LC
    vlc(:,i)=V10+dt*(1/Csnr*(Ib-Ina10-Ik10-Il10));
    q=q+dt*((qinf-q)./tq); 
end

figure
plot(vlc(1,:));

%fr=findfreq(vlc(1,:));
%disp(fr);