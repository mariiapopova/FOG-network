clear all
close all
addpath('gating')
%file for snr check 
Csnr=1;% pF

gl=[0.05 2.25 0.1 0.04]; El=[-70 -60 -65 -60];
gna=[3 37 120 30 35]; Ena=[50 55 55 45 50]; 
gk=[5 45 30 3.2 50]; Ek=[-75 -80 -80 -95 -90];
Casnrout=4; tmca6=0.5; thca6=18; %params for calcium currents for snr
gca=[0 2 0.15 0.7]; Eca=[0 140 120];
gnap6=0.175;
gsk6=0.4; %from fujita et al. maybe check it out more!
ksk=0.4; nsk=4;
alphaca=1*10^(-8); Casnrmin=5*10^(-8); tca=250/0.01;%in ms? %params for calcium currents for snr

%time variables
tmax=1000; %maximum time (ms)
dt=0.01; %timestep (ms)
t=0:dt:tmax; %time vector
n=10; %number of neurons in each nucleus (TH, STN, GPe, GPi)

v6=-62+randn(n,1)*5; %for SNr
vsnr=zeros(n,length(t)); %SNr membrane voltage
vsnr(:,1)=v6; %for SNr

M6=snr_minf(vsnr(:,1));%m for snr na
H6=snr_hinf(vsnr(:,1));%h for snr na
Sna6=snr_snainf(vsnr(:,1));%s for snr for na
Mnap6=snr_mnapinf(vsnr(:,1));%m for nap for snr
Hnap6=snr_hnapinf(vsnr(:,1));%h for nap for snr
Mk6=snr_mkinf(vsnr(:,1));%m for snr na
Hk6=snr_hkinf(vsnr(:,1));%h for snr na
Casksnr=Casnrmin; %for sure? check this!
Mca6=snr_mcainf(vsnr(:,1));% m for ca for snr
Hca6=snr_hcainf(vsnr(:,1));% h for ca for snr

for i=2:length(t)  
    V6=vsnr(:,i-1);
    
    m6=snr_minf(V6); h6=snr_hinf(V6); sna6=snr_snainf(V6); mk6=snr_mkinf(V6); hk6=snr_hkinf(V6);
    mnap6=snr_mnapinf(V6); hnap6=snr_hnapinf(V6); mca6=snr_mcainf(V6); hca6=snr_hcainf(V6);
    
    tm6=snr_taum(V6); th6=snr_tauh(V6); tsna6=snr_tausna(V6); tmk6=snr_taumk(V6); thk6=snr_tauhk(V6);
    tmnap6=snr_taumnap(V6); thnap6=snr_tauhnap(V6);
    
    %SNr cell currents
    Il6=gl(4)*(V6-El(4));
    Iappsnr=0.77; %healthy
    %Iappsnr=0.63; %parkinsonian
    Ina6=gna(5)*(M6.^3).*H6.*Sna6.*(V6-Ena(5));
    Msk=1./(1+((ksk./Casksnr).^nsk));
    Isk6= gsk6*Msk.*(V6-Ek(5));
    Inap6=gnap6*(Mnap6.^3).*Hnap6.*(V6-Ena(5));
    Ik6=gk(5)*(Mk6.^4).*Hk6.*(V6-Ek(5));
    Ecasnr=13.27*log(Casnrout./Casksnr);
    Ica6=gca(4)*Mca6.*Hca6.*(V6-Ecasnr);
    
    %SNr
    vsnr(:,i)=V6+dt*(1/Csnr*(-Ina6-Inap6-Ik6-Ica6-Isk6-Il6+Iappsnr));
    M6=M6+dt*((m6-M6)./tm6); 
    H6=H6+dt*((h6-H6)./th6); 
    Sna6=Sna6+dt*((sna6-Sna6)./tsna6); 
    Mnap6=Mnap6+dt*((mnap6-Mnap6)./tmnap6);
    Hnap6=Hnap6+dt*((hnap6-Hnap6)./thnap6);
    Mk6=Mk6+dt*((mk6-Mk6)./tmk6); 
    Hk6=Hk6+dt*((hk6-Hk6)./thk6); 
    Casksnr=Casksnr+dt*(-alphaca.*Ica6-((Casksnr-Casnrmin)./tca));
    Mca6=Mca6+dt*((mca6-Mca6)./tmca6); 
    Hca6=Hca6+dt*((hca6-Hca6)./thca6); 
end

figure
plot(vsnr(1,:));

fr=findfreq(vsnr(1,:));
disp(fr);